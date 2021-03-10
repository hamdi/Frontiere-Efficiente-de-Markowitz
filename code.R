packages <- c("nortest","urca","FactoMineR","factoextra","pairsD3","dHSIC","moments","partitions","zoo","rugarch","rmgarch",
              "PerformanceAnalytics","dplyr","tidyr","tseries","FinTS","FRAPO","fPortfolio")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(tidyr)
library(nortest)
library(urca)
library(tseries)
library(FactoMineR)
library(factoextra)
library(pairsD3)
library(dHSIC)
library(moments)
library(partitions)
library(zoo)
library(rugarch)
library(rmgarch)
library(PerformanceAnalytics)
library(FRAPO)
library(fPortfolio)
library(lattice)
library(FinTS)


# Chargement des données
df <- read.csv("donnees.csv", header = TRUE)
df_test <- df[209:312,]    # Données de test (2016 et 2017)
df <- df[1:208,]           # Données de construction du portefeuille (de 2012 à 2015)
n <- nrow(df)-1
p <- ncol(df)
nt <- nrow(df_test)-1
# Ajout d'un actif sans risque de rendement annuel 2.458% (égal à une OAT française de 2014 de maturité 10 ans)
df["AAA"] <- exp(0.0004670*seq(0,n))    # ln(1.02458)/52 = 0.0004670
df_test["AAA"] <- exp(0.0004670*seq(n+1,nt+n+1))

# Calcul des rendements logarithmiques
r <- log(df[2:(n+1),2:(p+1)]/df[1:n,2:(p+1)])  # données de construction du portefeuille
rt <- log(df_test[2:(nt+1),2:(p+1)]/df_test[1:nt,2:(p+1)])  # données de test

rownames(r) <- NULL
rownames(rt) <- NULL

plot(df[,"ACA"],type="l", main="Cours d'une action", ylab = "ACA")
plot(r[,"ACA"],type="l", main="Rendements hebdomadaires d'une action", ylab = "ACA")
# Calcul des volatilités simples
vol <- apply(r, 2, sd)
matrice_variance = var(r)


##### Normalité #####
# Test de Shapiro-Wilk
normalite <- apply(r, 2, shapiro.test) %>% lapply("[[", 2)
# Test de Lilliefors
normalite_2 <- apply(r, 2, lillie.test) %>% lapply("[[", 2)
# Test de Anderson-Darling
normalite_3 <- apply(r, 2, ad.test) %>% lapply("[[", 2)


print("Les rendements de ces actions ne suivent pas une loi normale selon l'un des tests:")
for (i in 1:p){
  if (min(unlist(c(normalite[i], normalite_2[i], normalite_3[i])))< 0.05) {
    print(names(normalite)[i])
  }
}

##### Stationnarité #####
# Test KPSS
stationnarite <- apply(r, 2, kpss.test) %>% lapply("[[", 3)
# Tous les rendements sont stationnaires (p-val > 0.05)


# Test ADF
tests=c()
for (i in 1:41) {
  adf=ur.df(r[,i],type="drift", selectlags = "AIC")
  tests=c(tests,adf@teststat[1])
}
names(tests) <- names(r)
barplot(tests, horiz=TRUE)
abline(v=-3.46) # valeur critique au niveau 1%
# Pour toutes les séries, la statistique du test est supérieure en module à la valeur critique
# Toutes les séries sont donc stationnaires


#### Sélection des actifs ####
# Calcul du ratio de Sharpe avec volatilité simple :
r_moy <- apply(r, 2, mean)
sharpe <- sqrt(52)*(r_moy-r_moy[p])/vol
plot(sharpe)
# La plupart des ratios de Sharpe sont entre 0 et 1.

# Matrice de corrélations
source("http://www.sthda.com/upload/rquery_cormat.r")
mat_corr <- rquery.cormat(r[,-41])   # on remarque que les actifs sont significativement corrélés au seuil de 0,1%

# ACP
res_acp <- PCA(r[,-41], graph=FALSE)
fviz_eig(res_acp, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(res_acp, col.var = "black") 
# Au regard du cercle de corrélation, on peut dire que la plupart des actifs sont fortement corrélés
# Cela est probablement dû à leur appartenance au même indice CAC 40, donc les traders achètent ou vendent ces actifs simultanément
# en tradant les trackers du CAC 40.

# Cos2 total des variables sur Dim.1 et Dim.2
fviz_cos2(res_acp, choice = "var", axes = 1:2) 
#les variables sont globalement bien representées sur les deux premiers axes : somme des cos2 assez élevée en moyenne

# Actifs retenus :
actifs=names(r)
retenus <- actifs[order(sharpe, decreasing=T)] # on retient le top 20 ratio de Sharpe
retenus <- retenus[1:20]


# Matrice de corrélation des actifs retenus
mat_corr <- rquery.cormat(r[retenus])
corr <- rquery.cormat(r[retenus], type="full", graph=FALSE)
corr <- apply(corr$r, 1, sum)
retenus <- names(corr[order(corr, decreasing=T)])[6:20]  # on garde les 15 les moins corrélés
rf <- r[1,p]
r_all <- r
r <- r[,retenus]
rt_all <- rt
rt <- rt[,retenus]
p <- length(retenus)

# ACP des actifs retenus
res_acp <- PCA(r, graph=FALSE)
fviz_pca_var(res_acp, col.var = "black") 
# Il reste 2 clusters corrélés: (OR,SAF,EL) mais ils sont répartis dans des industries différentes (cosmétique, aéronautique, luxe)
# et (RNO,VIE,CAP,DG,ML) qui représentent les industries automobile, énergétique, numérique, travaux, et pneumatique.
# Les actions gardées sont donc assez bien diversifiées


##  Indépendance des rendements
# Visualiser l'(in)dépendance entre les actifs avec un scatter matrix
pairsD3(r, big = TRUE, opacity = 0.9, cex = 1, width = 1200)
# Les rendements de la plupart des actifs sont assez corrélés -> non indépendants
dhsic.test(r)
# Le test dHSIC d'indépendance jointe des rendement rejette l'hypothése nulle d'indépendace avec une p-valeur < 0.01

# Skewness et kurtosis
skew <- apply(r, 2, skewness)
kurt <- apply(r, 2, kurtosis)
plot(skew, main="Skewness", pch=19)
plot(kurt, main="Kurtosis", pch=19)
mean(skew)
mean(kurt)
jb <- apply(r, 2, jarque.bera.test) %>% lapply("[[", 3)
print("Les rendements de ces actions ne suivent pas une loi normale selon Jarque-Bera:")
for (i in 1:p){
  if (jb[i]< 0.05) {
    print(retenus[i])
  }
}


# Etude de l'ensemble des rendements de toutes les actions
all_returns <- as.vector(unlist(r))
h <- hist(all_returns, breaks=20, prob=T, main="Histogramme des rendements", xlab="Rendements")
xfit <- seq(min(all_returns), max(all_returns), length = 80) 
yfit <- dnorm(xfit, mean = mean(all_returns), sd = sd(all_returns)) 
#yfit <- yfit * diff(h$mids[1:2]) * length(all_returns) 
lines(xfit, yfit, col = "cyan2", lwd = 2)
lines(density(all_returns), lty="dotted", col="coral1", lwd=2)
legend(x=0.1, y=11.5, legend=c("Densité normale", "Densité empirique"),
       col=c("cyan2", "coral1"), lty=1:2, cex=0.8)

skewness(all_returns)
# skewness proche de 0, donc les rendements sont symétriques
kurtosis(all_returns)
# kurtosis > 3 -> (leptokurtique)
qqnorm(all_returns)
qqline(all_returns)
# Les quantiles extremes sont significativement plus élevés que les quantiles de la loi normale
# La distribution des rendements contient donc beaucoup de valeurs extremes
# Pour mieux rendre compte de ce risque, il faudra utiliser des indicateurs de queue de distribution
# tels que la VaR et la CVaR

# Hétéroscédasticité
arch <- c()
for (i in 1:15){
  arch <- ArchTest(r[,i])$p.value
  if (arch>0.05){print(i)}
}
# Le test ne détecte pas d'homoscédasticité dans la plupart des actions


## Simulation des rendements et volatilités sur l'espace des poids possibles
# Attention : Prend environ 5 min pour s'exécuter
# Sans actif sans risque:
# w<-t(as.matrix(compositions(10,15))/10)
# cov_sim <- var(r)
# sim_r <- (as.matrix(w) %*% as.vector(colMeans(r)))
# sim_s <- apply(w, 1, function(x){sqrt(t(x) %*% cov_sim %*% x)})
# ggplot(data.frame(sim_s,sim_r), aes(sim_s, sim_r, color = sim_s, alpha = 0.4)) +
#   geom_point(shape = 16, size = 1, show.legend = FALSE) +
#   theme_light() +
#   scale_color_gradient(low = "#0091ff", high = "#f0650e") +
#   ggtitle("Rendements et volatilités empiriques pour différentes pondérations des actifs") +
#   xlab("Volatilité") + ylab("Espérance de rendement")

# Avec actif sans risque
# w <- t(as.matrix(compositions(10,16))/10)
# r_rf <- cbind(r, AAA=rep(rf,n))
# cov_sim <- var(r_rf)
# sim_r <- as.matrix(w) %*% as.vector(colMeans(r_rf))
# sim_s <- apply(w, 1, function(x){sqrt(t(x) %*% cov_sim %*% x)})
# ggplot(data.frame(sim_s,sim_r), aes(sim_s, sim_r, color = sim_s, alpha = 0.4)) +
#   geom_point(shape = 16, size = 1, show.legend = FALSE) +
#   theme_light() +
#   scale_color_gradient(low = "#0091ff", high = "#f0650e") +
#   ggtitle("Rendements et volatilités empiriques avec actif sans risque") +
#   xlab("Volatilité") + ylab("Espérance de rendement")


## Traçage de la frontière efficiente théorique
source("http://freakonometrics.free.fr/portfolio.r")
er <- apply(r, 2, mean)
matrice_variance <- var(r)
ef <- efficient.frontier(er, matrice_variance, alpha.min=-2, alpha.max=1.5, nport=20)
tp <- tangency.portfolio(er, matrice_variance,rf)
plot(ef, plot.assets=T)
abline(a=rf, b=((tp$er-rf)/tp$sd), lwd=2)
abline(v=0, col="lightgrey")
abline(h=0, col="lightgrey")
text(0.002, rf, labels="A.S.R.")



# Stabilité des portefeuilles efficients
# On calcule le portefeuille tangent sur les données de 2 ans et on le compare avec celui
# obtenu avec les données de 4 ans
r2 <- r[(n %/% 2):n,]
er1 <- apply(r, 2, mean)
matrice_variance1 <- var(r)
er2 <- apply(r2, 2, mean)
matrice_variance2 <- var(r2)
tp1 <- tangency.portfolio(er1, matrice_variance1,rf)$weights
tp2 <- tangency.portfolio(er2, matrice_variance2,rf)$weights
x <- seq(1,15)
d <- data.frame(tp1, tp2, x)
d <- pivot_longer(d, cols=c('tp1', 'tp2'), names_to='Periode', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Periode)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Pondérations du portefeuille tangent en fonction de la période de calcul") +
  scale_fill_discrete(name = "Période", labels = c("4 ans", "2 ans"))
# Certains poids restent stables mais la plupart sont modifiés



### VaR et CVaR ###
var <- c()
for (i in seq(1,p)){
  var <- c(var,VaR(as.vector(unlist(r[i])), p=0.95, method = "historical"))
}
mean(var)  # -0.0467
cvar <- c()
for (i in seq(1,p)){
  cvar <- c(cvar,CVaR(as.vector(unlist(r[i])), p=0.95, method = "historical"))
}
mean (cvar)  # -0.0642
# Interprétation : On peut affirmer avec 95% de confiance que les rendements sont supérieurs à -0.0467
# Dans le cas (de probabilité 5%) où le rendement est inférieur à -0.0467, on s'attend à un rendement de -0.0642
# Ces valeurs sont relativement élevées pour des rendements hebdomadaires
# On envisagera dans la suite comment intégrer la minimisation de la VaR dans le programme d'optimisation de Markowitz



#######  Backtest du portefeuille de tangence  #######
backtest <- function(weights,periode){
  # Vt = valeur du portefeuille à la date t = X.Pt avec X=nombre d'actions détenues de chaque titre du portefeuille
  if (periode==1){
    X <- as.vector(unlist(c(100*weights/df[1,retenus])))
    V <- apply(df[,retenus], 1, FUN = function(x){X %*% x})
    return(V <- V*100/V[1])
  }
  else{
    X <- as.vector(unlist(c(100*weights/df_test[1,retenus])))
    V <-apply(df_test[,retenus], 1, FUN = function(x){X %*% x})
    return(V <- V*100/V[1])
  }
}

backtest_all <- function(weights,periode){
  # Vt = valeur du portefeuille à la date t = X.Pt avec X=nombre d'actions détenues de chaque titre du portefeuille
  X <- as.vector(unlist(c(100*weights/df[1,-1])))
  if (periode==1){
    return(apply(df[,-1], 1, FUN = function(x){X %*% x}))
  }
  else{
    V <-apply(df_test[,-1], 1, FUN = function(x){X %*% x})
    return(V <- V*100/V[1])
  }
}
# On calcule le portefeuille tangent sur la période 2012-2015
er <- colMeans(r)
tp <- tangency.portfolio(er, var(r),rf, shorts=T)
tp_noshort <- tangency.portfolio(er, var(r),rf, shorts=F)
# On le teste sur la même période utilisée pour le calculer:
V <- backtest(tp$weights,1)
V_noshort <- backtest(tp_noshort$weights,1)
d <- data.frame(Date = as.Date(df[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_cartesian(ylim=c(50,350))+
  ggtitle("Valeur du portefeuille tangent sur la période d'entrainement")+
  theme(legend.title = element_blank())
# Rendement espéré :
exp(as.double( log(V[n+1]/V[1])/n*52))-1
exp(as.double( log(V_noshort[n+1]/V[1])/n*52))-1
#vol, cvar:
sd(log(V[2:(n+1)]/V[1:n]))
CVaR(log(V[2:(n+1)]/V[1:n]), p=0.95, method = "historical")
sd(log(V_noshort[2:(n+1)]/V_noshort[1:n]))
CVaR(log(V_noshort[2:(n+1)]/V_noshort[1:n]), p=0.95, method = "historical")



# On le teste sur une nouvelle période hors-échantillon :
V <- backtest(tp$weights,2)
V_noshort <- backtest(tp_noshort$weights,2)
d <- data.frame(Date = as.Date(df_test[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
d$variable <- factor(d$variable, levels = c("Sans.vad","Avec.vad"))
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  coord_cartesian(ylim=c(50,150))+
  ggtitle("Valeur du portefeuille tangent sur la période de test")+
  theme(legend.title = element_blank())
# Rendement réel sur la période de test :
exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
exp(as.double( log(V_noshort[nt+1]/V_noshort[1])/nt*52))-1
#vol, cvar:
sd(log(V[2:(nt+1)]/V[1:nt])) 
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")
sd(log(V_noshort[2:(nt+1)]/V_noshort[1:nt])) 
CVaR(log(V_noshort[2:(nt+1)]/V_noshort[1:nt]), p=0.95, method = "historical")



# Comparaison avec un portefeuille equi-réparti:
V <- backtest(1,1)
plot(V, type="l", main = "Valeur du portefeuille équi-réparti sur la période 2012-2015")
exp(as.double( log(V[n+1]/V[1])/n*52))-1
sd(log(V[2:(n+1)]/V[1:n]))
CVaR(log(V[2:(n+1)]/V[1:n]), p=0.95, method = "historical")
# En période de test :
V <- backtest(1/length(retenus),2)
plot(V, type="l", main = "Valeur du portefeuille équi-réparti sur la période 2016-2017")
d <- data.frame(Date = as.Date(df_test[,1]),"Valeur du portefeuille"=V) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.4)+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.4)+
  coord_cartesian(ylim=c(80,150))+
  ggtitle("Evolution de la valeur du portefeuille équi-réparti sur la période de test")+
  theme(legend.position = "none")
# Rendement :
exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
sd(log(V[2:(nt+1)]/V[1:nt]))
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")
# --> Sur la période de test, le portefeuille equi-réparti a fourni un rendement et une volatilité meilleurs que celui de Markowitz


####  EWMA  ####

# Calcul de la matrice de variance avec volatilité EWMA
ewma_var <- function(r,lambda=0.985){
  n <- nrow(r)
  Sigma <- cov(r)
  mu <- colMeans(r)
  rdev <- sweep(r, 2, mu, "-") # déviations de la moyenne
  for (i in 1:length(r[,1])) {
    ri <- as.double(rdev[i, ])
    r2 <- outer(ri, ri)
    Sigma <- (1 - lambda)/(1 - lambda^n)*r2 + lambda*Sigma
  }
  return(Sigma)
}

ewma_er <- function(r,lambda=0.985){
  mu <- colMeans(r)
  for (i in 1:length(r[,1])) {
    ri <- as.double(r[i, ])
    mu <- (1 - lambda)/(1 - lambda^n) * ri + lambda * mu
  }
  return(mu)
}

plot_ewma <- function(r,lambda=0.985){
  y <- c(mean(r[,1]))
  z <- c(mean(r[,1]))
  mu <- mean(r[,1])
  for (i in 1:length(r[,1])) {
    ri <- as.double(r[i, 1])
    mu <- (1 - lambda) * ri + lambda * mu
    y <- c(y, mu)
    z <- c(z, mean(r[1:i,1]))
  }
  plot(z, type="l", lwd=1.5, ylim=c(-0.01,0.02), main="Estimation du rendement par EWMA", xlab="Date", ylab="Rendement estimé")
  lines(y,type="l", col=colors()[92], lwd=2)
  lambda <- 0.9
  y <- c(mean(r[,1]))
  mu <- mean(r[,1])
  for (i in 1:length(r[,1])) {
    ri <- as.double(r[i, 1])
    mu <- (1 - lambda) * ri + lambda * mu
    y <- c(y, mu)
  }
  lines(y,type="l",col=colors()[144], lwd=2)
  lambda <- 0.995
  y <- c(mean(r[,1]))
  mu <- mean(r[,1])
  for (i in 1:length(r[,1])) {
    ri <- as.double(r[i, 1])
    mu <- (1 - lambda) * ri + lambda * mu
    y <- c(y, mu)
  }
  lines(y,type="l", col=colors()[507], lwd=2)
  legend("topright",
         legend = c("EWMA(0.9)", "EWMA(0.985)","EWMA(0.995)","Moyenne empirique"),
         col = c(colors()[144], colors()[92],colors()[507],"black"), lty = 1, lwd=c(2,2,2,1))
}
plot_ewma(r)

var_ewma <- ewma_var(r)
er_ewma <- ewma_er(r)

# Stabilité des portefeuilles efficients
er1 <- ewma_er(r)
er2 <- ewma_er(r2)
matrice_variance1 <- ewma_var(r)
matrice_variance2 <- ewma_var(r2)
tp1 <- tangency.portfolio(er1, matrice_variance1,rf)$weights
tp2 <- tangency.portfolio(er2, matrice_variance2,rf)$weights
x <- seq(1,15)
d <- data.frame(tp1, tp2, x)
d <- pivot_longer(d, cols=c('tp1', 'tp2'), names_to='Periode', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Periode)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Pondérations du portefeuille tangent en fonction de la période de calcul") +
  scale_fill_discrete(name = "Periode", labels = c("4 ans", "2 ans"))
# Certains poids restent stables mais la plupart sont modifiés



### Backtest EWMA ###
# Calcul du portefeuille tangent
tp <- tangency.portfolio(er_ewma, var_ewma,rf, shorts = T)
tp_noshort <- tangency.portfolio(er_ewma, var_ewma,rf, shorts = FALSE)
V <- backtest(tp$weights, 1)
V_noshort <- backtest(tp_noshort$weights, 1)
#rendement:
exp(as.double( log(V[n+1]/V[1])/n*52))-1
exp(as.double( log(V_noshort[n+1]/V[1])/n*52))-1
#vol, cvar:
sd(log(V_noshort[2:(n+1)]/V_noshort[1:n])) 
CVaR(log(V_noshort[2:(n+1)]/V_noshort[1:n]), p=0.95, method = "historical")

d <- data.frame(Date = as.Date(df[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_cartesian(ylim=c(50,450))+
  ggtitle("Valeur du portefeuille tangent sur la periode d'entrainement")
# En periode de test :
V <- backtest(tp$weights, 2)
V_noshort <- backtest(tp_noshort$weights, 2)
exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
exp(as.double( log(V_noshort[nt+1]/V_noshort[1])/nt*52))-1
#vol, cvar:
sd(log(V_noshort[2:(nt+1)]/V_noshort[1:nt])) 
CVaR(log(V_noshort[2:(nt+1)]/V_noshort[1:nt]), p=0.95, method = "historical")

d <- data.frame(Date = as.Date(df_test[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
d$variable <- factor(d$variable, levels = c("Sans.vad","Avec.vad"))
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  coord_cartesian(ylim=c(50,150))+
  ggtitle("Valeur du portefeuille tangent sur la période de test")

d <- data.frame(Date = as.Date(df_test[,1]),"Valeur du portefeuille"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.4)+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.4)+
  coord_cartesian(ylim=c(80,150))+
  ggtitle("Evolution de la valeur du portefeuille tangent EWMA sur la période de test")+
  theme(legend.position = "none")

####  GARCH  ####
# Calcul de la matrice de variance avec GARCH
garch_var <- function(r){
  # GARCH(1,1) normal univarié pour chaque actif 
  garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "norm") 
  
  # On utilise GARCH(1,1) pour les corrélations conditionnelles
  dcc.garch11.spec = dccspec(uspec = multispec( replicate(dim(r)[2], garch11.spec) ), dccOrder = c(1,1), distribution = "mvnorm") 
  dcc.fit = dccfit(dcc.garch11.spec, data = data.frame(r))
  
  # Prédiction de la volatilité conditionnelle et des corrélations
  dcc.fcst = dccforecast(dcc.fit, n.ahead=3) 
  
  # Prédiction de la matrice de covariance GARCH
  dcc.fcst@mforecast$H  
  m=dcc.fcst@mforecast$H
  cov.GARCH=m[[1]][,,1]
  rownames(cov.GARCH) <- colnames(r)
  colnames(cov.GARCH) <- colnames(r)
  return(cov.GARCH)
}

var_garch=garch_var(r)

er <- colMeans(r)
# Calcul du portefeuille tangent avec volatilité GARCH
tp <- tangency.portfolio(er, var_garch,rf, shorts = T)
tp_noshort <- tangency.portfolio(er, var_garch,rf, shorts = FALSE)
V <- backtest(tp$weights, 1)
V_noshort <- backtest(tp_noshort$weights, 1)
#rendement:
exp(as.double( log(V[n+1]/V[1])/n*52))-1
exp(as.double( log(V_noshort[n+1]/V_noshort[1])/n*52))-1
#vol, cvar:
sd(log(V_noshort[2:(n+1)]/V_noshort[1:n])) 
CVaR(log(V_noshort[2:(n+1)]/V_noshort[1:n]), p=0.95, method = "historical")
d <- data.frame(Date = as.Date(df[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_cartesian(ylim=c(50,400))+
  ggtitle("Valeur du portefeuille tangent sur la période d'entrainement")
# En période de test :
V <- backtest(tp$weights, 2)
V_noshort <- backtest(tp_noshort$weights, 2)
exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
exp(as.double( log(V_noshort[nt+1]/V_noshort[1])/nt*52))-1
#vol, cvar:
sd(log(V_noshort[2:(nt+1)]/V_noshort[1:nt])) 
CVaR(log(V_noshort[2:(nt+1)]/V_noshort[1:nt]), p=0.95, method = "historical")
d <- data.frame(Date = as.Date(df_test[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
d$variable <- factor(d$variable, levels = c("Sans.vad","Avec.vad"))
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  coord_cartesian(ylim=c(50,150))+
  ggtitle("Valeur du portefeuille tangent sur la période de test")

d <- data.frame(Date = as.Date(df_test[,1]),"Valeur du portefeuille"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.4)+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.4)+
  coord_cartesian(ylim=c(80,150))+
  ggtitle("Evolution de la valeur du portefeuille tangent GARCH sur la période de test")+
  theme(legend.position = "none")
# Rendement espéré :
tp$er*52
# Rendement réel :
exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
# Volatilité espérée :
tp$sd
# Volatilité réelle :
sd(log(V[2:(nt+1)]/V[1:nt]))


# Sensibilité GARCH
er1 <- colMeans(r)
er2 <- colMeans(r2)
matrice_variance1 <- garch_var(r)
matrice_variance2 <- garch_var(r2)
tp1 <- tangency.portfolio(er1, matrice_variance1,rf, shorts=FALSE)$weights
tp2 <- tangency.portfolio(er2, matrice_variance2,rf, shorts=FALSE)$weights
x <- seq(1,15)
d <- data.frame(tp1, tp2, x)
d <- pivot_longer(d, cols=c('tp1', 'tp2'), names_to='Periode', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Periode)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Pondérations du portefeuille tangent en fonction de la période de calcul") +
  scale_fill_discrete(name = "Periode", labels = c("4 ans", "2 ans"))


### Stabilité des poids du portefeuille tangent :
rolling_weights <-  function(x, voltype) {
  ch_dates <- x[,1]
  x <- x[,-1]
  class(x) <- "numeric"
  if(voltype=="GARCH"){
    ef = tangency.portfolio(colMeans(x), garch_var(x),rf,shorts = F)
    return(ef$weights)
  }
  
  else{ 
    if (voltype=="EWMA") {
      ef = tangency.portfolio(ewma_er(x,0.985), ewma_var(x,0.985),rf,shorts = F)
      return(ef$weights)
    }
    
    else {
      if (voltype=="CVAR") {
        returns <- timeSeries(x,charvec=ch_dates)
        ef <- tangencyPortfolio(data =returns, spec = cvar, constraints = "LongOnly")
      return(c(getWeights(ef)))
    }
      
      else { 
        if (voltype=="Contraction") {
          ef = tangency.portfolio(colMeans(x),(var(x)+diag(15))/2,rf,shorts = F)
          return(ef$weights)
        }
        
        else { 
      ef = tangency.portfolio(colMeans(x), var(x),rf,shorts = F)
      return(ef$weights)
    }
  }
}}
}
 plot_rolling_weights <- function(voltype){
  weights <- rollapply(cbind(dates,rbind(r,rt)), width=207, by.column=FALSE, align="right", FUN=rolling_weights, voltype)
  x <- as.vector(df_test[,1])
  d <- data.frame(x,weights)
  d <- gather(d, "Action",'Poids',-x)
  ggplot(d, aes(x=as.Date(x), y=Poids, fill=Action)) +
    geom_bar(stat='identity')+ scale_x_date() +
    labs(x = "Date",title = paste("Evolution des poids avec",voltype))
 }
dates <- c(paste(df[1:207,1]),paste(df_test[1:103,1]))
plot_rolling_weights("simple")
plot_rolling_weights("EWMA")  # Moins stable que Markowitz standard
# plot_rolling_weights("GARCH") # (prend qq min pour s'exécuter) Légèrement plus stable que le modèle Markowitz standard
plot_rolling_weights("CVAR")
plot_rolling_weights("Contraction")


# Détermination de la valeur de lambda
x <- c()
var_sd_2 <- var(rt)
for (i in seq(0.95,0.995, 0.001)){
  var_ewma_1 <- ewma_var(r,i)
  x <- c(x,Matrix::norm((var_ewma_1-var_sd_2), type="f"))
}
plot(seq(0.95,0.995, 0.001), x, main="Distance entre la matrice EWMA et la matrice de variance de test", xlab="lambda", ylab="distance")


er2 <- colMeans(rt)
x <- c()
for (i in seq(0.95,0.995, 0.001)){
  x <- c(x,norm((ewma_er(r, i)-er2), type="2"))
}
plot(seq(0.95,0.995, 0.001), x, main="Distance entre le rendement EWMA et le rendement moyen de test", xlab="lambda", ylab="distance")
# lambda = 0.985 semble être la valeur optimale


### Comparaison des matrices de covariance
var_sd_1 <- var(r)
var_sd_2 <- var(rt)

var_ewma_1 <- ewma_var(r)
var_ewma_2 <- ewma_var(rt)

var_garch_1 <- garch_var(r)
var_garch_2 <- garch_var(rt)

levelplot(var_sd_2, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0005,0.0015,length.out = 98),Inf))
Matrix::norm(var_sd_2,type="f")   # 0.00647
levelplot(var_sd_1, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0005,0.0015,length.out = 98),Inf))
Matrix::norm(var_sd_1,type="f")   # 0.00862
levelplot(var_ewma_1, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0005,0.0015,length.out = 98),Inf))
Matrix::norm(var_ewma_1,type="f")   # 0.00931
levelplot(var_garch_1, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0005,0.0015,length.out = 98),Inf))
Matrix::norm(var_garch_1,type="f")   # 0.00789


# Variations par rapport à la matrice de covariance standard de la période 2
# Variations relatives
levelplot((var_sd_1-var_sd_2)/var_sd_2, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-1.5,1.5,length.out = 98),Inf))
Matrix::norm((var_sd_1-var_sd_2)/var_sd_2,type="f")   #  9.12
levelplot((var_ewma_1-var_sd_2)/var_sd_2, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-1.5,1.5,length.out = 98),Inf))
Matrix::norm((var_ewma_1-var_sd_2)/var_sd_2, type="f")   # 12.77037
levelplot((var_garch_1-var_sd_2)/var_sd_2, col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-1.5,1.5,length.out = 98),Inf))
Matrix::norm((var_garch_1-var_sd_2)/var_sd_2, type="f")  # 8.797167

# Variations absolues
levelplot((var_sd_1-var_sd_2), col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0008,0.0008,length.out = 98),Inf))
Matrix::norm((var_sd_1-var_sd_2),type="f")/Matrix::norm(var_sd_2,type="f")   # 0.51
levelplot((var_ewma_1-var_sd_2), col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0008,0.0008,length.out = 98),Inf))
Matrix::norm((var_ewma_1-var_sd_2), type="f")/Matrix::norm(var_sd_2,type="f")   # 0.61
levelplot((var_garch_1-var_sd_2), col.regions=colorRampPalette(c("navyblue", "seagreen1", "firebrick1"))(100), at=c(-Inf,seq(-0.0008,0.0008,length.out = 98),Inf))
Matrix::norm((var_garch_1-var_sd_2), type="f")/Matrix::norm(var_sd_2,type="f")   # 0.44



# Comparaison des rendements espérés et des rendements réels
x <- seq(1,15)
d <- data.frame(colMeans(r), colMeans(rt), x)
colnames(d) <- c("er","ert","x") 
d <- pivot_longer(d, cols=c('er', 'ert'), names_to='Periode', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Periode)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Rendements moyens des actifs à chaque période") +
  scale_fill_discrete(name = "Periode", labels = c("1", "2"))

x <- seq(1,15)
d <- data.frame(er_ewma, colMeans(rt), x)
colnames(d) <- c("er","ert","x") 
d <- pivot_longer(d, cols=c('er', 'ert'), names_to='Periode', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Periode)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Rendements moyens des actifs à chaque période") +
  scale_fill_discrete(name = "Periode", labels = c("1", "2"))

# Différence avec la période de test
x <- seq(1,15)
d <- data.frame(er_ewma-colMeans(rt), colMeans(r)-colMeans(rt), x)
colnames(d) <- c("er","ert","x") 
d <- pivot_longer(d, cols=c('er', 'ert'), names_to='Méthode', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Méthode)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Différence avec la période de test") +
  scale_fill_discrete(name = "Méthode", labels = c("EWMA", "Moyenne"))

# Test Kolmogorov-Smirnov d'adéquation à la même loi
for (i in 1:15){
  suppressWarnings(print(ks.test(r[,i],rt[,i])$p.value))
}
# On ne rejette pas l'hypothèse que les rendements en période de test suivent la même loi que ceux de la première période


# Traçage des 3 frontières efficientes
var_ewma <- ewma_var(r)
er_ewma <- ewma_er(r)
er <- colMeans(r)
var_garch <- garch_var(r)
plot_EFs <- function()
{
  object1 <- efficient.frontier(er, var(r), alpha.min=-2, alpha.max=1.5, nport=20)
  object2 <- efficient.frontier(er_ewma, var_ewma, alpha.min=-2, alpha.max=1.5, nport=20)
  object3 <- efficient.frontier(er, var_garch, alpha.min=-2, alpha.max=1.5, nport=20)
  tp1 <- tangency.portfolio(er, var(r),rf)
  tp2 <- tangency.portfolio(er_ewma, var_ewma,rf)
  tp3 <- tangency.portfolio(er, var_garch,rf)
  y.lim=c(0,max(object1$er, object2$er, object3$er))
  x.lim=c(0,max(object1$sd, object2$sd, object3$sd))
  plot(object1$sd, object1$er, col = "green", type="b",xlim=x.lim, ylim=y.lim,
       xlab="Portfolio SD", ylab="Portfolio ER",
       main="Efficient Frontier")
  lines(object2$sd,object2$er, col = "dodgerblue2", type="b", pch=4, xlim=x.lim, ylim=y.lim,
        xlab="Portfolio SD", ylab="Portfolio ER")
  lines(object3$sd,object3$er, col = "coral2", type="b",xlim=x.lim, ylim=y.lim,
        xlab="Portfolio SD", ylab="Portfolio ER", pch=0)
  abline(a=rf, b=((tp1$er-rf)/tp1$sd), lwd=0.5, col="green")
  abline(a=rf, b=((tp2$er-rf)/tp2$sd), lwd=0.5, col="dodgerblue2")
  abline(a=rf, b=((tp3$er-rf)/tp3$sd), lwd=0.5, col="coral2")
  abline(v=0, col="lightgrey")
  abline(h=0, col="lightgrey")
  text(0.002, rf, labels="A.S.R.")
  legend("topleft",
         legend = c("EWMA", "GARCH","Markowitz"),
         col = c("dodgerblue2", "coral2","green"), lty = 1)
}
plot_EFs()


#### Optimisation Moyenne-CVaR ####


#donnees <- df[,retenus]
#rownames(donnees) <- as.vector(df[,1])
#donnees <- rbind(df[,retenus], df_test[,retenus])
#rownames(donnees) <- c(as.vector(df[,1]),as.vector(df_test[,1]))
returns <- timeSeries(r,charvec=paste(df[1:207,1]))

cvar <- portfolioSpec()
setType(cvar) <- "CVaR"
setAlpha(cvar) <- 0.25
suppressWarnings(setSolver(cvar) <- "solveRglpk.CVAR")
setRiskFreeRate(cvar) <- rf
cvarpf <- tangencyPortfolio(data =returns, spec = cvar, constraints = "LongOnly")

standard <- portfolioSpec()
setSolver(standard) = "solveRquadprog"
setRiskFreeRate(standard) <- rf
standardpf <- tangencyPortfolio(data = returns, spec = standard, constraints = "LongOnly")

wCVAR <- c(getWeights(cvarpf))

V_noshort <- backtest(wCVAR, 1)
#rendement:
exp(as.double( log(V_noshort[n+1]/V_noshort[1])/n*52))-1
#vol, cvar:
sd(log(V_noshort[2:(n+1)]/V_noshort[1:n])) 
CVaR(log(V_noshort[2:(n+1)]/V_noshort[1:n]), p=0.95, method = "historical")

d <- data.frame(Date = as.Date(df[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_cartesian(ylim=c(50,450))+
  ggtitle("Valeur du portefeuille tangent sur la période d'entrainement")
# En période de test :
V_noshort <- backtest(wCVAR, 2)
exp(as.double( log(V_noshort[nt+1]/V_noshort[1])/nt*52))-1
#vol, cvar:
sd(log(V_noshort[2:(nt+1)]/V_noshort[1:nt])) 
CVaR(log(V_noshort[2:(nt+1)]/V_noshort[1:nt]), p=0.95, method = "historical")

d <- data.frame(Date = as.Date(df_test[,1]),"Valeur du portefeuille"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.4)+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.4)+
  coord_cartesian(ylim=c(80,150))+
  ggtitle("Evolution de la valeur du portefeuille tangent CVaR sur la période de test")+
  theme(legend.position = "none")

wstandard <- c(getWeights(standardpf))
V_standard <- backtest(wstandard, 2)
V_standard <- V_standard*100/V_standard[1]
V <- backtest(wCVAR, 2)
plot(V, type="l", main = "Valeur du portefeuille tangent sur la période de test")
lines(V_standard, col="red")
legend("topleft",
       legend = c("Moyenne-Variance", "Moyenne-CVaR"),
       col = c("red", "black"), lty = 1)

# Rendement du portefeuille CVAR sur la période de test :
exp(log(V[nt+1]/V[1])/nt*52)-1
# Volatilité du portefeuille CVAR sur la période de test :
sd(log(V[2:(nt+1)]/V[1:nt]))

tp_noshort <- tangency.portfolio(er_ewma, var_ewma,rf, shorts = FALSE)
V_noshort <- backtest(tp_noshort$weights, 1)
d <- data.frame(Date = as.Date(df[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_cartesian(ylim=c(50,450))+
  ggtitle("Valeur du portefeuille tangent sur la période d'entrainement")
# En période de test :
V <- backtest(tp$weights, 2)
V_noshort <- backtest(tp_noshort$weights, 2)
d <- data.frame(Date = as.Date(df_test[,1]),"Avec vad"=V, "Sans vad"=V_noshort) %>%
  gather(key = "variable", value = "value", -Date)
d$variable <- factor(d$variable, levels = c("Sans.vad","Avec.vad"))
ggplot(d, mapping = aes(x = Date, y = value, fill = variable)) + 
  geom_area(alpha = 0.35, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_line(colour = 'black', size = 1, alpha = 0.35)+
  coord_cartesian(ylim=c(50,150))+
  ggtitle("Valeur du portefeuille tangent sur la période de test")

# Performance relative
perf_rel <- (V - V_standard) / V_standard * 100
plot(perf_rel, type = "h", xlab = "",
     ylab = "Pourcent",
     main = "Performance relative Moyenne-CVaR vs. Moyenne-Variance")
abline(h = 0, col = "grey")


CVaR(log(V_standard[2:(nt+1)]/V_standard[1:nt]), p=0.95, method = "historical")
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")


# Détermination du paramètre alpha optimal
res <- c()
res_cvar <- c()
for (alpha in seq(0.05,0.8,0.05)){
  setAlpha(cvar) <- alpha
  cvarpf <- tangencyPortfolio(data =returns, spec = cvar, constraints = "LongOnly")
  wCVAR <- c(getWeights(cvarpf))
  V_CVAR <- backtest(wCVAR, 1)
  res <- c(res, exp(log(V_CVAR[n+1]/V_CVAR[1])/n*52)-1)
  res_cvar <- c(res_cvar, CVaR(log(V_CVAR[2:(nt+1)]/V_CVAR[1:nt]), p=0.95, method = "historical") )
}
plot(seq(0.05,0.8,0.05), res,type="l", xlab="alpha", ylab="Rendement", main ="Rendement annualisé en fonction de alpha")
plot(seq(0.05,0.8,0.05), res_cvar,type="l", xlab="alpha", ylab="CVaR à 5%", main ="CVaR du rendement en fonction de alpha")


# Frontière efficiente
frontiere <- portfolioFrontier(returns, spec = cvar, constraints = "LongOnly",
                  include.mvl = TRUE)
frontierPlot(frontiere)
frontierPlot(frontiere, add=T, type="l")
tp <- tangencyPoints(frontiere, add=TRUE, pch=18, col="blue")
abline(a=rf, b=((tp[2]-rf)/tp[1]), lwd=0.5, col="dodgerblue2")

plot(1,1,main="Frontière efficiente Moyenne-CVaR")



### Méthode shrinkage ###

# 1. Shrinkage des espérances de rendement
# L'idée est de remplacer les estimations des rendements par (estimations+moyenne(esimations))/2
# pour atténuer les valeurs extremes et les rapprocher de la moyenne

# On calcule le portefeuille tangent sur la période 2012-2015
tp <- tangency.portfolio(colMeans(r), 1/2*var(r)+1/2*diag(15),rf, shorts=FALSE)
V <- backtest(tp$weights,1)
plot(V, type="l", main = "Valeur du portefeuille tangent sur la période 2016-2017")
exp(as.double( log(V[n+1]/V[1])/n*52))-1
#vol, cvar:
sd(log(V[2:(n+1)]/V[1:n])) 
CVaR(log(V[2:(n+1)]/V[1:n]), p=0.95, method = "historical")

# On le teste sur une nouvelle période :
V <- backtest(tp$weights,2)
plot(V, type="l", main = "Valeur du portefeuille tangent sur la période 2016-2017")
exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
#vol, cvar:
sd(log(V[2:(nt+1)]/V[1:nt])) 
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")

ef <- efficient.frontier(colMeans(r), 1/2*var(r)+1/2*diag(15), alpha.min=-2, alpha.max=1.5, nport=20)
plot(ef, plot.assets=T)

# Sensibilité Shrinkage
er1 <- colMeans(r)
er2 <- colMeans(r2)
matrice_variance1 <- 1/2*var(r)+1/2*diag(15)
matrice_variance2 <- 1/2*var(r2)+1/2*diag(15)
tp1 <- tangency.portfolio(er1, matrice_variance1,rf, shorts=FALSE)$weights
tp2 <- tangency.portfolio(er2, matrice_variance2,rf, shorts=FALSE)$weights
x <- seq(1,15)
d <- data.frame(tp1, tp2, x)
d <- pivot_longer(d, cols=c('tp1', 'tp2'), names_to='Période', 
                  values_to="Poids")
ggplot(d, aes(x=x, y=Poids, fill=Période)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Actions",title = "Pondérations du portefeuille tangent en fonction de la période de calcul") +
  scale_fill_discrete(name = "Période", labels = c("4 ans", "2 ans"))



######  Données du S&P500  ######
# On ajoute les actions du S&P500
# les données de cette partie sont suffixées par 'a' (comme données Augmentées)
dfa <- read.csv("cac_sp500.csv", header = TRUE)
dfa_test <- dfa[209:312,]    # Données de test: de 2016 à 2017
dfa <- dfa[1:208,]           # Données de construction du portefeuille: de 2012 à 2015
na <- nrow(dfa)-1
pa <- ncol(dfa)
nta <- nrow(dfa_test)-1

# Calcul des rendements logarithmiques
ra <- log(dfa[2:(na+1),2:pa]/dfa[1:na,2:pa])
rownames(ra) <- NULL
rta <- log(dfa_test[2:(nta+1),2:pa]/dfa_test[1:nta,2:pa])
rownames(rta) <- NULL
# Division en 105 portefeuilles en regroupant les actifs les moins corrélés
moins_correles <- rep(list(),210)
for (i in 1:210){
  correlations <- as.table(cor(ra))
  moins_correles[[i]] <- as.data.frame(correlations)[which.min(correlations),1:2]
  ra <- ra[,-c(as.numeric(paste(moins_correles[[i]])))]
}
ra1 <- log(dfa[2:(na+1),2:pa]/dfa[1:na,2:pa])
rta1 <- rta[,colnames(ra)]
ra <- lapply(seq(1,105),function(x){cbind(subset(ra,select=106-x),ra1[,as.vector(unlist(moins_correles[[x]]))],ra1[,as.vector(unlist(moins_correles[[x+105]]))])})
rta <- lapply(seq(1,105),function(x){cbind(subset(rta1,select=106-x),rta[,as.vector(unlist(moins_correles[[x]]))],rta[,as.vector(unlist(moins_correles[[x+105]]))])})

backtest_a <- function(weights,stocks){
  # Vt = valeur du portefeuille à la date t = X.Pt avec X=nombre d'actions détenues de chaque titre du portefeuille
  X <- as.vector(unlist(c(100*weights/dfa[1,stocks])))
  return(apply(dfa_test[,stocks], 1, FUN = function(x){X %*% x}))
}

optim_portefeuille <- function(er,ev){
  tp = tryCatch({
    tangency.portfolio(er, ev,0, shorts=FALSE)
  }, error = function(e) {
    return(efficient.portfolio(er, ev,target.return=0, shorts=FALSE))
  })
  # On le teste sur une nouvelle période :
  V <- backtest_a(tp$weights,names(er))
  # Rendement sur la période de test :
  rendement <- exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
  # Volatilité sur la période de test :
  volatilite <- sd(log(V[2:(nt+1)]/V[1:nt]))
  # CVaR :
  cvar_95 <- CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical") 
  return(c(rendement, volatilite, cvar_95))
}

rendements <- matrix(nrow=105, ncol=6)
colnames(rendements) <- c("Markowitz", "EWMA","GARCH","Shrinkage","CVaR","Equi-réparti")
volatilites <- matrix(nrow=105, ncol=6)
colnames(volatilites) <- colnames(rendements)
cvars <- matrix(nrow=105, ncol=6)
colnames(cvars) <- colnames(rendements)

cvar <- portfolioSpec()
setType(cvar) <- "CVaR"
setAlpha(cvar) <- 0.2
setSolver(cvar) <- "solveRglpk.CVAR"
setRiskFreeRate(cvar) <- 0

for (i in 1:105){
  # Markowitz standard
  er <- colMeans(ra[[i]])
  ev <- var(ra[[i]])
  res <- optim_portefeuille(er, ev)
  rendements[i,1] <- res[1]
  volatilites[i,1] <- res[2]
  cvars[i,1] <- res[3]
  # EWMA
  er <- ewma_er(ra[[i]])
  ev <- ewma_var(ra[[i]])
  res <- optim_portefeuille(er, ev)
  rendements[i,2] <- res[1]
  volatilites[i,2] <- res[2]
  cvars[i,2] <- res[3]
  # GARCH
  er <- colMeans(ra[[i]])
  ev <- garch_var(ra[[i]])
  res <- optim_portefeuille(er, ev)
  rendements[i,3] <- res[1]
  volatilites[i,3] <- res[2]
  cvars[i,3] <- res[3]
  # Shrinkage
  er <- colMeans(ra[[i]])
  ev <- (var(ra[[i]])+diag(5))/2
  res <- optim_portefeuille(er, ev)
  rendements[i,4] <- res[1]
  volatilites[i,4] <- res[2]
  cvars[i,4] <- res[3]
  # CVaR
  returns <- timeSeries(ra[[i]],charvec=paste(dfa[1:207,1]))
  cvarpf <- tangencyPortfolio(data =returns, spec = cvar, constraints = "LongOnly")
  wCVAR <- c(getWeights(cvarpf))
  V <- backtest_a(wCVAR,names(wCVAR))
  rendements[i,5] <- exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
  volatilites[i,5] <- sd(log(V[2:(nt+1)]/V[1:nt]))
  cvars[i,5] <- CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")
  # Equi-réparti
  V <- backtest_a(1/5,names(er))
  rendements[i,6] <- exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
  volatilites[i,6] <- sd(log(V[2:(nt+1)]/V[1:nt]))
  cvars[i,6] <- CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")
  print(i)
}

colMeans(rendements)
colMeans(volatilites)
colMeans(cvars)

lillie.test(rendements[,1])
t.test(rendements[,6], rendements[,4], alternative="greater")

# On inverse les espérances de rendement

for (i in 1:105){
  # Markowitz standard
  colmeans <- colMeans(ra[[i]])
  inverted_means <- 2*mean(colmeans)-colmeans
  er <- inverted_means
  ev <- var(ra[[i]])
  res <- optim_portefeuille(er, ev)
  rendements[i,1] <- res[1]
  volatilites[i,1] <- res[2]
  cvars[i,1] <- res[3]
  # EWMA
  er <- 2*mean(ewma_er(ra[[i]])) - ewma_er(ra[[i]])
  ev <- ewma_var(ra[[i]])
  res <- optim_portefeuille(er, ev)
  rendements[i,2] <- res[1]
  volatilites[i,2] <- res[2]
  cvars[i,2] <- res[3]
  # GARCH
  er <- inverted_means
  ev <- garch_var(ra[[i]])
  res <- optim_portefeuille(er, ev)
  rendements[i,3] <- res[1]
  volatilites[i,3] <- res[2]
  cvars[i,3] <- res[3]
  # Shrinkage
  er <- inverted_means
  ev <- (var(ra[[i]])+diag(5))/2
  res <- optim_portefeuille(er, ev)
  rendements[i,4] <- res[1]
  volatilites[i,4] <- res[2]
  cvars[i,4] <- res[3]
  # CVaR
  returns <- timeSeries(2*colMeans(ra[[i]])-ra[[i]],charvec=paste(dfa[1:207,1]))
  cvarpf <- tangencyPortfolio(data =returns, spec = cvar, constraints = "LongOnly")
  wCVAR <- c(getWeights(cvarpf))
  V <- backtest_a(wCVAR,names(wCVAR))
  rendements[i,5] <- exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
  volatilites[i,5] <- sd(log(V[2:(nt+1)]/V[1:nt]))
  cvars[i,5] <- CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")
  # Equi-réparti
  V <- backtest_a(1/5,names(er))
  rendements[i,6] <- exp(as.double( log(V[nt+1]/V[1])/nt*52))-1
  volatilites[i,6] <- sd(log(V[2:(nt+1)]/V[1:nt]))
  cvars[i,6] <- CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical")
  print(i)
}

colMeans(rendements)
colMeans(volatilites)
colMeans(cvars)

# Portefeuille des 15 actifs retenus
colmeans <- colMeans(r)
inverted_means <- 2*mean(colmeans)-colmeans
er <- inverted_means
er <- colMeans(r)
ev <- var(r)
tp <- tangency.portfolio(er, ev,rf, shorts=FALSE)
# On le teste sur une nouvelle période :
V <- backtest(tp$weights,2)
# Rendement sur la période de test :
as.double( log(V[nt+1]/V[1])/nt*52 )  #0.136 avec shrinkage, # 0.126 sans shrinkage,  # 0.107 avec vad
# Volatilité sur la période de test :
sd(log(V[2:(nt+1)]/V[1:nt]))  # 0.0184, # 0.0204 sans shrinkage # 0.0228 avec vad
# CVaR :
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical") 

# EWMA
er <- 2*mean(ewma_er(r)) - ewma_er(r)
ev <- ewma_var(r)
tp <- tangency.portfolio(er, ev,rf, shorts=FALSE)
# On le teste sur une nouvelle période :
V <- backtest(tp$weights,2)
# Rendement sur la période de test :
as.double( log(V[nt+1]/V[1])/nt*52 )  #0.136 avec shrinkage, # 0.126 sans shrinkage,  # 0.107 avec vad
# Volatilité sur la période de test :
sd(log(V[2:(nt+1)]/V[1:nt]))  # 0.0184, # 0.0204 sans shrinkage # 0.0228 avec vad
# CVaR :
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical") 


# GARCH
er <- inverted_means
ev <- garch_var(r)
tp <- tangency.portfolio(er, ev,rf, shorts=FALSE)
# On le teste sur une nouvelle période :
V <- backtest(tp$weights,2)
# Rendement sur la période de test :
as.double( log(V[nt+1]/V[1])/nt*52 )  #0.136 avec shrinkage, # 0.126 sans shrinkage,  # 0.107 avec vad
# Volatilité sur la période de test :
sd(log(V[2:(nt+1)]/V[1:nt]))  # 0.0184, # 0.0204 sans shrinkage # 0.0228 avec vad
# CVaR :
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical") 

# Shrinkage
er <- (mean(inverted_means)+inverted_means)/2
ev <- (var(r)+diag(15))/2
tp <- tangency.portfolio(er, ev,rf, shorts=FALSE)
# On le teste sur une nouvelle période :
V <- backtest(tp$weights,2)
# Rendement sur la période de test :
as.double( log(V[nt+1]/V[1])/nt*52 )  #0.136 avec shrinkage, # 0.126 sans shrinkage,  # 0.107 avec vad
# Volatilité sur la période de test :
sd(log(V[2:(nt+1)]/V[1:nt]))  # 0.0184, # 0.0204 sans shrinkage # 0.0228 avec vad
# CVaR :
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical") 

#equi-réparti
# On le teste sur une nouvelle période :
V <- backtest(1/15,2)
# Rendement sur la période de test :
as.double( log(V[nt+1]/V[1])/nt*52 )  #0.136 avec shrinkage, # 0.126 sans shrinkage,  # 0.107 avec vad
# Volatilité sur la période de test :
sd(log(V[2:(nt+1)]/V[1:nt]))  # 0.0184, # 0.0204 sans shrinkage # 0.0228 avec vad
# CVaR :
CVaR(log(V[2:(nt+1)]/V[1:nt]), p=0.95, method = "historical") 


# Meilleur paramètre de shrinkage
rendements <- c()
volatilites <- c()
cvars <- c()
for (lambda in seq(0,1,0.05)){
  rendement <- c()
  volatilite <- c()
  Cvar <- c()
  for (i in 1:105){
    inverted_means <- 2*mean(colMeans(ra[[i]]))-colMeans(ra[[i]])
    er <- lambda* mean(inverted_means) + (1-lambda)*inverted_means
    ev <- (1-1/2)*var(ra[[i]])+1/2*diag(5)
    res <- optim_portefeuille(er, ev)
    rendement <- c(rendement, res[1])
    volatilite<- c(volatilite, res[2])
    Cvar <- c(Cvar,res[3])
  }
  rendements <- c(rendements, mean(rendement))
  volatilites <- c(volatilites, mean(volatilite))
  cvars <- c(cvars,mean(Cvar))
}
plot(rendements)
plot(cvars)



