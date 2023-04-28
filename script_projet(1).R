#### Instalation et chargement des Packages utiles
install.packages("zoo")
install.packages("tseries")
install.packages("fUnitRoots")#tests de racine unitaire plus modulables
install.packages("dplyr")
install.packages("Rmisc")
install.packages("ellipse")
install.packages("ellipsis")
install.packages("car")
install.packages("ggplot2")
install.packages("corpcor")
install.packages("plotrix")
install.packages("geomorph")

require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
require(fUnitRoots)#tests de racine unitaire plus modulables
require(dplyr)
require(ellipse)
require(ellipsis)
require(car)
library(ellipse)
library(Rmisc)
library(forecast)
library(ggplot2)
library(corpcor)
library(plotrix)
library(geomorph)

#### Partie I: Les donnees

# Importation des donnees
path <- "C:/Users/Diego Renaud/OneDrive/Bureau/ENSAE/SerieTemp"
setwd(path)

datafile <- "data1.csv"
data <- read.csv(datafile, sep=";")

# Formatage base de donnees
dates_char <- as.character(data$Periode)
dates_char[1] 
tail(dates_char,1) 
dates <- as.yearmon(seq(from=2000+1/12, to=2023+1/12, by=1/12)) 
xm.source <- zoo(data$Source, order.by=dates)
xm.source <- rev(xm.source)
T <- length(xm.source)
xm <- head(xm.source, T - 4)
dates <- head(dates, T - 4)

# Verifier si la serie temporelle contient des valeurs manquantes (NA)
if (any(is.na(xm))) {
  print("La série temporelle contient des valeurs manquantes")
  
  # Remplacer les valeurs manquantes par la moyenne de la serie temporelle
  xm[is.na(xm)] <- mean(xm, na.rm = TRUE)
} else {
  print("La série temporelle ne contient pas de valeurs manquantes")
}

# Verifier si la série temporelle contient des valeurs non numeriques
if (any(!is.numeric(xm))) {
  print("La série temporelle contient des valeurs non numériques")
  
  # Convertir les valeurs non numériques en valeurs numeriques
  xm <- as.numeric(as.character(xm))
} else {
  print("La série temporelle ne contient pas de valeurs non numériques")
}

# Representation graphique serie source xm 
par(mfrow=c(2,1))  
par(mar=c(2,2,2,2))  # ajuster les marges à gauche, droite, haut et bas
plot(xm)
acf(xm)
# Ca n'a pas l'air statio mais on test quand meme

summary(lm(xm ~ dates))

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
series <- xm; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}
adf <- adfTest_valid(xm,24,adftype="ct")
adf
# On a une p-value de 0.17 donc pas stationnaire 

# Analyse dxm (serie differenciee)
dxm <- diff(xm,1)
summary(lm(dxm ~ dates[-1]))
adf <- adfTest_valid(dxm,24,"nc")
adf
# La p-value est 0.01 donc notre serie differenciee est bien stationnaire



# Partie II: Modeles ARIMA(p,d,q)
dev.off()
par(mfrow=c(1,2))
pacf(dxm,24);acf(dxm,24) #on regarde jusqu'a deux ans de retard

pmax=4;qmax=5

signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

## fonction pour estimer un arima et en verifier l'ajustement et la validite
modelchoice <- function(p,q,data=dxm, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}


armamodels <- armamodelchoice(pmax,qmax) #estime tous les arima (patienter...)

selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec
### On a 2 modeles bien ajustes et valides: ARMA(4,0) et ARMA(4,4)

pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #cree une liste des ordres p et q des modeles candidats
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les elements de la liste
models <- lapply(pqs, function(pq) arima(dxm,c(pq[["p"]],0,pq[["q"]]))) #cree une liste des modeles candidats estimes
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et BIC des modeles candidats
### L'ARMA(4,0) minimise les criteres d'information AIC et BIC



### Question 5: 

# Question 5
# Ajuster le modèle ARIMA(4,1,0)
model <- arima (xm ,c(4,1,0),include.mean=F)
model
# Calculer la somme des carrés des résidus et de la variance totale
ss_res <- sum(model$residuals^2)
ss_tot <- sum((xm-mean(xm))^2)

# Calculer le nombre de paramètres dans le modèle
p <- length(model$coef)

# Calculer le R carré ajusté
adj_r2 <- 1 - (ss_res/(length(xm)-p))/(ss_tot/(length(xm)-1))
adj_r2


# Partie 3: Prévisions
# Question 7
arma40 <- arima(dxm,c(4,0,0),include.mean=F)
arma40
par(mar = c(3, 3, 2, 1))

plot(density(arma40$residuals), xlim=c(-10, 10), main="Densité des résidus", xlab="Valeurs prises", lwd=0.5)

mu <- mean(arma40$residuals)
sigma <- sd(arma40$residuals)
x <- seq(-10, 10)
y <- dnorm(x, mu, sigma)
lines(x, y, lwd=0.5, col="blue")

# Extraction des coeffs du modèle et de la variance des résidus
arma40$coef
phi_1 <- as.numeric(arma40$coef[1])
phi_2 <- as.numeric(arma40$coef[2])
phi_3 <- as.numeric(arma40$coef[3])
phi_4 <- as.numeric(arma40$coef[4])
sigma <- as.numeric(arma40$sigma2)
phi_1
phi_2
phi_3
phi_4
sigma

# Créez le polynôme caractéristique
polynomial_coeffs <- c(1, -phi_1, -phi_2, -phi_3, -phi_4)

# Calculez les racines du polynôme
roots <- polyroot(polynomial_coeffs)

# Affichez les racines
print(roots)

root_modulus <- Mod(roots)

# Vérifiez si toutes les racines sont à l'extérieur du cercle unité
stationary <- all(root_modulus > 1)

# Affichez le résultat
if (stationary) {
  cat("Le processus est stationnaire.\n")
} else {
  cat("Le processus n'est pas stationnaire.\n")
}

#Question 8
pred <- predict(arma40, n.ahead = 2)
XT1 <- pred$pred[1]
XT2 <- pred$pred[2]

fore <- forecast(arma40, h = 5, level = 95)
par(mfrow = c(1, 1))
plot(fore, col = 1, fcol = 2, shaded = TRUE, xlab = "Temps", ylab = "Valeur", main = "Prévision pour la série originale S")

lower_bound <- pred$se * qnorm(0.025) # Bornes inférieures
upper_bound <- pred$se * qnorm(0.975) # Bornes supérieures

# Créer un data.frame pour les prévisions et les intervalles de confiance
df <- data.frame(Time = c("XT+1", "XT+2"),
                 Value = c(XT1, XT2),
                 Lower = c(XT1 + lower_bound[1], XT2 + lower_bound[2]),
                 Upper = c(XT1 + upper_bound[1], XT2 + upper_bound[2]))

# Tracer les prévisions futures (XT+1, XT+2) et les intervalles de confiance à 95%
library(ggplot2)
plot_confidence_intervals <- ggplot(df, aes(x = Time, y = Value)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, size = 1) +
  ggtitle("Prévision pour la série originale S et intervalles de confiance à 95%") +
  xlab("Temps") +
  ylab("Valeur") +
  theme_bw()

# Afficher le graphique
print(plot_confidence_intervals)

if (!requireNamespace("ellipse", quietly = TRUE)) {
  install.packages("ellipse")
}
library(ellipse)

# Obtenir les prévisions pour XT+1 et XT+2
XT1 <- predict(arma40, n.ahead = 2)$pred[1]
XT2 <- predict(arma40, n.ahead = 2)$pred[2]

# Estimer les paramètres du modèle AR(4)
arma_coef <- coef(arma40)
phi_1 <- arma_coef[1]
phi_2 <- arma_coef[2]
phi_3 <- arma_coef[3]
phi_4 <- arma_coef[4]
sigma2 <- arma40$sigma2

# Construire la matrice de covariance pour un modèle AR(4)
Sigma <- matrix(c(sigma2*(1+phi_1^2+phi_2^2+phi_3^2+phi_4^2), sigma2*(phi_1+phi_1*phi_2+phi_1*phi_3+phi_1*phi_4),
                  sigma2*(phi_1+phi_1*phi_2+phi_1*phi_3+phi_1*phi_4), sigma2*(1+phi_1^2+phi_2^2+phi_3^2+phi_4^2)), ncol = 2)

# Tracer l'ellipse de confiance bivariée à 95% pour les prévisions futures (XT+1, XT+2) d'un modèle AR(4)
plot(XT1, XT2, xlim = c(XT1 - 10, XT1 + 10), ylim = c(XT2 - 10, XT2 + 10), xlab = "Prévision de X(T+1)",
     ylab = "Prévision de X(T+2)", main = "Région de confiance bivariée 95%")
lines(ellipse(Sigma, centre = c(XT1, XT2), level = 0.95), type = "l", col = "red", xlab = "Xt+1", ylab = "Xt+2",
      main = "Ellipse de confiance pour (Xt+1, Xt+2)")
abline(h = XT1, v = XT2)











