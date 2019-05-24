


library(ade4)
library(ape)
library(tseries)
library(lmtest)
library(caper)




###########################################
# Functions arguments
############################################

# Arguments

# phy 	an object of class "phylo"

# X 	a vector of mode numeric giving the values of a first trait for the tips.
# 	Its values are taken to be in the same order than the tip labels of phy

# 	Its values are taken to be in the same order than the tip labels of phy



##########################################
# Trend-detection tests
############################################


SR_trend.test <- function(phy, X)
{
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	AT <- AllTimes[(nT+1),1:nT]
	names(AT) <- phy$tip.label
	resTR <- lm(X ~ AT)
	temp <- summary(resTR)$coef[2,4]
 	JB <-  jarque.bera.test(summary(resTR)$res)$p.val
	DW <- dwtest(resTR, order.by=TL)$p.val
	HMC <- hmctest(resTR, order.by=TL, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("SR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for SR trend-detection test with conditions checking")
	return(vec)
}


hR_trend.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	AT <- AllTimes[(nT+1),1:nT]
	names(AT) <- phy$tip.label
	hk <- pic(AT, phy)
	MUchapX <- sum(pic.X* hk) / sum(hk^2)
	reshR <- lm(pic.X ~ hk - 1)
	temp <- c(MUchapX, summary(reshR)$coef[1,4])
	JB <-  jarque.bera.test(summary(reshR)$res)$p.val
	DW <- dwtest(reshR, order.by=hk)$p.val
	HMC <- hmctest(reshR, order.by=hk, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("trend-estimate", "hR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("trend estimation and p-values for hR trend-detection test with conditions checking")
	return(vec)
}



	






##########################################
# Correlation tests functions 
############################################

SR_corr.test <- function(phy, X, Y)
{
	resTR <- lm(Y ~ X)
	temp <- summary(resTR)$coef[2,4]
	JB <-  jarque.bera.test(summary(resTR)$res)$p.val
	DW <- dwtest(resTR, order.by=X)$p.val
	HMC <- hmctest(resTR, order.by=X, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("SR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for SR correlation test with conditions checking")
	return(vec)
}


IC_corr.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	resCR <- lm(pic.Y ~ pic.X - 1)
	temp <- summary(resCR)$coef[1,4]
	JB <-  jarque.bera.test(summary(resCR)$res)$p.val
	DW <- dwtest(resCR, order.by=pic.X)$p.val
	HMC <- hmctest(resCR, order.by=pic.X, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("IC", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for IC correlation test with conditions checking")
	return(vec)
}


PGLS_corr.test <- function(phy, X, Y)
{
  DF <-  data.frame(Species=phy$tip.label, X=X, Y=Y)
  rownames(DF) <- phy$tip.label
  nT <- length(phy$tip.label)
  Global <- comparative.data(phy, DF, "Species")
  resPR <- pgls(Y ~ X, data=Global)
  temp <- summary(resPR)$coef[2,4]
  JB <-  jarque.bera.test(summary(resPR)$res)$p.val
  DW <- dwtest(resPR, order.by=X)$p.val
  HMC <- hmctest(resPR, order.by=X, nsim=10000)$p.val
  vec <- c(temp, JB, DW, HMC)
  names(vec) <- c("PGLS", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
  print("p-values for PR correlation test with conditions checking")
  return(vec)
}  

DR_corr.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	AT <- AllTimes[(nT+1),1:nT]
	names(AT) <- phy$tip.label
	hk <- pic(AT, phy)
	MUchapX <- sum(pic.X* hk) / sum(hk^2)
	MUchapY <- sum(pic.Y* hk) / sum(hk^2)
	pic.X.mod <- pic.X - MUchapX * hk
	pic.Y.mod <- pic.Y - MUchapY * hk
	resDCR <- lm(pic.Y.mod ~ pic.X.mod)
	temp <- summary(resDCR)$coef[2,4]
	JB <-  jarque.bera.test(summary(resDCR)$res)$p.val
	DW <- dwtest(resDCR, order.by=pic.X.mod)$p.val
	HMC <- hmctest(resDCR, order.by=pic.X.mod, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("DR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for DR correlation test with conditions checking")
	return(vec)
}



MR_corr.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	AT <- AllTimes[(nT+1),1:nT]
	names(AT) <- phy$tip.label
	hk <- pic(AT, phy)
	resMCR <- lm(pic.Y ~ pic.X + hk - 1)
	temp <- summary(resMCR)$coef[1,4]
	JB <-  jarque.bera.test(summary(resMCR)$res)$p.val
	DFuh <- data.frame(UX=pic.X, hk=hk)
	acppc <- dudi.pca(DFuh, center=T, scale=F, scannf=F, nf=ncol(DFuh))
	nbObs <- nrow(DFuh)
	scores <- rep(NA, length=nbObs)
	for (il in 1 :nbObs)
	{
		scores[il] <- sum(acppc$eig*acppc$li[il,]^2)/nbObs
	}
	DW <- dwtest(resMCR, order.by=scores)$p.val
	HMC <- hmctest(resMCR, order.by=scores, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("MR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for MR correlation test with conditions checking")
	return(vec)
}


PGLSt_corr.test <- function(phy, X, Y)
{
  DF <-  data.frame(Species=phy$tip.label, X=X, Y=Y)
  rownames(DF) <- phy$tip.label
  nT <- length(phy$tip.label)
  Global <- comparative.data(phy, DF, "Species")
  AllTimes <- dist.nodes(phy)
  AT <- AllTimes[(nT+1),1:nT]
  resMPR <- pgls(Y ~ X + AT, data=Global)
  temp <- summary(resMPR)$coef[2,4]
  JB <-  jarque.bera.test(summary(resMPR)$res)$p.val
  DFxt <- data.frame(X=X, Time=AT)
  acppc <- dudi.pca(DFxt, center=T, scale=F, scannf=F, nf=ncol(DFxt))
  nbObs <- nrow(DFxt)
  scores <- rep(NA, length=nbObs)
  for (il in 1 :nbObs)
  {
    scores[il] <- sum(acppc$eig*acppc$li[il,]^2)/nbObs
  }
  DW <- dwtest(resMPR, order.by=scores)$p.val
  HMC <- hmctest(resMPR, order.by=scores, nsim=10000)$p.val
  vec <- c(temp, JB, DW, HMC)
  names(vec) <- c("PGLSt", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
  print("p-values for PGLSt correlation test with conditions checking")
  return(vec)
}  	




