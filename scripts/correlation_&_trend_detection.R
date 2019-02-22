

library(ape)
library(ade4)
library(lmtest)
library(tseries)
library(caper)



###########################################
# Functions arguments
############################################

# Arguments

# phy 	an object of class "phylo"

# X 	a vector of mode numeric giving the values of a first trait for the tips.
# 	Its values are taken to be in the same order than the tip labels of phy

# Y 	a vector of mode numeric giving the values of a second trait for the tips.
# 	Its values are taken to be in the same order than the tip labels of phy



##########################################
# Trend-detection tests
############################################


SR_trend.test <- function(phy, X)
{
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	TL <- AllTimes[(nT+1),1:nT]
	names(TL) <- phy$tip.label
	resSR <- lm(X ~ TL)
	temp <- summary(resSR)$coef[2,4]
 	JB <-  jarque.bera.test(summary(resSR)$res)$p.val
	DW <- dwtest(resSR, order.by=TL)$p.val
	HMC <- hmctest(resSR, order.by=TL, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("SR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for SR trend-detection test and for conditions checking")
	return(vec)
}


hR_trend.test <- function(phy, X)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	TL <- AllTimes[(nT+1),1:nT]
	names(TL) <- phy$tip.label
	hk <- pic(TL, phy)
	MUchapX <- sum(pic.X* hk) / sum(hk^2)
	reshR <- lm(pic.X ~ hk - 1)
	temp <- c(MUchapX, summary(reshR)$coef[1,4])
	JB <-  jarque.bera.test(summary(reshR)$res)$p.val
	DW <- dwtest(reshR, order.by=hk)$p.val
	HMC <- hmctest(reshR, order.by=hk, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("trend-estimate", "hR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("trend estimation and p-values for hR trend-detection test and for conditions checking")
	return(vec)
}



	






##########################################
# Correlation tests functions 
############################################

SR_corr.test <- function(phy, X, Y)
{
	resSR <- lm(Y ~ X)
	temp <- summary(resSR)$coef[2,4]
	JB <-  jarque.bera.test(summary(resSR)$res)$p.val
	DW <- dwtest(resSR, order.by=X)$p.val
	HMC <- hmctest(resSR, order.by=X, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("SR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for SR correlation test and for conditions checking")
	return(vec)
}


IC_corr.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	resIC <- lm(pic.Y ~ pic.X - 1)
	temp <- summary(resIC)$coef[1,4]
	JB <-  jarque.bera.test(summary(resIC)$res)$p.val
	DW <- dwtest(resIC, order.by=pic.X)$p.val
	HMC <- hmctest(resIC, order.by=pic.X, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("IC", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for IC correlation test and for conditions checking")
	return(vec)
}


DC_corr.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	TL <- AllTimes[(nT+1),1:nT]
	names(TL) <- phy$tip.label
	hk <- pic(TL, phy)
	MUchapX <- sum(pic.X* hk) / sum(hk^2)
	MUchapY <- sum(pic.Y* hk) / sum(hk^2)
	pic.X.mod <- pic.X - MUchapX * hk
	pic.Y.mod <- pic.Y - MUchapY * hk
	resDC <- lm(pic.Y.mod ~ pic.X.mod)
	temp <- summary(resDC)$coef[2,4]
	JB <-  jarque.bera.test(summary(resDC)$res)$p.val
	DW <- dwtest(resDC, order.by=pic.X.mod)$p.val
	HMC <- hmctest(resDC, order.by=pic.X.mod, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("DC", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for DC correlation test and for conditions checking")
	return(vec)
}


MR_corr.test <- function(phy, X, Y)
{
	pic.X <- pic(X, phy)
	pic.Y <- pic(Y, phy)
	nT <- length(phy$tip.label)
	AllTimes <- dist.nodes(phy)
	TL <- AllTimes[(nT+1),1:nT]
	names(TL) <- phy$tip.label
	hk <- pic(TL, phy)
	resMR <- lm(pic.Y ~ pic.X + hk - 1)
	temp <- summary(resMR)$coef[1,4]
	JB <-  jarque.bera.test(summary(resMR)$res)$p.val
	DFuh <- data.frame(UX=pic.X, hk=hk)
	acppc <- dudi.pca(DFuh, center=T, scale=F, scannf = F, nf = ncol(DFuh))
	nbObs <- nrow(DFuh)
	scores <- rep(NA, length=nbObs)
	for (il in 1 :nbObs)
	{
		scores[il] <- sum(acppc$eig*acppc$li[il,]^2)/nbObs
	}
	DW <- dwtest(resMR, order.by=scores)$p.val
	HMC <- hmctest(resMR, order.by=scores, nsim=10000)$p.val
	vec <- c(temp, JB, DW, HMC)
 	names(vec) <- c("MR", "Jarque-Bera", "Durbin-Watson", "Harrison-McCabe")
	print("p-values for MR correlation test and for conditions checking")
	return(vec)
}

	




