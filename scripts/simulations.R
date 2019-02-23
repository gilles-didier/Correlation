
library(ape)
library(ade4)


VerRes <- "v19"
# alpha 5%
NivAl <- 0.05


###########################################
# Simulation parameters
############################################

nb <- 50000
rhoS <- seq(0,1, by=0.1)
nrho <- length(rhoS)
muAS <- c(0, 0, 0, 5, 1, 3, 3, 5, 5, 10)/10
muBS <- c(0, 2, 5, 0, 2, 3, 5, 5, 10, 10)/10
nmu <- length(muAS)
#
sigma1 <- 3/10
sigma2 <- 3/10
#


###########################################
# Tree used for simulations
############################################


new <- read.delim("Hominin_tree_DMMC.txt", h=F)
new2 <- as.character(new[,1])
conv <- newick2phylog(new2)
arbre <- as.phylo(conv)
nT <- length(arbre$tip.label)
nN <- arbre$Nnode
#
AllTimes <- dist.nodes(arbre)
TL <- AllTimes[(nT+1),1:nT]
names(TL) <- arbre$tip.label
hk <- pic(TL, arbre)




##########################################
# Correlated ABM simulation
############################################"



NomsCol <- c("rho", "muA", "muB", "sigmaA", "sigmaB", "PIC", "PIC_Elliott_avecC", 
             "PIC_Elliott_sansC", "TIPS", "PIC_New", "Tend_YX", "Tend_X", "Tend_Y", 
             "PIC_New_Corr1", "PIC_New_Corr2")
pourc <- matrix(NA, ncol=length(NomsCol), nrow=nmu*nrho)
colnames(pourc) <- NomsCol
#
for (j in 1:nrho)
{
  for (i in 1:nmu)
  {
    rho <- rhoS[j]
    muA <- muAS[i]
    muB <- muBS[i]
    MatV <- matrix(c(sigma1^2, sigma1*sigma2*rho, sigma1*sigma2*rho, sigma2^2), ncol=2)
    #
    mod <- function(x, l) 
    {
      Bt <- rnorm(2, m=0, sd=sqrt(l))
      y1 <- muA*l + sigma1 * Bt[1] 
      y2 <- muB*l + sigma2 * rho * Bt[1] + sigma2 * sqrt(1-rho^2) * Bt[2]
      x+c(y1, y2)
    }
    #
    pvalF <- rep(0,nb)
    pvalK <- rep(0,nb)
    pvalS <- rep(0,nb)
    pvalNX_Y_Xh <- rep(0,nb)
    # 
    for(k in 1:nb)
    {
      x <- rTraitMult(arbre, mod, 2, ancestor = TRUE)
      nf <- Ntip(arbre)
      X <- x[(1:nf),1]
      X <- as.vector(X)
      Y <- x[(1:nf),2]
      Y <- as.vector(Y)
      #
      pic.X <- pic(X, arbre)
      pic.Y <- pic(Y, arbre)
      #
      MUchapX <- sum(pic.X* hk) / sum(hk^2)
      MUchapY <- sum(pic.Y* hk) / sum(hk^2)
      pic.X.mod <- pic.X - MUchapX * hk
      pic.Y.mod <- pic.Y - MUchapY * hk
      #
      resF0 <- lm(pic.Y ~ pic.X - 1)
      pvalF[k] <- summary(resF0)$coef[1,4]
      #
      resK <- lm(pic.Y.mod ~ pic.X.mod)
      pvalK[k] <- summary(resK)$coef[2,4]
      #
      resS <- lm(Y ~ X)
      pvalS[k] <- summary(resS)$coef[2,4]
      #
      resNX <- lm(pic.Y ~ pic.X + hk - 1)
      pvalNX_Y_Xh[k] <- summary(resNX)$coef[1,4]
    }
    #
    nomF <- paste("Pvalsigma", 10*sigma1, "rho", 10*rho, "muA", muA, "muB", muB, VerRes, ".csv", sep="")
    DF <- data.frame(PIC=pvalF, PIC_Elliott_avecC= pvalK, PIC_Elliott_sansC= pvalK0, 
                     TIPS=pvalS, PIC_New=pvalNX_Y_Xh, PIC_New_Corr1=pvalNX_Y_Xh_Corr1, 
                     PIC_New_Corr2=pvalNX_Y_Xh_Corr2 )
    write.csv2(DF, nomF, quote=F, row.names=F)
  }
}


