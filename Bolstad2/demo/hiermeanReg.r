priorTau <- list(tau0 = 0, v0 = 1000)
priorPsi <- list(psi0 = 500, eta0 = 1)
priorVar <- list(s0 = 500, kappa0 = 1)
priorBeta <- list(b0 = c(0,0), bMat = matrix(c(1000,100,100,1000), nc = 2))

data(hiermeanRegTest.df)
data.df <- hiermeanRegTest.df
design <- list(y = data.df$y, group = data.df$group,
               x = as.matrix(data.df[,3:4]))

r<-hierMeanReg(design, priorTau, priorPsi, priorVar, priorBeta)
