Sim.Data.STS <- function(N.Total=2000, R.Indiv.Target=.8, Means=c(0, 0, 0, 0), Seed=sample(1:1000, size=1)) { 

Smat <- diag(2)
Smat[1,2] <- Smat[2,1] <- R.Indiv.Target 
set.seed(Seed)
errors <- mvrnorm(N.Total, rep(0,2), Smat)
Treat <- sample(x=rep(c(-1, 1), each=(ceiling(N.Total/2))), N.Total, replace=FALSE)
Pat_ID <- 1:N.Total
supp <- data.frame(cbind(errors, Treat, Pat_ID))
colnames(supp) <- c("Surr", "True", "Treat", "Pat.ID")
supp$Surr[supp$Treat==-1] <- supp$Surr[supp$Treat==-1]+Means[1] 
supp$Surr[supp$Treat==1] <- supp$Surr[supp$Treat==1]+Means[2] 
supp$True[supp$Treat==-1] <- supp$True[supp$Treat==-1]+Means[3] 
supp$True[supp$Treat==1] <- supp$True[supp$Treat==1]+Means[4] 

Data.Observed.STS <- NULL
Data.Observed.STS <<- supp

fit <- list(Data.Observed.STS=Data.Observed.STS)

}


