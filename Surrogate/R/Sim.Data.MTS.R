Sim.Data.MTS <- function(N.Total=2000, N.Trial=50, R.Trial.Target=.8, R.Indiv.Target=.8,  
                     Fixed.Effects=c(0, 0, 0, 0), D.aa=10, D.bb=10, Seed=sample(1:1000, size=1), Model=c("Full")) {   

if ((Model==c("Full") | Model==c("Reduced") | Model==c("SemiReduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\"), Model=c(\"Reduced\"), or Model=c(\"SemiReduced\").")}     

N.patients <- N.Total/N.Trial 
if ((N.patients%%1==0)==FALSE) {
  cat("\nNOTE: The number of patients per trial requested in the function call equals ", N.patients, " (=N.Total/N.Trial), which is not a whole number. ", sep="")
  cat("\nThe number of patients per trial was rounded to ")
  N.patients <- rounded <- ceiling(N.patients)
  N.Total.old <- N.Total
  N.Total <- rounded * N.Trial 
  cat(rounded, " to generate the dataset. Data.Observed.MTS thus contains a total of ", N.Total, " patients \nrather than the requested ", N.Total.old, " in the function call.", sep="")
}

R2.trial.target <- R.Trial.Target**2  

if (Model==c("Full")|Model==c("SemiReduced")){ 
Dmat <- diag(4)
Dmat[3,3] <- D.aa   
Dmat[4,4] <- D.bb 
Dmat[3,4] <- Dmat[4,3] <- sqrt(R2.trial.target * (Dmat[3,3] * Dmat[4,4])) 
set.seed(Seed)
ran.eff <- mvrnorm(N.Trial, rep(0,4), Dmat) 
Smat <- diag(2)
Smat[1,2] <- Smat[2,1] <- R.Indiv.Target 
set.seed(Seed)
errors <- mvrnorm(N.Total, rep(0,2), Smat)

Z <- Trial_ID <- Surr <- True <- NULL
for (i in 1: N.Trial){
  Z_temp <- sample(x=rep(c(-1, 1), each=(ceiling(N.patients/2))), N.patients, replace=FALSE) 
  Trial_ID_temp <- rep(i, N.patients)
  Z <- append(x=Z, values=Z_temp)
  Trial_ID <- append(Trial_ID, Trial_ID_temp)
}
supp <- data.frame(cbind(Z, Trial_ID))

for (i in 1: N.Total){
  Surr_temp <- (Fixed.Effects[1]) + ran.eff[supp$Trial_ID[i],1] + (((Fixed.Effects[3]) + ran.eff[supp$Trial_ID[i],3])*Z[i]) + errors[i,1]
  True_temp <- (Fixed.Effects[2]) + ran.eff[supp$Trial_ID[i],2] + (((Fixed.Effects[4]) + ran.eff[supp$Trial_ID[i],4])*Z[i]) + errors[i,2] 
  Surr <- append(x=Surr, values=Surr_temp)
  True <- append(x=True, values=True_temp)
     }
  }

if (Model==c("Reduced")){ 
  Dmat <- diag(2)
  Dmat[1,1] <- D.aa   
  Dmat[2,2] <- D.bb 
  Dmat[1,2] <- Dmat[2,1] <- sqrt(R2.trial.target * (Dmat[1,1] * Dmat[2,2])) 
  set.seed(Seed)
  ran.eff <- mvrnorm(N.Trial, rep(0,2), Dmat) 
  Smat <- diag(2)
  Smat[1,2] <- Smat[2,1] <- R.Indiv.Target 
  set.seed(Seed)
  errors <- mvrnorm(N.Total, rep(0,2), Smat)
  
  Z <- Trial_ID <- Surr <- True <- NULL
  for (i in 1: N.Trial){
    Z_temp <- sample(x=rep(c(-1, 1), each=(ceiling(N.patients/2))), N.patients, replace=FALSE) 
    Trial_ID_temp <- rep(i, N.patients)
    Z <- append(x=Z, values=Z_temp)
    Trial_ID <- append(Trial_ID, Trial_ID_temp)
  }
  supp <- data.frame(cbind(Z, Trial_ID))
  
  for (i in 1: N.Total){
    Surr_temp <- (((Fixed.Effects[1]) + ran.eff[supp$Trial_ID[i],1])*Z[i]) + errors[i,1]
    True_temp <- (((Fixed.Effects[2]) + ran.eff[supp$Trial_ID[i],2])*Z[i]) + errors[i,2] 
    Surr <- append(x=Surr, values=Surr_temp)
    True <- append(x=True, values=True_temp)
  }
}


Pat_ID <- 1:N.Total
Data.Observed.MTS <- cbind(supp, Surr, True, Pat_ID)
names(Data.Observed.MTS) <- c("Treat", "Trial.ID", "Surr", "True", "Pat.ID")
Data.Observed.MTS <<- Data.Observed.MTS

fit <- list(Data.Observed.MTS=Data.Observed.MTS)

}
