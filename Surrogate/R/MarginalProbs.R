MarginalProbs <- function(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat) {
  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  
  Dataset <- na.exclude(data.frame(cbind(Surr, True, Treat)))
  if (length(unique(Dataset$Treat))!=2) stop("Please make sure that the treatment variable has only 2 levels.")
  if ((sort(unique(Dataset$Treat))[1]==c(-1) & sort(unique(Dataset$Treat))[2]==c(1))==FALSE)
    stop("Please make sure that the treatment is coded as control = -1 and experimental = 1.")
  
  cont <- table(Dataset$Surr[Dataset$Treat==-1], Dataset$True[Dataset$Treat==-1])
  exp <- table(Dataset$Surr[Dataset$Treat==1], Dataset$True[Dataset$Treat==1])
  
  try(Theta_T0S0 <- (cont[1,1]*cont[2,2])/(cont[2,1]*cont[1,2]), silent=TRUE)
  try(Theta_T1S1 <- (exp[1,1]*exp[2,2])/(exp[2,1]*exp[1,2]), silent=TRUE)  
  
  if (exists("Theta_T0S0")==FALSE){Theta_T0S0="NA"}
  if (exists("Theta_T1S1")==FALSE){Theta_T1S1="NA"}
    
  pi1_1_ <- (dim(subset(x=Dataset, subset=(Dataset$True==1 & Dataset$Surr==1 & Dataset$Treat==-1)))[1])/dim(Dataset[Dataset$Treat==-1,])[1] 
  pi0_1_ <- (dim(subset(x=Dataset, subset=(Dataset$True==0 & Dataset$Surr==1 & Dataset$Treat==-1)))[1])/dim(Dataset[Dataset$Treat==-1,])[1]
  pi1_0_ <- (dim(subset(x=Dataset, subset=(Dataset$True==1 & Dataset$Surr==0 & Dataset$Treat==-1)))[1])/dim(Dataset[Dataset$Treat==-1,])[1]
  pi0_0_ <- (dim(subset(x=Dataset, subset=(Dataset$True==0 & Dataset$Surr==0 & Dataset$Treat==-1)))[1])/dim(Dataset[Dataset$Treat==-1,])[1] 
  
  pi_1_1 <- (dim(subset(x=Dataset, subset=(Dataset$True==1 & Dataset$Surr==1 & Dataset$Treat==1)))[1])/dim(Dataset[Dataset$Treat==1,])[1] 
  pi_1_0 <- (dim(subset(x=Dataset, subset=(Dataset$True==1 & Dataset$Surr==0 & Dataset$Treat==1)))[1])/dim(Dataset[Dataset$Treat==1,])[1]
  pi_0_1 <- (dim(subset(x=Dataset, subset=(Dataset$True==0 & Dataset$Surr==1 & Dataset$Treat==1)))[1])/dim(Dataset[Dataset$Treat==1,])[1]
  pi_0_0 <- (dim(subset(x=Dataset, subset=(Dataset$True==0 & Dataset$Surr==0 & Dataset$Treat==1)))[1])/dim(Dataset[Dataset$Treat==1,])[1]
  
  
  fit <- list(Theta_T0S0=Theta_T0S0, Theta_T1S1=Theta_T1S1, Freq.Cont=cont, Freq.Exp=exp, pi1_1_=pi1_1_, pi0_1_=pi0_1_, pi1_0_=pi1_0_, pi0_0_=pi0_0_,
              pi_1_1=pi_1_1, pi_1_0=pi_1_0, pi_0_1=pi_0_1, pi_0_0=pi_0_0)   
  
  class(fit) <- "MarginalProbs"
  fit
    
}
