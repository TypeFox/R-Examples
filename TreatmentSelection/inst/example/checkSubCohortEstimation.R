library(TreatmentSelection)

n = 1e4
expit = function(x) exp(x)/(1+exp(x))

bootstraps = 0

if(marker.type ==1){ 
  marker.type <- "weak"
  #alpha.onc <- c(-2.236,.4, .037, -0.0192)
  a<-c(-1.23, -0.09, .6, -0.5)  
  mu = 0; sd = 1
  square = FALSE
  
  
}else if(marker.type ==0){
  marker.type <- "strong"
  #alpha.strong <- c( -1.2402598, -0.6910426,  0.6, -2.25) 
  a <-  c( -1.24, -0.29,  0.6, -1.5)
  mymean <- 0; mysd <- 1
  square <- FALSE
}




a0 <- a[1] 
a1 <- a[2]
a2 <- a[3]
a3 <- a[4]

Sims = 100
VD <- data.frame("cohort" = rep(NA, Sims), 
                 "casecontrol" = rep(NA, Sims), 
                 "strat" = rep(NA, Sims))
for(s in 1:Sims){
  Y <- rnorm(n, mean = mu, sd = sd)
  if(square) Y <- Y^2 
  trt <- rbinom(n,size=1,prob=0.5) #randomly assign treatment/non treatment
  D   <- numeric(n)
  D[trt == 0] <- rbinom(sum(trt==0),size=1,prob=expit(a0      + a2*Y[trt==0]))
  D[trt==1]   <- rbinom(sum(trt==1),size=1,prob=expit(a0 + a1 + a2*Y[trt==1] + a3*Y[trt==1]))
  
 
cohort.dat <- data.frame("event" = D, "trt" = trt, "Y" = Y)
mytrtsel <- trtsel( event = "event", trt = "trt", marker = "Y",data = cohort.dat,study.design = "randomized cohort")

out.cohort <- eval.trtsel(mytrtsel, bootstraps = bootstraps)




  ### 1:2 Case Control
  
  nmatch = 2
  S <- NULL
  S[D==1] <- 1 #select all cases
  numcontrols <- sum(D==1)*nmatch
  S[D==0] <- sample(rep(c(1,0), c(numcontrols, sum(D==0) - numcontrols) )) #sample(c(rep(1,numcontrols),rep(0,sum(D==0)-numcontrols)))
  
  tmpD<-D[S==1]; tmptrt<-trt[S==1]; tmpY<-Y[S==1]
  tmpData <- data.frame("D" = tmpD, "T" = tmptrt, "Y" = tmpY)
  
  tmptrtsel <- trtsel( event = "D", trt = "T", marker = "Y", data = tmpData,
                       study.design = "nested case-control",
                       cohort.attributes = c(n, mean(trt), mean(D), 1))
  
  
  out.cc2 <- eval.trtsel(tmptrtsel, bootstraps = bootstraps)
 

  ### 1:2 stratified case-control
  nmatch = 2
  S <- NULL
  S[D==1] <- 1 #select all cases
  numcontrols <- min(sum(D==1 & trt==0)*nmatch, sum(D==0 & trt==0)) #number of cases in trt==0
  S[D==0 & trt==0] <- sample( rep(c(1,0), c(numcontrols, sum(D==0 & trt==0)-numcontrols)))
  
  numcontrols <- min(sum(D==1 & trt==1)*nmatch, sum(D==0 & trt==1)) #number of cases in trt==1
  S[D==0 & trt==1] <- sample( rep(c(1,0), c(numcontrols, sum(D==0 & trt==1)-numcontrols)))
  
  tmpD<-D[S==1]; tmptrt<-trt[S==1]; tmpY<-Y[S==1]
  
 
tmpData <- data.frame("D" = tmpD, "T" = tmptrt, "Y" = tmpY)

  tmptrtsel <- trtsel( event = "D", trt = "T", marker = "Y", data = tmpData, 
                       study.design = "stratified nested case-control",
                       cohort.attributes = c(n, mean(trt==0 & D==0), 
                                             mean(trt==0 & D==1), 
                                             mean(trt==1 & D==0), 1, 1))
  out.scc2 <- eval.trtsel(tmptrtsel, bootstraps = bootstraps)


#round(rbind(out.cohort$estimates, out.cc2$estimates, out.scc2$estimates), 4)

VD[s,] <- c(out.cohort$estimates$Var.Delta, 
            out.cc2$estimates$Var.Delta, 
            out.scc2$estimates$Var.Delta)
print(s)
print(round(VD[s,], 4))

}


colMeans(VD)
VD.long <- melt(VD)
ggplot(VD.long, aes(value, color = variable, fill = variable)) + geom_histogram(alpha = .2) + facet_grid(variable~.)


