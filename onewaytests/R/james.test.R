james.test <- function(y, group, alpha=0.05, na.rm = TRUE) {
  
if (na.rm){
    completeObs <- complete.cases(y, group)
    y <- y[completeObs]
    group <- group[completeObs]
  }
  df <- data.frame(Response = y, Group = group)
  DNAME <- "y vs group"
  METHOD <- "James's Second-Order Test"

  n <- length(y)
  x.levels <- levels(factor(group))
  J<-length(x.levels)
  c<-qchisq(1-alpha, J-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  y.vars <- y.means <- y.n <- y.standarderror <- a <- Ybar <- t <- T <-NULL
  R10 <- R11 <- R12 <- R20 <- R21 <- R22 <- R23 <- CV<- NULL
  
 for (i in x.levels) {
    
    y.vars[i] <- var(y[group==i])
    
    y.means[i] <- mean(y[group==i])
    
    y.n[i] <- length(y[group==i])
    
    y.standarderror[i] <- sqrt(y.vars[i]/y.n[i])
    
  }
  
  for (j in x.levels) {
    
    a[j] <- (1/(y.standarderror[j])^2)/(sum(1/(y.standarderror)^2))
    
    Ybar[j] <-a[j]*y.means[j]
    
    T[j]<- ((1-a[j])^2)/(y.n[j]-1)
    
    R10[j]<-(a[j]^0)/((y.n[j]-1)^1)
    R11[j]<-(a[j]^1)/((y.n[j]-1)^1)
    R12[j]<-(a[j]^2)/((y.n[j]-1)^1)
    R20[j]<-(a[j]^0)/((y.n[j]-1)^2)
    R21[j]<-(a[j]^1)/((y.n[j]-1)^2)
    R22[j]<-(a[j]^2)/((y.n[j]-1)^2)
    R23[j]<-(a[j]^3)/((y.n[j]-1)^2)
  }
   
  R10<-sum(R10)
  R11<-sum(R11)
  R12<-sum(R12)
  R20<-sum(R20)
  R21<-sum(R21)
  R22<-sum(R22)
  R23<-sum(R23)
  
Tsum<-sum(T)
Ybarsum<-sum(Ybar)
  
  for (k in x.levels) {
        
    t[k] <- (y.means[k]-Ybarsum)/y.standarderror[k]
    
  }

Jtest=sum(t^2)
  
####### critical value calculation ####

chi2<- (c^1)/(J+2-3)
chi4<- (c^2)/((J+2-3)*(J+4-3))
chi6<- (c^3)/((J+2-3)*(J+4-3)*(J+6-3))
chi8<- (c^4)/((J+2-3)*(J+4-3)*(J+6-3)*(J+8-3))

CV <- c+((1/2)*(3*chi4+chi2)*Tsum)+((1/16)*((3*chi4+chi2)^2)*(1-((J-3)/c))*(Tsum^2))+((1/2)*(3*chi4+chi2))*((8*R23-10*R22+4*R21-6*R12^2+8*R12*R11-4*R11^2)+(2*R23-4*R22+2*R21-2*R12^2+4*R12*R11-2*R11^2)*(chi2-1)+(1/4)*(-(R12^2)+4*R12*R11-2*R12*R10-4*R11^2+4*R11*R10-R10^2)*(3*chi4-2*chi2-1))+(R23-3*R22+3*R21-R20)*(5*chi6+2*chi4+chi2)+(3/16)*(R12^2-4*R23+6*R22-4*R21+R20)*(35*chi8+15*chi6+9*chi4+5*chi2)
+(1/16)*(-2*R22+4*R21-R20+2*R12*R10-4*R11*R10+R10^2)*(9*chi8-3*chi6-5*chi4-chi2)+(1/4)*(-R22+R11^2)*(27*chi8+3*chi6+chi4+chi2)+(1/4)*(R23-R12*R11)*(45*chi8+9*chi6+7*chi4+3*chi2)


result<-( if (Jtest >= CV) "Reject H_0" else "Fail to reject H_0")


names(Jtest) <- "Jtest"
names(CV) <- c("CriticalValue")
  
 if(sum(!completeObs) > 0){
if(sum(!completeObs)==1){cat("\n", paste("NOTE: ", sum(!completeObs), " of ", sum(!completeObs)+length(y), "observations was removed due to missingness."), "\n")
}else  {cat("\n", paste("NOTE: ", sum(!completeObs), " of ", sum(!completeObs)+length(y), "observations were removed due to missingness."), "\n")}
}

structure(list(statistic = Jtest, parameter= CV, method = METHOD, data.name = DNAME), class = "htest")

}

