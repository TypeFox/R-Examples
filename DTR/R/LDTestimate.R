###################################################
### Reference:
### Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment 
### policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunkfunction
###################################################

#Subfunction to calculate among each arm (A1 or A2)
sub.LDTestimate <- function(pdata, t, L) {
  
  #Retrieve data
  n <- nrow(pdata)
  R <- pdata$R
  Z <- pdata$Z # Z=0 for B1, 1 for B2
  U <- pdata$U
  delta <- pdata$delta
  cens <- 1 - delta
  
  #Calculate probability of Z given R  
  pi.z <- sum(R*Z) / sum(R)
  
  #Calculate weighting Q  
  Q1 <-(1-R) + R*(1-Z)/(1-pi.z) # weighting for A1B1
  Q2 <-(1-R) + R*Z/(pi.z) # weighting for A1B2
  
  #Kaplan-Meier estimate of the censoring survival curve (weighting for censoring) 
  cfit <- summary(survfit(Surv(U,cens)~1))  
  K <- rep(0, n)
  for (i in 1:n) { if(round(U[i],4) < round(min(cfit$time),4)) { K[i] <- 1 } else { dt <- round(cfit$time,4) - round(U[i],4); K[i] <- cfit$surv[which(dt==max(dt[dt<=0]))[1]] } }  
  
  #Combine two forms of inverse probability weighting   
  w1 <- w2 <- rep(0, n)  
  w1[which(K!=0)] <- delta[which(K!=0)] * Q1[which(K!=0)] / K[which(K!=0)] # weighting for A1B1
  w2[which(K!=0)] <- delta[which(K!=0)] * Q2[which(K!=0)] / K[which(K!=0)] # weighting for A1B2
  
  #Calculate the survival estimates irrespective of the arms
  s <- rep(0, n)  
  for (i in 1:n) { sind <- as.numeric(U <= U[i]); s[i] <- 1 - sum(delta[which(K!=0)] * sind[which(K!=0)] / K[which(K!=0)]) / sum(delta[which(K!=0)] / K[which(K!=0)]) }
   
  #Calculate the estimate and standard error for each time point  
  SURV1 <- SURV2 <- rep(NA, length(t))  
  SE1 <- SE2 <- COV12 <- rep(NA, length(t))  
  
  for(j in 1:length(t)) {
    
    ind <- as.numeric(U <= t[j])
    
    #Calculate the survival estimates
    SURV1[j] <- 1 - sum(w1*ind) / sum(w1)
    SURV2[j] <- 1 - sum(w2*ind) / sum(w2)
    
    #Calculate the variance estimates
    #Define G1' and G2'
    G1 <- G2 <- rep(0, n)
    #Define E1 and E12
    E1 <- E2 <- E12 <- rep(0, n)
    #Define Y
    Y <- rep(0, n)
    
    for (k in 1:n) {
      
      pind <- as.numeric(U >= U[k])
      
      #Calculate individual G1' and G2'    
      if (s[k] != 0) G1[k] <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * (ind[which(K!=0)]-1+SURV1[j]) * pind[which(K!=0)] / K[which(K!=0)]) / (n * s[k])
      if (s[k] != 0) G2[k] <- sum(delta[which(K!=0)] * Q2[which(K!=0)] * (ind[which(K!=0)]-1+SURV2[j]) * pind[which(K!=0)] / K[which(K!=0)]) / (n * s[k])		
      
      #Calculate individual E1, E2 and E12
      E1[k] <- sum(delta[K!=0] * (Q1[K!=0]*(ind[which(K!=0)]-1+SURV1[j]) - G1[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
      E2[k] <- sum(delta[K!=0] * (Q2[K!=0]*(ind[which(K!=0)]-1+SURV2[j]) - G2[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
      E12[k] <- sum(delta[K!=0] * (Q1[K!=0]*(ind[which(K!=0)]-1+SURV1[j]) - G1[k]) * (Q2[K!=0]*(ind[which(K!=0)]-1+SURV2[j]) - G2[k]) * pind[which(K!=0)] / K[which(K!=0)]) / n
      
      #Calculate individual Y	
      Y[k] <- sum(pind)
      
    }
    
    #Calcualte variance and covariance    
    v1 <- sum(delta[which(K!=0)] * Q1[which(K!=0)]^2 * (ind[which(K!=0)]-1+SURV1[j])^2 / K[which(K!=0)]) / (n^2) + sum(E1[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; SE1[j] <- sqrt(v1)
    v2 <- sum(delta[which(K!=0)] * Q2[which(K!=0)]^2 * (ind[which(K!=0)]-1+SURV2[j])^2 / K[which(K!=0)]) / (n^2) + sum(E2[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; SE2[j] <- sqrt(v2)		  
    COV12[j] <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * Q2[which(K!=0)] * (ind[which(K!=0)]-1+SURV1[j]) * (ind[which(K!=0)]-1+SURV2[j]) / K[which(K!=0)]) / (n^2) + sum(E12[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n
    
  }
  
  #Return subfunction results
  est <- data.frame(t, SURV1, SURV2, SE1, SE2, COV12)
  return(est)
  
}


LDTestimate <- function(data, # A data frame representing the data from sequentially randomized designs
                              # data = data frame {X, R, Z, U, delta}
                          L=.Machine$double.xmax # Optional restricted survival time L
) {

  #Chek for input errors
  if (is.null(data$X)) stop("X can not be empty") 
  if (is.null(data$R)) stop("R can not be empty")
  if (is.null(data$Z)) stop("Z can not be empty")  
  if (is.null(data$U)) stop("U can not be empty")  
  if (is.null(data$delta)) stop("delta can not be empty")    

  if (is.null(L)) stop("L can not be empty")
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a positive value")
  
  #Times to be assessed
  t <- unique(data$U[which(data$delta==1)])
  #Order times
  t <- t[order(t)]
   
  #Number at risk
  n.risk <- apply(as.array(t), 1, function(x) sum(as.numeric(data$U >= x)))
  
  #Number event
  n.event <- apply(as.array(t), 1, function(x) length(which(data$U==x & data$delta==1)))
  
  #Run sub.LDTestimate to obtain the defined estimates A1 arm
  cat ("Estimating for A1 arm... \n")
  est1 <- sub.LDTestimate(pdata=data[which(data$X==0),],t,L)

  #Run sub.LDTestimate to obtain the defined estimates A2 arm
  cat ("Estimating for A2 arm... \n")
  est2 <- sub.LDTestimate(pdata=data[which(data$X==1),],t,L)
  
  #Return class
  results <- list(Call=match.call(),
                  DTR=c("A1B1", "A1B2", "A2B1", "A2B2"),
                  records=c(length(which((data$X==0 & data$R==0) | (data$X==0 & data$R==1 & data$Z==0))),
                            length(which((data$X==0 & data$R==0) | (data$X==0 & data$R==1 & data$Z==1))),
                            length(which((data$X==1 & data$R==0) | (data$X==1 & data$R==1 & data$Z==0))),
                            length(which((data$X==1 & data$R==0) | (data$X==1 & data$R==1 & data$Z==1)))),
                  events=c(sum(data$delta[which((data$X==0 & data$R==0) | (data$X==0 & data$R==1 & data$Z==0))]), 
                           sum(data$delta[which((data$X==0 & data$R==0) | (data$X==0 & data$R==1 & data$Z==1))]),
                           sum(data$delta[which((data$X==1 & data$R==0) | (data$X==1 & data$R==1 & data$Z==0))]),
                           sum(data$delta[which((data$X==1 & data$R==0) | (data$X==1 & data$R==1 & data$Z==1))])),
                  censorDTR=c(rep("A1B1", length(which((data$X==0 & data$R==1 & data$Z==0 & data$delta==0) | (data$X==0 & data$R==0 & data$delta==0)))),
                              rep("A1B2", length(which((data$X==0 & data$R==1 & data$Z==1 & data$delta==0) | (data$X==0 & data$R==0 & data$delta==0)))),
                              rep("A2B1", length(which((data$X==1 & data$R==1 & data$Z==0 & data$delta==0) | (data$X==1 & data$R==0 & data$delta==0)))),
                              rep("A2B2", length(which((data$X==1 & data$R==1 & data$Z==1 & data$delta==0) | (data$X==1 & data$R==0 & data$delta==0))))),
                  censortime=c(data$U[which((data$X==0 & data$R==1 & data$Z==0 & data$delta==0) | (data$X==0 & data$R==0 & data$delta==0))],
                               data$U[which((data$X==0 & data$R==1 & data$Z==1 & data$delta==0) | (data$X==0 & data$R==0 & data$delta==0))],
                               data$U[which((data$X==1 & data$R==1 & data$Z==0 & data$delta==0) | (data$X==1 & data$R==0 & data$delta==0))],
                               data$U[which((data$X==1 & data$R==1 & data$Z==1 & data$delta==0) | (data$X==1 & data$R==0 & data$delta==0))]),
                  censorsurv=c(apply(as.array(data$U[which((data$X==0 & data$R==1 & data$Z==0 & data$delta==0) | (data$X==0 & data$R==0 & data$delta==0))]), 1, 
                                     function(x) { if(x<min(est1$t)) 1 else est1$SURV1[which(abs(est1$t-x)==min(abs(est1$t-x)[which(est1$t<=x)]))] }),
                               apply(as.array(data$U[which((data$X==0 & data$R==1 & data$Z==1 & data$delta==0) | (data$X==0 & data$R==0 & data$delta==0))]), 1, 
                                     function(x) { if(x<min(est1$t)) 1 else est1$SURV2[which(abs(est1$t-x)==min(abs(est1$t-x)[which(est1$t<=x)]))] }),
                               apply(as.array(data$U[which((data$X==1 & data$R==1 & data$Z==0 & data$delta==0) | (data$X==1 & data$R==0 & data$delta==0))]), 1, 
                                     function(x) { if(x<min(est2$t)) 1 else est2$SURV1[which(abs(est2$t-x)==min(abs(est2$t-x)[which(est2$t<=x)]))] }),
                               apply(as.array(data$U[which((data$X==1 & data$R==1 & data$Z==1 & data$delta==0) | (data$X==1 & data$R==0 & data$delta==0))]), 1, 
                                     function(x) { if(x<min(est2$t)) 1 else est2$SURV2[which(abs(est2$t-x)==min(abs(est2$t-x)[which(est2$t<=x)]))] })),
                  time=t, n.risk=n.risk, n.event=n.event,
                  SURV11=est1$SURV1, SURV12=est1$SURV2, SURV21=est2$SURV1, SURV22=est2$SURV2,
                  SE11=est1$SE1, SE12=est1$SE2, COV1112=est1$COV12,
                  SE21=est2$SE1, SE22=est2$SE2, COV2122=est2$COV12)
  class(results) <- "DTR"
	return(results)

}




