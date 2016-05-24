###################################################
### Reference:
### Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage 
### randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

require(survival)

###################################################
### code chunk number 2: chunkWRSE
###################################################

#Subfunction to calculate among each arm (A1 or A2)
sub.WRSEestimate <- function(pdata, t) {
  
  #Retrieve data
  n <- nrow(pdata)
  TR <- pdata$TR
  R <- pdata$R
  Z <- pdata$Z # Z=0 for B1, 1 for B2
  U <- pdata$U
  delta <- pdata$delta
  
  #Calculate probability of Z given R
  pi.z <- sum(R*Z) / sum(R)
  
  #Calculate time-varying weight function for each observed Ui
  #W1 for B1; W2 for B2
  #Time as rows, individual as columns
  w1 <- w2 <- matrix(0, nrow=n, ncol=n)
  #Calculate summation of Wj(Ui)Yj(Ui)
  s1 <- s2 <- rep(0, n)
  #Calculate Yi(Uk) for second part of variance sigma
  #Each individual time as one row
  kind <- matrix(0, nrow=n, ncol=n)
  
  for(i in 1:n) {
    
    #Calculate I(TR <= Ui)
    rind <- as.numeric(TR <= U[i])
    
    #Calculate individual weights for each time Ui
    w1[i,] <- 1 - R*rind + R*rind*(1-Z)/(1-pi.z) # weighting for A1B1
    w2[i,] <- 1 - R*rind + R*rind*Z/pi.z # weighting for A1B2
    
    #Calculate Yj(ui)
    yind <- as.numeric(U >= U[i])
    
    #Caculate summation Wj(Ui)Yj(Ui) for each time Ui
    s1[i] <- sum(w1[i,]*yind)
    s2[i] <- sum(w2[i,]*yind)
    
    #Calculate Yi(Uk) for second part of variance sigma
    kind[i,] <- as.numeric(U[i] >= U)
    
  }
  
  #Calculate the estimate and standard error for each time point  
  SURV1 <- SURV2 <- rep(NA, length(t))  
  SE1 <- SE2 <- COV12 <- rep(NA, length(t))
  
  for(j in 1:length(t)) {
    
    ind <- as.numeric(U <= t[j])
    
    #Calculate the survival estimates    
    SURV1[j] <- exp(-sum(diag(w1)[which(s1!=0)]*delta[which(s1!=0)]*ind[which(s1!=0)]/s1[which(s1!=0)]))
    SURV2[j] <- exp(-sum(diag(w2)[which(s2!=0)]*delta[which(s2!=0)]*ind[which(s2!=0)]/s2[which(s2!=0)]))
    
    #Calculate first part (pI) of the variance sigma
    #Calculate second part (pII) of the variance sigma
    pI1 <- pI2 <- rep(0, n)
    pII1 <- pII2 <- rep(0, n)
    sdp1 <- sdp2 <- sdp12 <- rep(0, n)
    
    for(k in 1:n) {
      
      #Calculate first part for each individual
      if(s1[k] != 0) pI1[k] <- w1[k,k]*delta[k]*ind[k]/s1[k]
      if(s2[k] != 0) pI2[k] <- w2[k,k]*delta[k]*ind[k]/s2[k]
      
      #Calculate second part for each individual
      pII1[k] <- sum(w1[which(s1!=0),k]*kind[k,][which(s1!=0)]*diag(w1)[which(s1!=0)]*delta[which(s1!=0)]*ind[which(s1!=0)]/s1[which(s1!=0)]^2)
      pII2[k] <- sum(w2[which(s2!=0),k]*kind[k,][which(s2!=0)]*diag(w2)[which(s2!=0)]*delta[which(s2!=0)]*ind[which(s2!=0)]/s2[which(s2!=0)]^2)
      
      #Calcualte the square of the difference between first part and second part
      sdp1[k] <- (pI1[k] - pII1[k])^2
      sdp2[k] <- (pI2[k] - pII2[k])^2
      sdp12[k] <- (pI1[k] - pII1[k])*(pI2[k] - pII2[k])
      
    }
    
    #Calculate the variance and covariace sigma
    sig1 <- n*sum(sdp1)
    sig2 <- n*sum(sdp2)
    sig12 <- n*sum(sdp12)
    
    #Calculate the variance and covariace
    v1 <- SURV1[j]^2*sig1/n; SE1[j] <- sqrt(v1)
    v2 <- SURV2[j]^2*sig2/n; SE2[j] <- sqrt(v2)
    COV12[j] <- SURV1[j]*SURV2[j]*sig12/n
    
  }
  
  #Return results  
  est <- data.frame(t, SURV1, SURV2, SE1, SE2, COV12)  
  return(est)  
  
}

WRSEestimate <- function(data # A data frame representing the data from sequentially randomized designs
                              # data = data frame {X, TR, R, Z, U, delta}
) {
  
  #Chek for errors
  if (is.null(data$X)) stop("X can not be empty") 
  if (is.null(data$TR)) stop("TR can not be empty")  
  if (is.null(data$R)) stop("R can not be empty")
  if (is.null(data$Z)) stop("Z can not be empty")  
  if (is.null(data$U)) stop("U can not be empty")  
  if (is.null(data$delta)) stop("delta can not be empty") 
  
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
  est1 <- sub.WRSEestimate(pdata=data[which(data$X==0),],t)
  
  #Run sub.LDTestimate to obtain the defined estimates A2 arm
  cat ("Estimating for A2 arm... \n")
  est2 <- sub.WRSEestimate(pdata=data[which(data$X==1),],t)
  
  #Return class
  censorDTR <- censortime <- NULL
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

