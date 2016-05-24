####################################################################################
#THEORETICAL
####################################################################################

HBH <- function(p1,p2,L,N)
{
  if(p1 >= 1 | p1 <= 0) stop("p1 must be a real number in (0,1)")
  if(p2 >= 1 | p2 <= 0) stop("p2 must be a real number in (0,1)")
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")

  cat("   ","\n")
  cat("###########################################","\n")
  cat("HOP BY HOP - THEORETICAL RESULTS","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Data success probability         p1 = ",p1),"\n")
  cat(paste("ACK success probability          p2 = ",p2),"\n")
  cat(paste("Maximum number of transmissions  L  = ",L),"\n")
  cat(paste("Number of Hops                   N  = ",N),"\n")
  cat("   ","\n")

  P = function(y,p1,p2,L)
  {
    pp = p1*p2
    return(pp*(1-pp)^(y-1) + ifelse(y==L,(1-pp)^L,0))
  }

  if(L==Inf)
  {
    expect1 = 1/(p1*p2)
    expect2 = p1*expect1
    ET1 = (1+p1)/(p1*p2)

    REC_expect1 = 1/p2
    REC_expect2 = 1
    REC_ET1 = (1+p2)/p2
  }else
  {
    y = seq(1,L)
    expect1 = sum(y*sapply(X = y,FUN = P,p1,p2,L))
    expect2 = p1*expect1
    ET1  = expect1 + expect2

    REC_expect0 = sum(y*sapply(X = y,FUN = P,p1,p2,L))
    REC_expect1 = p1*REC_expect0
    REC_expect2 = p2*REC_expect1
    REC_ET1  = REC_expect1 + REC_expect2
  }

  Pr1    = 1-(1-p1)^L
  PrS    = Pr1^N
  PrSV   = Pr1^seq(1,N)
  w      = Pr1^seq(0,N-1)
  geo    = sum(w)

  ETData = expect1*geo
  ETACK  = expect2*geo
  ETS    = ET1*geo

  ETDataV = expect1*w
  ETACKV  = expect2*w
  ETSV    = ET1*w

  REC_ETData = REC_expect1*geo
  REC_ETACK  = REC_expect2*geo
  REC_ETS    = REC_ET1*geo

  REC_ETDataV = REC_expect1*w
  REC_ETACKV  = REC_expect2*w
  REC_ETSV    = REC_ET1*w

  res    = round(matrix(data = c(c(PrSV,PrS),c(ETDataV,ETData),c(ETACKV,ETACK),c(ETSV,ETS),c(REC_ETDataV,REC_ETData),c(REC_ETACKV,REC_ETACK),c(REC_ETSV,REC_ETS)),nrow = 7,ncol = N+1,byrow = T),3)
  rownames(res)=c("Success Probability","Expected Data Transmissions","Expected ACK Transmissions","Expected Total Transmissions","Expected Data Receptions","Expected ACK Receptions","Expected Total Receptions")
  colnames(res)= c(paste("Hop ",seq(1,N),"/",N,sep = ""),"Total")
  return(res)
}

####################################################################################
#MONTE CARLO
####################################################################################

MCHBH = function(p1,p2,L,N,M=5000)
{
  if(p1 >= 1 | p1 <= 0) stop("p1 must be a real number in (0,1)")
  if(p2 >= 1 | p2 <= 0) stop("p2 must be a real number in (0,1)")
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")
  if(M%%1!=0 | M<1) stop("N must be a positive integer")

  cat("   ","\n")
  cat("###########################################","\n")
  cat("HOP BY HOP - MONTE CARLO SIMULATION RESULTS","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Data success probability         p1 = ",p1),"\n")
  cat(paste("ACK success probability          p2 = ",p2),"\n")
  cat(paste("Maximum number of transmissions  L  = ",L),"\n")
  cat(paste("Number of Hops                   N  = ",N),"\n")
  cat(paste("Monte Carlo Simulations          M  = ",M),"\n")
  cat("   ","\n")

  prog2 = function(p1,p2,L,N)
  {
    pos=1
    transD = 0
    transA = 0
    RtransD = 0
    RtransA = 0
    paso = TRUE
    TX = matrix(data = 0,nrow = N,ncol = 3)
    RX = matrix(data = 0,nrow = N,ncol = 3)

    while(pos<=N && paso == TRUE)
    {
      trial = 1
      paso = FALSE
      resp = FALSE
      while(trial<=L && resp == FALSE)
      {
        transD = transD +1
        TX[pos,1] = TX[pos,1] + 1
        u1 = runif(1)
        if(u1<p1)
        {
          paso = TRUE
          RtransD = RtransD +1
          RX[pos,1] = RX[pos,1] + 1
          transA = transA +1
          TX[pos,2] = TX[pos,2] + 1
          u2 = runif(1)
          if(u2<p2)
          {
            RtransA = RtransA +1
            RX[pos,2] = RX[pos,2] + 1
            resp=TRUE
          }
        }
        trial=trial+1
      }
      if(paso==TRUE){pos = pos + 1}
    }
    TX[,3] = TX[,1] + TX[,2]
    RX[,3] = RX[,1] + RX[,2]
    return(list(pos=pos,tDATA = transD,tACK = transA,trans = transD+transA,TX = TX,RDATA = RtransD,RACK = RtransA,Rtrans = RtransD+RtransA,RX = RX))
    cat("   ","\n")
  }

  resul = matrix(data = NA,nrow = M,ncol = 7)
  resulA = resulB = array(data = NA,dim = c(M,N,3))
  resulB = array(data = NA,dim = c(M,N,3))
  pb <- txtProgressBar(min = 0, max = M, style = 3)

  for(k in 1:M)
  {
    run = prog2(p1,p2,L,N)
    resul[k,1]=run$pos
    resul[k,2]=run$tDATA
    resul[k,3]=run$tACK
    resul[k,4]=run$trans
    resul[k,5]=run$RDATA
    resul[k,6]=run$RACK
    resul[k,7]=run$Rtrans
    resulA[k,,] = run$TX
    resulB[k,,] = run$RX
    setTxtProgressBar(pb, k)
  }
  close(pb)

  counts = matrix(data = 0,nrow = N+1,ncol = 1)
  counts[sort(unique(resul[,1]))] = c(as.data.frame(table(resul[,1]))$Freq)

  success = rev(cumsum(rev(counts)))[-1]/M
  mean2 = apply(X = resulA,MARGIN = 2:3,FUN = mean)
  mean = apply(X = resul[,2:4],MARGIN = 2,FUN = mean)
  Rmean2 = apply(X = resulB,MARGIN = 2:3,FUN = mean)
  Rmean = apply(X = resul[,5:7],MARGIN = 2,FUN = mean)
  success[is.na(success)] = 1

  table = round(rbind(c(success,success[N]),cbind(t(mean2),mean),cbind(t(Rmean2),Rmean)),3)
  rownames(table)=c("MC Success Probability","MC Mean Data Transmissions","MC Mean ACK Transmissions","MC Mean Total Transmissions","MC Mean Data Receptions","MC Mean ACK Receptions","MC Mean Total Receptions")
  colnames(table)= c(paste("Hop ",seq(1,N),"/",N,sep = ""),"Total")
  return(table)
  cat("   ","\n")
}

# #RUN
#
# #Theoretical
# HBH(p1=0.65,p2=0.4,L=7,N=5)
# #MonteCarlo Simulations
# MCHBH(p1=0.65,p2=0.4,L=7,N=5,M=1000)
