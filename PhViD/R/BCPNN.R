`BCPNN` <-
function(DATABASE, RR0 = 1, MIN.n11 = 1, DECISION = 1,
                    DECISION.THRES = 0.05, RANKSTAT = 1, MC = FALSE,
                    NB.MC = 10000) {

DATA <- DATABASE$data
N <- DATABASE$N
L <- DATABASE$L

n11 <- DATA[,1]
n1. <- DATA[,2] 
n.1 <- DATA[,3] 
n10 <- n1. - n11
n01 <- n.1 - n11
n00 <- N - (n11+n10+n01)
E <- n1. * n.1 / N # les counts attendus

if(MIN.n11 > 1) {
  E <- E[n11 >= MIN.n11]
  n1. <- n1.[n11 >= MIN.n11]
  n.1 <- n.1[n11 >= MIN.n11]
  n10 <- n10[n11 >= MIN.n11]
  n01 <- n01[n11 >= MIN.n11]
  n00 <- n00[n11 >= MIN.n11]
#  LL <- data.frame(drugs=L[,1],events=L[,2],n11)
#  LL1 <- LL[,1][n11 >= MIN.n11]
#  LL2 <- LL[,2][n11 >= MIN.n11]
#  rm(list="L")
#  L <- data.frame(LL1,LL2)
  L <- L[n11 >= MIN.n11,]
  n11 <- n11[n11 >= MIN.n11]
}

Nb.Cell <- length(n11)

if (MC == FALSE) {
  post.H0 <- matrix(nrow=Nb.Cell,ncol=length(RR0))
  p1  <- 1 + n1.
  p2  <- 1 + N - n1.
  q1  <- 1 + n.1
  q2  <- 1 + N - n.1
  r1  <- 1 + n11
  r2b <- N - n11 -1 + (2+N)^2/(q1*p1)
  EICb <- log(2)^(-1) * (digamma(r1) - digamma(r1+r2b) - (digamma(p1) - digamma(p1+p2) + digamma(q1) - digamma(q1+q2)))
  VICb <- log(2)^(-2) * (trigamma(r1) - trigamma(r1+r2b) + (trigamma(p1) - trigamma(p1+p2) + trigamma(q1) - trigamma(q1+q2)))
  post.H0 <- pnorm(log(RR0),EICb,sqrt(VICb))
  # Calculation of the Lower Bound
  LB <- qnorm(0.025,EICb,sqrt(VICb))
}

if (MC == TRUE) { # Advanced option MC
  require(MCMCpack)
  n1. <- n11 + n10
  n.1 <- n11 + n01
  Nb_Obs <- length(n11)

  ## Nouvelles priors
  q1. <- (n1. +.5)/(N +1)
  q.1 <- (n.1 +.5)/(N +1)
  q.0 <- (N - n.1 +.5)/(N +1)
  q0. <- (N - n1. +.5)/(N +1)
  
  a.. <- .5/(q1.*q.1) ## le .5 devrait pouvoir être changé
  
  a11 <- q1.*q.1* a..
  a10 <- q1.*q.0* a..
  a01 <- q0.*q.1* a..
  a00 <- q0.*q.0* a..
  
  g11 <- a11 + n11
  g10 <- a10 + n10
  g01 <- a01 + n01
  g00 <- a00 + n00
  g1. <- g11 + g10
  g.1 <- g11 + g01
  
  post.H0 <- vector(length=length(n11))
  LB <- vector(length=length(n11))
  quantile <- vector("numeric",length=length(n11))
  for (m in 1 : length(n11)){
  	p <- rdirichlet(NB.MC,c(g11[m],g10[m],g01[m],g00[m]))
  	p11 <- p[,1]
  	p1. <- p11 + p[,2]
  	p.1 <- p11 + p[,3]	
  	IC_monte <- log(p11/(p1.* p.1))
  	temp <- IC_monte < log(RR0)
  	post.H0[m] <- sum(temp)/NB.MC
  	LB[m] <- sort(IC_monte)[round(NB.MC * 0.025)]
  }
  rm(p11,p1.,p.1,IC_monte,temp)
  gc()
}

if (RANKSTAT==1) RankStat <- post.H0
if (RANKSTAT==2) RankStat <- LB

if (RANKSTAT==1) {
FDR <- (cumsum(post.H0[order(post.H0)]) / (1:length(post.H0)))
FNR <- rev(cumsum((1-post.H0)[order(1-post.H0)])) / (Nb.Cell - 1:length(post.H0))
Se <- cumsum((1-post.H0)[order(post.H0)]) / (sum(1-post.H0))
Sp <- rev(cumsum(post.H0[order(1-post.H0)])) / (Nb.Cell - sum(1-post.H0))
}

if (RANKSTAT==2) {
FDR <- (cumsum(post.H0[order(LB,decreasing=TRUE)]) / (1:length(post.H0)))
FNR <- rev(cumsum((1-post.H0)[order(1-LB,decreasing=TRUE)])) / (Nb.Cell - 1:length(post.H0))
Se <- cumsum((1-post.H0)[order(LB,decreasing=TRUE)]) / (sum(1-post.H0))
Sp <- rev(cumsum(post.H0[order(1-LB,decreasing=TRUE)])) / (Nb.Cell - sum(1-post.H0))
}

# Calculation of the number of signals according to the decision rule (pp/FDR/Nb of Signals)
if (DECISION == 1) Nb.signaux <- sum(FDR <= DECISION.THRES)
if (DECISION == 2) Nb.signaux <- min(DECISION.THRES,Nb.Cell)
if (DECISION == 3) {
    if (RANKSTAT==1) Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT==2) Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
}


############################ SORTIE DE LA FONCTION #############################
RES <- vector(mode="list")

#RES$LIBEL <- L
#colnames(RES$LIBEL) <- c("DRUG","EVENT")

#RES$param.D <- param.D
RES$INPUT.PARAM <- data.frame(RR0,MIN.n11,DECISION,DECISION.THRES,RANKSTAT)
#colnames(RES$INPUT.PARAM) <- c("RR0","MIN.n11","DECISION","DECISION.THRES","RANKSTAT")

# Presentation of the statistics calculated for each couple
##RES$STATISTIC <- data.frame(n11,post.H0,LB)
##rownames(RES$STATISTIC) <- paste(L[,1],L[,2]) # liste des libellés ingénue
##colnames(RES$STATISTIC) <- c("count","postH0","Lower Bound")

# STATISTICAL VALUE TO BE CONSIDERED (used in function compare)
#RES$COMPARE <- vector(mode="list")
#RES$COMPARE$RANKSTAT <- RANKSTAT
#RES$COMPARE$STAT <- RankStat

# SIGNALS RESULTS and presentation
if (RANKSTAT==1) {
  RES$ALLSIGNALS <- data.frame( L[,1][order(RankStat)],
                                L[,2][order(RankStat)],
                                n11[order(RankStat)],
                                E[order(RankStat)],
                                RankStat[order(RankStat)],
                                (n11/E)[order(RankStat)],
                                n1.[order(RankStat)],
                                n.1[order(RankStat)],
                                FDR, FNR, Se, Sp )
  colnames(RES$ALLSIGNALS) <- c("drug code","event effect","count","expected count","post.H0",
                                "n11/E","drug margin","event margin","FDR","FNR","Se","Sp")
}
if (RANKSTAT==2) {
  RES$ALLSIGNALS <- data.frame( L[,1][order(RankStat, decreasing=TRUE)],
                                L[,2][order(RankStat, decreasing=TRUE)],
                                n11[order(RankStat, decreasing=TRUE)],
                                E[order(RankStat,decreasing=TRUE)],
                                RankStat[order(RankStat,decreasing=TRUE)],
                                (n11/E)[order(RankStat,decreasing=TRUE)],
                                n1.[order(RankStat,decreasing=TRUE)],
                                n.1[order(RankStat,decreasing=TRUE)],
                                FDR, FNR, Se, Sp,
                                post.H0[order(RankStat,decreasing=TRUE)] )
  colnames(RES$ALLSIGNALS) <- c("drug code","event effect","count","expected count","Q_0.025(log(IC))",
                                "n11/E","drug margin","event margin","FDR","FNR","Se","Sp","postH0")
}

# List of Signals generated according to the DECISION.THRES
RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux,]

# FDR,FNR,Se,Sp
#RES$OpChar <- data.frame(FDR,FNR,Se,Sp)
#rownames(RES$OpChar) <- paste(RES$ALLSIGNALS[,1],RES$ALLSIGNALS[,2])

# Number of signals
RES$NB.SIGNALS <- Nb.signaux
RES

}

