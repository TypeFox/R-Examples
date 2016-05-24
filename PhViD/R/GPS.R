`GPS` <-
function(DATABASE, RR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES = 0.05, RANKSTAT = 1, TRONC = FALSE, TRONC.THRES = 1, PRIOR.INIT = c(alpha1= 0.2, beta1= 0.06, alpha2=1.4, beta2=1.8, w=0.1), PRIOR.PARAM = NULL){

DATA <- DATABASE$data
N <- DATABASE$N
L <- DATABASE$L          

n11 <- DATA[,1] # les nij ou n11 (c'est pareil)
n1. <- DATA[,2] # les marges lignes (effets indésirables)
n.1 <- DATA[,3] # les marges colonnes (médicaments)
E <- DATA[,2] * DATA[,3] / N # les effectifs attendus      

P_OUT <- TRUE

if (is.null(PRIOR.PARAM)) {    
  P_OUT <- FALSE   
  f<-function(p,n11,E){
	sum(-log((p[5] * dnbinom(n11, size=p[1], prob=p[2]/(p[2]+E)) + (1-p[5]) * dnbinom(n11, size=p[3], prob=p[4]/(p[4]+E)))))
  }

  if (TRONC == FALSE){  
    # Dumouchel situation (we will consider the whole contingency table with 0)
    data_cont<-xtabs(DATA[,1]~L[,1]+L[,2])
    n1._mat <- apply(data_cont,1, sum) # on recalcule les marges des lignes...
    n.1_mat <- apply(data_cont,2, sum) # ...et des colonnes
    n1._c <- rep(n1._mat, times=length(n.1_mat))
    n.1_c <- rep(n.1_mat, each=length(n1._mat))
    E_c <- n1._c * n.1_c / N
    n11_c <- as.vector(data_cont)

    p_out <-nlm(f, p=PRIOR.INIT, n11=n11_c, E=E_c, iterlim=500)
  }
  
  # alternative tenant compte de la troncature
  if (TRONC == TRUE){
    tronc <- TRONC.THRES - 1  # si l'utilisateur a décidé de tronquer sa base
    lik <-function(p,n11,E,tronc){
      sum(-log(
        (p[5] * dnbinom(n11, size=p[1], prob=p[2]/(p[2]+E))
      + (1-p[5]) * dnbinom(n11, size=p[3], prob=p[4]/(p[4]+E)))
      / (1-(p[5] * pnbinom(tronc, size=p[1], prob=p[2]/(p[2]+E))
      + (1-p[5]) * pnbinom(tronc, size=p[3], prob=p[4]/(p[4]+E))))
      ))
    }
  
  p_out <-nlm(lik, p=PRIOR.INIT, n11=n11[n11>=TRONC.THRES], E=E[n11>=TRONC.THRES], tronc,iterlim=500)
  }
  
  PRIOR.PARAM <-p_out$estimate
  code.convergence <- p_out$code
}

# Algorithm allowing to conserve only the couples presenting the minimal nb of notifications entered by the user
if(MIN.n11 > 1) {
  E <- E[n11 >= MIN.n11]
  n1. <- n1.[n11 >= MIN.n11]
  n.1 <- n.1[n11 >= MIN.n11]
  LL <- data.frame(drugs=L[,1],events=L[,2],n11)
  LL1 <- LL[,1][n11 >= MIN.n11]
  LL2 <- LL[,2][n11 >= MIN.n11]
  rm(list="L")
  L <- data.frame(LL1,LL2)
  n11 <- n11[n11 >= MIN.n11]
}

Nb.Cell <- length(n11)
post.H0 <- vector(length=Nb.Cell)


# Posterior probability of the null hypothesis
Q <- PRIOR.PARAM[5]*dnbinom(n11,size=PRIOR.PARAM[1],prob=PRIOR.PARAM[2]/(PRIOR.PARAM[2]+E)) / 
    (PRIOR.PARAM[5]*dnbinom(n11,size=PRIOR.PARAM[1],prob=PRIOR.PARAM[2]/(PRIOR.PARAM[2]+E)) + (1-PRIOR.PARAM[5])* 
     dnbinom(n11,size=PRIOR.PARAM[3],prob=PRIOR.PARAM[4]/(PRIOR.PARAM[4]+E)))

post.H0 <- Q*pgamma(RR0,PRIOR.PARAM[1]+n11,PRIOR.PARAM[2]+E) +(1-Q)*pgamma(RR0,PRIOR.PARAM[3]+n11,PRIOR.PARAM[4]+E) # proba a posteriori de H0

# Posterior Expectation of Lambda
postE <- log(2)^(-1)*(Q*(digamma(PRIOR.PARAM[1]+n11)-log(PRIOR.PARAM[2]+E))+(1-Q)*(digamma(PRIOR.PARAM[3]+n11)-log(PRIOR.PARAM[4]+E)))


# Algorithme allowing the calculation of the CI Lower Bound (at Seuil%)

# Algorithm Emmanuel Roux
	QuantileDuMouchel<-function(Seuil,Q,a1,b1,a2,b2) {
		m<-rep(-100000,length(Q))
		M<-rep(100000,length(Q))
		x<-rep(1,length(Q))
		Cout<-FCoutQuantileDuMouchel(x,Seuil,Q,a1,b1,a2,b2)
		while (max(round(Cout*1e4))!=0)	{
			S<-sign(Cout)
			xnew<-(1+S)/2*((x+m)/2)+(1-S)/2*((M+x)/2)
			M<-(1+S)/2*x+(1-S)/2*M
			m<-(1+S)/2*m+(1-S)/2*x
			x<-xnew
			Cout<-FCoutQuantileDuMouchel(x,Seuil,Q,a1,b1,a2,b2)
		}
		x
	}
	FCoutQuantileDuMouchel<-function(p,Seuil,Q,a1,b1,a2,b2) {
		Q*pgamma(p,shape=a1,rate=b1)+(1-Q)*pgamma(p,shape=a2,rate=b2)-Seuil
	}

# Calculation of the Lower Bound. (VALUEbis is former LBQ)
LB <- QuantileDuMouchel(0.05,Q,PRIOR.PARAM[1]+n11,PRIOR.PARAM[2]+E,PRIOR.PARAM[3]+n11,PRIOR.PARAM[4]+E)


# Assignment based on the method (pp/postE/LB)
if (RANKSTAT==1) RankStat <- post.H0
if (RANKSTAT==2) RankStat <- LB
if (RANKSTAT==3) RankStat <- postE

# FDR, FNR, Se and Sp based on the method (pp/postE/LB)
if (RANKSTAT==1) {
  FDR <- (cumsum(post.H0[order(RankStat)]) / (1:length(post.H0)))
  FNR <- rev(cumsum((1-post.H0)[order(1-RankStat)])) / (Nb.Cell - 1:length(post.H0))
  Se <- cumsum((1-post.H0)[order(RankStat)]) / (sum(1-post.H0))
  Sp <- rev(cumsum(post.H0[order(1-RankStat)])) / (Nb.Cell - sum(1-post.H0))
}
if (RANKSTAT==2 | RANKSTAT==3) {
  FDR <- (cumsum(post.H0[order(RankStat,decreasing=TRUE)]) / (1:length(post.H0)))
  FNR <- rev(cumsum((1-post.H0)[order(1-RankStat,decreasing=TRUE)])) / (Nb.Cell - 1:length(post.H0))
  Se <- cumsum((1-post.H0)[order(RankStat,decreasing=TRUE)]) / (sum(1-post.H0))
  Sp <- rev(cumsum(post.H0[order(1-RankStat,decreasing=TRUE)])) / (Nb.Cell - sum(1-post.H0))
}

# Number of signals according to the decision rule (pp/FDR/Nb of Signals)

if (DECISION == 1) Nb.signaux <- sum(FDR <= DECISION.THRES)
if (DECISION == 2) Nb.signaux <- min(DECISION.THRES,Nb.Cell)
if (DECISION == 3) {
    if (RANKSTAT==1) Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT==2 | RANKSTAT==3) Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
} 


############################ Output #############################
RES <- vector(mode="list")

#RES$LIBEL <- L
#colnames(RES$LIBEL) <- c("DRUG","EVENT")

# list of the parameters used
RES$INPUT.PARAM <- data.frame(RR0,MIN.n11,DECISION,DECISION.THRES,RANKSTAT,TRONC,TRONC.THRES)
#colnames(RES$INPUT.PARAM) <- c("RR0","MIN.n11","TRONC","TRONC.THRES","DECISION","DECISION.THRES","RANKSTAT")

RES$PARAM <- vector(mode="list")
# vector of the final a priori parameters (if P_OUT=TRUE)
if (P_OUT==TRUE) RES$PARAM$PRIOR.PARAM <- data.frame(PRIOR.PARAM)
# vector of the initial a priori and final a priori parameters (if P_OUT=FALSE)
if (P_OUT==FALSE) {
  RES$PARAM$PRIOR.INIT <- data.frame(PRIOR.INIT)
  RES$PARAM$PRIOR.PARAM <- PRIOR.PARAM
  RES$PARAM$CONVERGENCE <- code.convergence
}
#RES$PARAM$Q <- Q

# Presentation of the statistics calculated for each couple
##RES$STATISTIC <- data.frame(n11,post.H0,LB,postE)
##rownames(RES$STATISTIC) <- paste(L[,1],L[,2]) 
##colnames(RES$STATISTIC) <- c("Effectif","postH0","Lower Bound","postE")

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
    colnames(RES$ALLSIGNALS) <- c("drug","event","count","expected count","postH0",
                                  "n11/E","drug margin","event margin","FDR","FNR","Se","Sp")
}
if (RANKSTAT==2 | RANKSTAT==3) {
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
  if (RANKSTAT==2) colnames(RES$ALLSIGNALS) <- c("drug","event","count","expected count",
                            "Q_0.05(lambda)","n11/E","drug margin","event margin","FDR","FNR","Se","Sp","postH0")
  if (RANKSTAT==3) colnames(RES$ALLSIGNALS) <- c("drug","event","count","expected count",
                            "post E(Lambda)","n11/E","drug margin","event margin","FDR","FNR","Se","Sp","postH0")

}

# List of Signals generated according to the method 
RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux,]

# FDR,FNR,Se,Sp
#RES$OpChar <- data.frame(FDR,FNR,Se,Sp)
#rownames(RES$OpChar) <- paste(RES$ALLSIGNALS[,1],RES$ALLSIGNALS[,2])

# Number of signals
RES$NB.SIGNALS <- Nb.signaux

RES  
}

