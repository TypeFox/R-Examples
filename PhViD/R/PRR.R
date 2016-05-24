`PRR` <-
function(DATABASE, RR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES = 0.05, RANKSTAT = 1) {

# DATABASE :  object. It is the object returned by the function transform_data. It contains :
#             DATABASE$PARAM : the parameters used when calling the function transform_data
#             DATABASE$data :  matrix. The first column of DATA must contain the number of notifications n11, the second
#                              column the row marges n10 and the third column the column marges n01
#             DATABASE$N :     Nb de notifications total
#             DATABASE$L :     LIBELLES

# RR0 :        positive double. The value of the risk you want to consider. By default, RR0=1

# MIN.n11 : You can choose a minimum number of notifications contraint for a couple may be considered as a possible
#             signal. By default, MIN.n11 = 1. (ex: if MIN.n11=5, only couples with at least 5 notifications will
#             be able to be declared as an alert)

# DECISION :  You can choose which rule to use to determine the list of signals
#                 1 = FDR (ex : 0.05) : This decision is only available when the P-values are chosen as ranking stat SEE RANKSTAT
#                 2 = Number of Signals (ex : 1000)
#                 3 = Ranking statistic See RankStat 

# DECISION.THRES : The value of the FDR / Number of signals, Ranking statistic, considering the DECISION rule you choose

# RankStat:   You can choose which statistic to use to order the signals :
# (expert)        1 = pVALUE
#                 2 = CI Lower Bound (95%) of ln(PRR)


require("LBE")
if (RANKSTAT==2 & DECISION == 1) stop("The FDR can't be used as decision rule with this ranking Statistic")
# Initialization              
DATA <- DATABASE$data
N <- DATABASE$N
L <- DATABASE$L

n11 <- DATA[,1]
n1. <- DATA[,2] # marginal drug counts
n.1 <- DATA[,3] # marginal AE counts
n10 <- n1. - n11
n01 <- n.1 - n11
n00 <- N - (n11+n10+n01)
E <- n1. * n.1 / N 

if(MIN.n11 > 1) {
  E <- E[n11 >= MIN.n11]
  n1. <- n1.[n11 >= MIN.n11]
  n.1 <- n.1[n11 >= MIN.n11]
  n10 <- n10[n11 >= MIN.n11]
  n01 <- n01[n11 >= MIN.n11]
  n00 <- n00[n11 >= MIN.n11]
  LL <- data.frame(drugs=L[,1],events=L[,2],n11)
  LL1 <- LL[,1][n11 >= MIN.n11]
  LL2 <- LL[,2][n11 >= MIN.n11]
  rm(list="L")
  L <- data.frame(LL1,LL2)
  n11 <- n11[n11 >= MIN.n11]
}

Nb.Cell <- length(n11)

#logPRR <- log( (n11 / (n11 + n01)) /  (n10 / (n10 + n00)) )
#var.logPRR <- 1/n11 - 1/(n11 + n01) + 1/n10 - 1/(n10 + n00)
logPRR <- log( (n11 / (n11 + n10)) /  (n01 / (n01 + n00)) )
var.logPRR <- 1/n11 - 1/(n11 + n10) + 1/n01 - 1/(n01 + n00)
pval.logPRR.uni<- 1-pnorm(logPRR,log(RR0),sqrt(var.logPRR))
petit_rankstat <- (logPRR-log(RR0))/sqrt(var.logPRR) 
pval.uni       <- pval.logPRR.uni

pval.uni[pval.uni>1] <-1
pval.uni[pval.uni<0] <-0

PVAL.UNI <- pval.uni
LBE.res <- LBE(2 * apply(cbind(pval.uni, 1-pval.uni),1,min),plot.type="none")
pi.c <- LBE.res$pi0
fdr <- pi.c  * sort(pval.uni[pval.uni <= .5]) / (c(1:sum(pval.uni <= .5)) / Nb.Cell)
fdr <- c(fdr,
         pi.c /(2 *((sum(pval.uni <= .5)+1) : Nb.Cell)/ Nb.Cell)
        + 1 
        - sum(pval.uni <= .5)/((sum(pval.uni <= .5)+1):Nb.Cell)
        )
FDR <- apply(cbind(fdr,1),1,min)

if (RANKSTAT==2){FDR <- rep(NaN,length(n11))}

# Calculation of the Lower Bound
LB <- qnorm(0.025,logPRR,sqrt(var.logPRR))


if (RANKSTAT==1) RankStat <- PVAL.UNI
if (RANKSTAT==2) RankStat <- LB

# Calculation of the number of signals according to the decision rule (pval/FDR/Nb of Signals)
if (DECISION == 1 & RANKSTAT==1) Nb.signaux <- sum(FDR <= DECISION.THRES)
if (DECISION == 2) Nb.signaux <- min(DECISION.THRES,Nb.Cell)
if (DECISION == 3) {
    if (RANKSTAT==1) Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT==2) Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
}


############################ SORTIE DE LA FONCTION #############################
RES <- vector(mode="list")

#RES$LIBEL <- L
#colnames(RES$LIBEL) <- c("DRUG","EVENT")

RES$INPUT.PARAM <- data.frame(RR0, MIN.n11, DECISION, DECISION.THRES, RANKSTAT)
#colnames(RES$INPUT.PARAM) <- c("RR0","notification threshold")

# Presentation of the statistics calculated for each couple
##RES$STATISTIC <- data.frame(n11,PVAL.UNI,LB)
##rownames(RES$STATISTIC) <- paste(L[,1],L[,2]) # liste des libellés ingénue
##colnames(RES$STATISTIC) <- c("Effectif","pvalue","Lower Bound")

# STATISTICAL VALUE TO BE CONSIDERED (used in function compare)
#RES$COMPARE <- vector(mode="list")
#RES$COMPARE$RANKSTAT <- RANKSTAT
#RES$COMPARE$STAT <- RankStat

# SIGNALS RESULTS and presentation
#if (RANKSTAT==1) {
    RES$ALLSIGNALS <- data.frame( L[,1][order(petit_rankstat, decreasing=TRUE)],
                                  L[,2][order(petit_rankstat, decreasing=TRUE)],
                                  n11[order(petit_rankstat, decreasing=TRUE)],
                                  E[order(petit_rankstat, decreasing=TRUE)],
                                  RankStat[order(petit_rankstat, decreasing=TRUE)],
                                  exp(logPRR)[order(petit_rankstat, decreasing=TRUE)],
                                  n1.[order(petit_rankstat, decreasing=TRUE)],
                                  n.1[order(petit_rankstat, decreasing=TRUE)],
                                  FDR )
    colnames(RES$ALLSIGNALS) <- c("drug code","event effect","count","expected count","p-value",
                                  "PRR","drug margin","event margin","FDR")
if (RANKSTAT==2){colnames(RES$ALLSIGNALS)[5] <- "LB95(log(PRR))"}  
#}
#if (RANKSTAT==2 & DECISION != 1) {
#    RES$ALLSIGNALS <- data.frame( L[,1][order(RankStat, decreasing=TRUE)],
#                                  L[,2][order(RankStat, decreasing=TRUE)],
#                                  n11[order(RankStat, decreasing=TRUE)],
#                                  E[order(RankStat,decreasing=TRUE)],
#                                  RankStat[order(RankStat,decreasing=TRUE)],
#                                  exp(logPRR)[order(RankStat,decreasing=TRUE)],
#                                  n1.[order(RankStat,decreasing=TRUE)],
#                                  n.1[order(RankStat,decreasing=TRUE)],FDR)
#    colnames(RES$ALLSIGNALS) <- c("drug code","event effect","count","expected count","CI Lower Bound",
#                                  "ROR","drug margin","event margin")
#}
# List of Signals generated according to the DECISION.THRES 
RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux,]

# FDR,FNR,Se,Sp
#RES$OpChar <- as.matrix(FDR)
#if (RANKSTAT==2){RES$OpChar <- matrix(NaN,nr=length(n11),nc=1)}
#rownames(RES$OpChar) <- paste(RES$ALLSIGNALS[,1],RES$ALLSIGNALS[,2])
#colnames(RES$OpChar) <- "FDR"

# Number of signals
RES$NB.SIGNALS <- Nb.signaux

RES
}

