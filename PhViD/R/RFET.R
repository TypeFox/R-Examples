`RFET` <-
function(DATABASE, OR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES=0.05, MID.PVAL = FALSE) {

# DATABASE :  object. It is the object returned by the function transform_data. It contains :
#             DATABASE$PARAM : the parameters used when calling the function transform_data
#             DATABASE$data :  matrix. The first column of DATA must contain the number of notifications n11, the second
#                              column the row marges n10 and the third column the column marges n01
#             DATABASE$N :     Nb de notifications total
#             DATABASE$L :     LIBELLES

# OR0 :  positive double. The value of the risk you want to consider. By default, OR0=1

# MIN.n11 : You can choose a minimum number of notifications contraint for a couple may be considered as a possible
#             signal. By default, MIN.n11 = 1. (ex: if MIN.n11=5, only couples with at least 5 notifications will
#             be able to be declared as an alert)

# DECISION :  You can choose which rule to use to determine the list of signals
#                 1 = FDR (ex : 0.05)
#                 2 = Number of Signals (ex : 1000)
#                 3 = pVALUE 

# DECISION.THRES : The value of the FDR / Number of signals, Ranking statistic, considering the DECISION rule you choose

# MID.PVAL :  Default = False

#OR0 <-c(1,2,5) ## seuil pour la definition des associations

require("LBE")


# Initialization              
DATA <- DATABASE$data
N <- DATABASE$N
L <- DATABASE$L

n11 <- DATA[,1]
n1. <- DATA[,2] # les marges lignes (effets indésirables)
n.1 <- DATA[,3] # les marges colonnes (médicaments)
n10 <- n1. - n11
n01 <- n.1 - n11
n00 <- N - (n11+n10+n01)
E <- n1. * n.1 / N # les effectifs attendus

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

logRFET <- log(n11 * n00 /(n10 * n01))


  pval.fish.uni <- vector(length=Nb.Cell)
  for (p in 1 : Nb.Cell) {
	 pval.fish.uni[p] <-  fisher.test(matrix(c(n11[p],n10[p],n01[p],n00[p]),ncol=2,byrow=TRUE),or=OR0,alternative="g")$p.value
  }

  if (MID.PVAL == TRUE) { # option MID.PVAL
    require(MCMCpack)
    for (p in 1 : Nb.Cell) {
	   pval.fish.uni[p] <- pval.fish.uni[p] - 0.5 * dnoncenhypergeom(x = n11[p], n1 = n11[p] + n01[p] , n2 = n10[p] + n00[p], m1 = n11[p] + n10[p], psi = OR0)
    }
  }
  
  pval.uni <- pval.fish.uni
  pval.uni[pval.uni>1] <-1
  pval.uni[pval.uni<0] <-0

  RankStat <- pval.uni
  LBE.res <- LBE(2 * apply(cbind(RankStat, 1-RankStat),1,min),plot.type="none")
  pi.c <- LBE.res$pi0

  fdr <- pi.c  * sort(RankStat[RankStat <= .5]) / (c(1:sum(RankStat <= .5)) / Nb.Cell)
  fdr <- c(fdr,
           pi.c /(2 *((sum(RankStat <= .5)+1) : Nb.Cell)/ Nb.Cell)
	         + 1 
	         - sum(RankStat <= .5)/((sum(RankStat <= .5)+1):Nb.Cell)
          )
  FDR <- apply(cbind(fdr,1),1,min)

# Calculation of the number of signals according to the decision rule (pp/FDR/Nb of Signals)
if (DECISION == 1) Nb.signaux <- sum(FDR <= DECISION.THRES)
if (DECISION == 2) Nb.signaux <- min(DECISION.THRES,Nb.Cell)
if (DECISION == 3) Nb.signaux <- sum(RankStat <= DECISION.THRES)


############################ SORTIE DE LA FONCTION #############################
RES <- vector(mode="list")

#RES$LIBEL <- L
#colnames(RES$LIBEL) <- c("DRUG","EVENT")

RES$INPUT.PARAM <- data.frame(OR0, MIN.n11, DECISION, DECISION.THRES, MID.PVAL,RANKSTAT=1)
#colnames(RES$INPUT.PARAM) <- c("OR0","notification threshold")

# Presentation of the statistics calculated for each couple
##RES$STATISTIC <- data.frame(n11,RankStat)
##rownames(RES$STATISTIC) <- paste(L[,1],L[,2]) # liste des libellés ingénue
##colnames(RES$STATISTIC) <- c("Effectif","pvalue")

# STATISTICAL VALUE TO BE CONSIDERED (used in function compare)
#RES$COMPARE <- vector(mode="list")
#RES$COMPARE$RANKSTAT <- 1
#RES$COMPARE$STAT <- RankStat

# SIGNALS RESULTS and presentation
RES$ALLSIGNALS <- data.frame( L[,1][order(RankStat)],
                              L[,2][order(RankStat)],
                              n11[order(RankStat)],
                              E[order(RankStat)],
                              RankStat[order(RankStat)],
                              exp(logRFET)[order(RankStat)],
                              n1.[order(RankStat)],
                              n.1[order(RankStat)],
                              FDR )
colnames(RES$ALLSIGNALS) <- c("drug code","event effect","count","expected count","p-value",
                              "ROR","drug margin","event margin","FDR")

# List of Signals generated according to the DECISION.THRES 
RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux,]

# FDR,FNR,Se,Sp
#RES$OpChar <- as.matrix(FDR)
#rownames(RES$OpChar) <- paste(RES$ALLSIGNALS[,1],RES$ALLSIGNALS[,2])
#colnames(RES$OpChar) <- "FDR"

# Number of signals
RES$NB.SIGNALS <- Nb.signaux

RES
}

