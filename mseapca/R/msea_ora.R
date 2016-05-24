msea_ora <-
function (SIG, ALL, M){

ALL <- as.character(as.matrix(ALL))
SIG <- as.character(as.matrix(SIG))
num_all <- length(ALL)
num_sig <- length(SIG)

# --------------------------------------------------------------
#Generating label matrix for detected metabolites 
# --------------------------------------------------------------

Lall0 <- setlabel(ALL,M)

# delete metabolite set
Lall <- Lall0[,colSums(Lall0)!=0]

# error handling
if (ncol(Lall) < 2){
# error message
stop("more than two metabolite set are necessary")
# stop function
return()
}

# --------------------------------------------------------------
#Generating label matrix for significant metabolites 
# --------------------------------------------------------------

Lsig <- setlabel(SIG,M)
Lsig <- Lsig[,colSums(Lall0)!=0]

l <- colSums(Lall0)!=0

# -------------------------------
#Calculating  ORA
# -------------------------------

P<-NaN;
for (i in 1:sum(l)){

# ------------------------------------
#Generating 2~2 table
# -------------------------------------
a1 <- sum(Lsig[,i])# significant and including pathway
a2 <- sum(Lall[,i])-sum(Lsig[,i])# not significant and including pathway
a3 <- length(SIG)-a1# significant and not including pathway
a4 <- (length(ALL)-length(SIG))-a2# not significant and not including pathway

tab <- t(matrix(c(a1,a2,a3,a4),2));

# ----------------------------------
# Fisher's exact test
# ----------------------------------
resfish <- fisher.test(tab, alternative="greater")
P[i] <- resfish$p.value

}

# -----------------------
#q-value
# -----------------------
Q <- p.adjust(P, method="BH")

# --------------------------------------------------------
#significant metabolites for metabolite set
# --------------------------------------------------------
LES <- NaN
for (i in 1:ncol(Lsig)){
les <- SIG[Lsig[,i]==1]
LES[i] <- list(les)

}
names(LES) <- colnames(Lsig)

# ----------------------
#Result
# ----------------------
PQ <- cbind(P,Q)
rownames(PQ) <- colnames(Lsig)
colnames(PQ) <- c("p.value","q.value")

RES <- list(PQ, LES)
names(RES) <- c("Result of MSEA(ORA)","significant metabolites")

# -------------------
#Return
# -------------------

return(RES)

}
