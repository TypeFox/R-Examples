msea_sub <-
function (M, D, y, maxiter=1000){

# ---------------------
# Data
# ---------------------
M_ID <- D[,1]# annotation
X <- as.matrix(D[,2:ncol(D)])# data matrix

# ---------------------------
#error process
# ---------------------------
X[is.na(X)] <- 0 # impute 0
index <- sd(t(X))!=0 # checking

X <- X[index,]
M_ID <- M_ID[index]

S_name <- names(M) # metabolite set name

# -----------------------------------------------------
#Generating [set~metabolite ID] matrix L
# -----------------------------------------------------
L <- setlabel(M_ID,M)

# delete metabolite set
L <- L[,colSums(L)!=0]

# error handling
if (ncol(L) < 2){
# error message
stop("more than two metabolite set are necessary")

# stop function
return()
}

# -----------------------------------
# ES for observed data
# -----------------------------------

ESMAT <- NaN

# -------------------------
# correlation
# -------------------------

ro <- NaN
for (i in 1:nrow(X)){
        R <- cor.test(X[i,],y)
ro[i] <- R$estimate
}
        
# descending order
d <- sort.list(ro,decreasing=TRUE) 
r <- sort(ro, decreasing=TRUE)

# ----------------------------
# Enrichment score
# ----------------------------

Lo <- L[d,]# label matrix
Mo <- M_ID[d]   # metabolite IDs

# //  Calculation of Phit, Pmiss
p <- 1; NR <- (abs(r)^(p))%*%Lo

ESo <- NaN
Z <- NaN
for (S in 1:length(NR)){

NH <- colSums(Lo);
        N <- nrow(Lo);
        
Phit <- NaN; Pmiss <- NaN
for (i in 1:length(r)){
            Phit[i] <- (abs(r[1:i])^(p))%*%Lo[1:i,S]/NR[S]
            Pmiss[i] <- sum(Lo[1:i,S]==0)*(1/(N-NH[S]))
        }

        z <- Phit-Pmiss
Z <- cbind(Z,z)

maxdev <- max(abs(Phit-Pmiss))
maxdev_index <- which.max(abs(Phit-Pmiss))

        ESo[S]=sign(Phit[maxdev_index]-Pmiss[maxdev_index])*maxdev
    }

Z <- Z[,-1]

# ------------------------------------
# ES for null hypothesis
# ------------------------------------

for (iter in 1:maxiter){

        # --------------------------------------
        #  randomized of label matrix
        # --------------------------------------

Lr <- NaN
        for (i in 1:ncol(L)){
             Lr <- cbind(Lr, L[sample(d),i])
        }

Lr <- Lr[,-1] # randmized label matrix
        
# ------------------------------
        #Enrichment score
# ------------------------------

        NR <- (abs(r)^(p))%*%Lr

ESr <- NaN
for (S in 1:length(NR)){
NH <- colSums(Lr)
N <- nrow(Lr)

Phit <- NaN; Pmiss <- NaN
for (i in 1:length(r)){
                Phit[i] <- (abs(r[1:i])^(p))%*%Lr[1:i,S]/NR[S]
                Pmiss[i] <- sum(Lr[1:i,S]==0)*(1/(N-NH[S]))
}

maxdev <- max(abs(Phit-Pmiss))
maxdev_index <- which.max(abs(Phit-Pmiss))

ESr[S] <- sign(Phit[maxdev_index]-Pmiss[maxdev_index])*maxdev

}

ESMAT <- cbind(ESMAT, ESr)
}

ESMAT <- ESMAT[,-1]

# ----------------------
# p-value
# ----------------------

P <- NaN
for (i in 1:nrow(ESMAT)){

esmat <- ESMAT[i,]
                
        # Positive value of NES
        if (ESo[i] >= 0){
        p_esmat <- esmat[esmat >= 0]
            P[i] <- sum(p_esmat >= ESo[i])/length(p_esmat)
        }
        
        # Negative value of NES
        if (ESo[i] < 0){
            n_esmat <- esmat[esmat <= 0]
            P[i] <- sum(n_esmat <= ESo[i])/length(n_esmat)
        }
    }

# --------------------------
# Normalized ES
# --------------------------

NESMAT <- NaN
nesmat <- NaN
NES <- NaN
for (i in 1:nrow(ESMAT)){
        
        esmat <- ESMAT[i,]
        
        p_esmat <- esmat[esmat>=0]
        n_esmat <- esmat[esmat<=0]
        
        # // Normalized enrichment score for null distribution
        nesmat[esmat >= 0] <- p_esmat/mean(p_esmat);
        nesmat[esmat < 0]  <- n_esmat/abs(mean(n_esmat));

NESMAT <- cbind(NESMAT,nesmat)

        # // Positive value of NES
        if (ESo[i] >= 0){
            NES[i] <- ESo[i]/mean(p_esmat);
        }
        
        # // Negative value of NES
        if (ESo[i] < 0){
        NES[i] <- ESo[i]/abs(mean(n_esmat));
        }
}

NESMAT <- NESMAT[,-1]

# -------------------
# q-value
# -------------------

Q1 <- NaN;Q2 <- NaN
K <- length(NES)
for (l in 1:length(NES)){

# ----------------
#  Q1
# ----------------
# // Positive value of NES
        if (ESo[l] >= 0){Q1[l] <- sum(sum((NES[l] <= NESMAT)))/sum(sum((NESMAT >= 0)))}

        # // Negative value of NES
        if (ESo[l] < 0){Q1[l] <- sum(sum((NES[l] >= NESMAT)))/sum(sum((NESMAT <= 0)))}

# ----------------
#  Q2
# -----------------
# // Positive value of NES
if (NES[l] >= 0){Q2[l] <- sum(NES[l] <= NES)/sum(NES >= 0)}

        # // Negative value of NES
        if(NES[l] < 0){Q2[l] <- sum(NES[l] >= NES)/sum(NES <= 0)}
}

Q <- Q1/Q2
Q[Q > 1] <- 1

PQ <- cbind(NES, P,Q)
colnames(PQ) <- c("normalized enrichment score","p-value","q-value")
rownames(PQ) <- colnames(Lo)

# ----------------------------------
#Leading edge subset
# ----------------------------------

M <- NaN
for (i in 1:ncol(Z)){

    if (ESo[i] > 0){
        m <- max(Z[,i])
n <- which.max(Z[,i])
    }
    if (ESo[i] < 0){
        m <- min(Z[,i])
n <- which.min(Z[,i])
    }
    M[i] <- n # // position

}

# // metabolite set
LES <- NaN
for (i in 1:ncol(Lo)) {
    
# // Positive value of NES
    if (ESo[i] > 0){
m <- Mo[1:M[i]] # metabolites
les <- m[Lo[1:M[i],i]==1]
    }

    # // Negative value of NES
    if (ESo[i] < 0){
        m <- Mo[-(1:M[i]-1)] 
les <- m[Lo[M[i]:nrow(Lo),i]==1]
    }

LES[i] <- list(les)
}

names(LES) <- colnames(Lo)

# ---------------------
#return
# ---------------------

g <- list(PQ,LES)
names(g) <- c("Result of MSEA", "leading edge subset")

return(g)

}
