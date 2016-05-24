pca_scaled <-
function(D){

# --------------------
# data
# --------------------
Z0 <- as.matrix(D[,2:ncol(D)])
M_ID <- D[,1]

# ---------------------------
#error handling
# ---------------------------
Z0[is.na(Z0)] <- 0 # impute 0
index <- apply(Z0,1,sd)!=0 # checking

Z <- Z0[index,]
Z <- t(Z)

# ---------------------------
#autoscaling
# ---------------------------
m <- apply(Z,2,mean)# mean
s <- apply(Z,2,sd)# sd
n <- nrow(Z) # num. of samples
X <- (Z-(matrix(1,nrow=n)%*%m))/(matrix(1,nrow=n)%*%s)

# ------------------------------------------
#singular value decomposition
# ------------------------------------------
a <- svd(X/sqrt(n-1))
w <- a$v # eigenvector
t <- X%*%w # PC score

# ---------------------------------------------
#factor loading, p-value, q-value
# ---------------------------------------------
FL <- NaN;P <- NaN;Q <- NaN
fn <- NaN;pn <- NaN;qn <- NaN;cn <- NaN;tn <- NaN
for (i in 1:nrow(t)){
r <- NaN; p <- NaN
for (j in 1:ncol(X)){
R <- cor.test(t[,i],X[,j], method="pearson")
r[j] <- R$estimate
p[j] <- R$p.value
}
FL <- cbind(FL,r)

# p-value
P <- cbind(P, p)

# q-value
q <- p.adjust(p,"BH")
Q <- cbind(Q, q)

# -------------------
#label
# -------------------
tn[i] <- paste("PC score (PC", i, ")", sep="")
fn[i] <- paste("factor loading", "(PC", i ,")" , sep="")
pn[i] <- paste("p-value" , "(PC", i ,")" , sep="")
qn[i] <- paste("q-value" , "(PC", i, ")" , sep="")
cn[i] <- paste("Contribution ratio (PC", i, ")", sep="")
}

# ---------------------------------
#Contribution ratio
# ---------------------------------
b <- apply(t,2,var)
contribution.ratio <- 100*b/sum(b)
names(contribution.ratio) <- cn

# --------------------
#Result
# --------------------

# ---------------------
#PC score
# ---------------------
T <- data.frame(t)
colnames(T) <- tn

# ----------------------------
#factor loading
# ----------------------------
FL <- FL[,-1]
FL0 <- matrix(NA, nrow=nrow(Z0), ncol=ncol(Z0))
FL0[index,] <- FL
factor.loading <- data.frame(M_ID, FL0)
names(factor.loading) <- c("Metabolite ID", fn)

# ----------------------
#p-value
# ----------------------
P <- P[,-1]
P0 <- matrix(NA, nrow=nrow(Z0), ncol=ncol(Z0))
P0[index,] <- P
p.value <- data.frame(M_ID, P0)
names(p.value) <- c("Metabolite ID", pn)

# ---------------------
#q-value
# ---------------------
Q <- Q[,-1]
Q0 <- matrix(NA, nrow=nrow(Z0), ncol=ncol(Z0))
Q0[index,] <- Q
q.value <- data.frame(M_ID, Q0)
names(q.value) <- c("Metabolite ID", qn)

# ----------------------
# Return
# ----------------------
res <- list(T, factor.loading,p.value,q.value, contribution.ratio)
names(res) <- c("score", "factor.loading","p.value","q.value", "contribution.ratio")
return(res)
}
