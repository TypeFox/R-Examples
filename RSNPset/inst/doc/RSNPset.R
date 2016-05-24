## ----setup2, include=FALSE------------------------------------------
options(width=70)  # make the printing fit on the page
set.seed(1121)     # make the results repeatable
stdt<-date()

## ----pts------------------------------------------------------------
set.seed(123)

n <- 100
status <- rbinom(n, 1, 0.5) 
table(status)

LDL <- rnorm(n, mean=115, sd=35)
quantile(LDL)

time <- rexp(n, 1/10)
event <- rbinom(n, 1, 0.9)
quantile(time)
table(event)

## ----covars---------------------------------------------------------
X <- matrix(c(rnorm(n), rbinom(n, 4, .25)), nrow=n)
summary(X[,1])
table(X[,2])

## ----snps-----------------------------------------------------------
m <- 500

G <- matrix(as.double(rbinom(n*m, 2, runif(n*m,.1,.9))), n, m)
dim(G)

## ----ids------------------------------------------------------------
rsIDs <- paste0("rs100",1:m)
colnames(G) <- rsIDs

G[1:5,1:5]

## -------------------------------------------------------------------
K <- 10
genes <- paste0("XYZ",1:K)
geneSets <- lapply(sample(3:50, size=K, replace=TRUE), sample, x=rsIDs)
names(geneSets) <- genes

unlist(lapply(geneSets, length))

## ----loadpkg--------------------------------------------------------
library(RSNPset)
set.seed(456)

## ----ccstat---------------------------------------------------------
ccres <- rsnpset(Y=status, G=G, snp.sets=geneSets, score="binomial", 
                 B=10, r.method="permutation", ret.rank=TRUE, v.permute=TRUE)

## ----ccres1---------------------------------------------------------
ccres[["Observed"]]

## ----ccsum----------------------------------------------------------
summary(ccres)

## ----ccpval---------------------------------------------------------
rsnpset.pvalue(ccres, pval.transform=TRUE)

## ----ldlstat--------------------------------------------------------
ldlres <- rsnpset(Y=LDL, G=G, X=X, snp.sets=geneSets, score="gaussian", B=10)

## ----ldlres2--------------------------------------------------------
ldlres[["Observed"]]
ldlres[["Replication.1"]]

## ----ldlpval--------------------------------------------------------
rsnpset.pvalue(ldlres)

## ----ttestat--------------------------------------------------------
tteres <- rsnpset(Y=time, delta=event, G=G, snp.sets=geneSets, score="cox", 
                  B=10, r.method="permutation", pinv.check=TRUE)

## ----tteres1--------------------------------------------------------
pinv.diag <- summary(tteres)
pinv.diag[["Observed"]]

unlist(lapply(pinv.diag, max))

## ----ttepval--------------------------------------------------------
ttepvals <- rsnpset.pvalue(tteres)
ttepvals

## ----ttepvalsum1----------------------------------------------------
summary(ttepvals, verbose=TRUE)

## ----ttepvalsum2----------------------------------------------------
ttesum <- summary(ttepvals, sort="pB", nrows=5, dropcols=c("m","Q","QB"))
ttesum

