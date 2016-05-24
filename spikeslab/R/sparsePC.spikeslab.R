####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

sparsePC <- function(...)
{
  sparsePC.spikeslab(...)
}

sparsePC.spikeslab <-
  function(x=NULL,             #gene expressions
           y=NULL,             #class labels
           n.rep=10,           #no. replicates
           n.iter1=150,        #no. burn-in values
           n.iter2=100,        #no. Gibbs sampled values (following burn-in)
           n.prcmp=5,          #no. principal components
           max.genes=100,      #max genes allowed in signature
           ntree=1000,         #no. trees (forest classifier)
           nodesize=1,         #nodesize (forest classifier)
           verbose=TRUE,       #diagnostic output
            ...    )
{

### -------------------useful functions---------------------


permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)     
   nfolds= max(nfolds, 2)
   # makes no sense to have more folds than the max class size
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y) 
         # nice we to get the ids in a list, split by class
   ###Make a big matrix, with enough rows to get in all the folds per class
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
   ### Now do a clever sort to mix up the NAs
   smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
         # the "clever" unlist doesn't work when there are no NAs
         #       apply(smallmat, 2, function(x)
         #        x[!is.na(x)])
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
}

class.error <- function(y, ytest, pred) {
  cl <- sort(unique(y))
  err <- rep(NA, length(cl))
  for (k in 1:length(cl)) {
    cl.pt  <- (ytest == cl[k])
    if (sum(cl.pt) > 0) {
        err[k] <- mean(ytest[cl.pt] != pred[cl.pt])
    }
  }
  err
}

get.pc <- function(X, Y, n.prcmp) {
  cat("Getting prcmp...\n")
  #preliminary stuff
  nclass <- length(unique(Y))
  n <- nrow(X)
  p <- ncol(X)
  o.r <- order(Y)
  Y <- Y[o.r]
  Y.freq <- tapply(Y, Y, length)
  X <- X[o.r, ]
  #standardize
  X.mean <-  as.double(c(apply(X, 2, mean.center, center = TRUE)))
  X.sd <-  as.double(c(apply(X, 2, sd.center, center = TRUE)))
  ss.X <- scale(X, center = X.mean, scale = X.sd) 
  #get PC from SVD of X
  pc.out <- svd(ss.X)
  U <- pc.out$u
  D <- pc.out$d
  UD <- t(t(U)*D)
  #find pc with most variability across groups using RF-R
  #UD is better choice for "Y" than U
  imp <- apply(UD, 2, function(p){
     rf.out <- randomForest(cbind(Y), p, importance = TRUE)
     rf.out$importance[1,1]})
  o.r <- order(imp, decreasing = TRUE)
  #sort pc according to RF importance
  #enhance pc by drawing from within class N(m,0.1*s) 
  tot.var <- cumsum(D[o.r]/sum(D))
  cat("total variation for top prcmp(s):", round(100 * tot.var[n.prcmp]), "%", "\n")
  Y.prcmp <- as.matrix(UD[, o.r[1:n.prcmp]])
  Y.prcmp <- as.matrix(apply(Y.prcmp, 2, function(p) {
    p.mean <- rep(tapply(p, Y, mean), Y.freq)
    p.sd <- rep(tapply(p, Y, sd), Y.freq)
    rnorm(n, p.mean, 0.1 * p.sd) }
  ))
  colnames(ss.X) <- paste("X.", 1:p, sep="")
  return(list(Y.prcmp=Y.prcmp, ss.X=ss.X, tot.var=tot.var))
}


### -------------------get data ------------------------------
### convert Y to a factor
### dimensioning
### coherence checks

gene.signature <- NULL
X <- as.matrix(x)
Y <- y
if (!is.factor(y)) Y <- factor(y)
n.genes <- ncol(X)
n.data <- nrow(X)
if (any(is.na(Y))) stop("Missing values not allowed in y")
if (n.data != length(Y)) stop("number of rows of x should match length of y")


### -------------------Monte Carlo Loop----------------------------------
### repeat balanced CV (2/3, 1/3)
### sparse pc/ spikeslab analysis
### create centroid classifier
### compute PE over hold-out data (1/3)

n.class <- length(unique(Y))
dim.results <-  pred.results <- rep(0, n.rep)
pred.class.results <- matrix(NA, n.rep, n.class)

### standardize the expression data
X <- scale(X, center = TRUE, scale = TRUE)


for (k in 1:n.rep) {

  if (verbose) cat("\n ---> Monte Carlo Replication:", k, "\n")
  
  #cv sample
  if (n.rep > 1) {
    cv.sample <- balanced.folds(Y, 3)
    train.pt <- c(cv.sample[[1]], cv.sample[[2]])
    test.pt <- cv.sample[[3]]
  }
  else {
    train.pt <- test.pt <- 1:n.data
  }

  #calculate pc
  #call spikeslab
  #get significant genes
  n.prcmp <- min(n.prcmp, length(train.pt))
  pc.out <- get.pc(X[train.pt, ], Y[train.pt], n.prcmp = n.prcmp)  
  sig.genes <- NULL
  signal <- rep(0, n.genes)
  for (p in 1:n.prcmp) {
    if (verbose) cat("fitting principal component:", p, "(", round(100*pc.out$tot.var[p]), "%)", "\n")
    ss.out <- spikeslab(x = pc.out$ss.X, y = pc.out$Y.prcmp[, p],
                n.iter1 = n.iter1,
                n.iter2 = n.iter2,
                max.var = max.genes, 
                bigp.smalln.factor = max(1, round(max.genes/nrow(pc.out$ss.X))),
                bigp.smalln = (nrow(pc.out$ss.X) < ncol(pc.out$ss.X)))
    sig.genes.p <- as.double(which(abs(ss.out$gnet) > .Machine$double.eps))
    if (length(sig.genes.p) > 0) {
      sig.genes <- c(sig.genes, sig.genes.p)
      signal[sig.genes.p] <- signal[sig.genes.p] + abs(ss.out$gnet)[sig.genes.p]
    }
  }
  if (length(sig.genes) == 0) {
    sig.genes <- as.double(which(signal == max(signal, na.rm = TRUE)))[1]
  }
  sig.genes <- sort(unique(sig.genes))
  
  #forest classifier using spikeslab genes
  P <- min(max.genes, length(sig.genes))
  o.r <- order(signal[sig.genes], decreasing = TRUE)
  sig.genes.k <- sig.genes[o.r][1:P]
  rf.data.x <- as.matrix(X[, sig.genes.k])
  colnames(rf.data.x) <- paste("x.", 1:length(sig.genes.k))
  Y.train <- Y[train.pt]
  Y.test <- Y[test.pt]
  rf.out <- randomForest(x=as.matrix(rf.data.x[train.pt, ]), y=Y.train,
                   importance = TRUE, ntree=ntree, nodesize=nodesize)
  gene.signature <- c(gene.signature, sig.genes.k)
  dim.results[k] <- length(sig.genes.k)
  if (n.rep > 1) {
    rf.pred <- predict(rf.out, newdata = as.matrix(rf.data.x[test.pt, ]))
    pred.results[k] <- mean(as.character(Y.test) != rf.pred)
    pred.class.results[k, ] <- class.error(as.character(Y), as.character(Y.test), rf.pred)
  }
  
  #verbose
  if (verbose & (n.rep > 1)) {
    cat("\n", "PE:", round(pred.results[k], 3), "dimension:", dim.results[k], "\n")  
  }

}

### ------------------- spikeslab signature ------------------------------
### keep top mean(dim) occuring genes

gene.signature.all <- gene.signature
gene.signature.freq <- tapply(gene.signature, gene.signature, length)
gene.signature <- as.double(names(gene.signature.freq)[rev(order(gene.signature.freq))][1:mean(dim.results)])

### ------------------- final forest call ------------------------------

if (verbose) cat("growing the forest classifier...\n")
rf.data.x <- as.matrix(X[, gene.signature])
colnames(rf.data.x) <- paste("x.", gene.signature)
rf.out <- randomForest(x=rf.data.x, y=Y, importance = TRUE, ntree=ntree, nodesize=nodesize)


### ------------------- Output ------------------------------
###
###

cat("\n\n")
cat("-----------------------------------------------------------\n")
cat("no. prcmps           :", n.prcmp, "\n")
cat("no. genes            :", n.genes, "\n")
cat("max genes            :", max.genes, "\n")
cat("no. samples          :", nrow(X), "\n")
cat("no. classes          :", n.class, "\n")
cat("class freq           :", as.double(tapply(Y, Y, length)), "\n")
cat("class names          :", levels(Y), "\n")
cat("replicates           :", n.rep, "\n")
cat("model size           :", round(mean(dim.results), 4), "+/-", round(sd(dim.results), 4), "\n")
if (n.rep > 1) {
  cat("misclass             :", round(mean(100*pred.results), 4), "+/-", round(sd(100*pred.results), 4), "\n")
  for (j in 1:n.class) {
     cat(paste("  class #         ", j,  ":"), round(mean(100*pred.class.results[, j], na.rm=TRUE), 4), "\n")
  }
}
cat("\n")
cat("Gene Signature:\n")
print(gene.signature)
cat("-----------------------------------------------------------\n")


invisible(list(gene.signature=gene.signature,
               gene.signature.all=gene.signature.all,
               rf.object=rf.out))

}
