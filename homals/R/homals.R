`homals` <-
function(data, ndim = 2, rank = ndim, level = "nominal", sets = 0, active = TRUE, 
eps = 1e-6, itermax = 1000, verbose = 0)
{
#data ... data frame
#sets ...  list of vectors of set indices 
#level ... which measurement level (either single string or vector
#ndim ... number of dimensions
#active ... which variables are active (single TRUE means all)
#rank ... which category quantification ranks (default all ndim)
#eps ... iteration precision eigenvalues (default 1e-6)

#----------------------------- constants --------------------------------

dframe <- data
name <- deparse(substitute(dframe))		# frame name
nobj <- nrow(dframe)					# number of objects
nvar <- ncol(dframe)					# number of variables
vname <- names(dframe)					# variable names
rname <- rownames(dframe)				# object names

#-----------------------------convert to factors, check data -------------------

for (j in 1:nvar) {
	dframe[, j] <- as.factor(dframe[, j])
	levfreq <- table(dframe[,j])
	if (any(levfreq == 0)) {
    newlev <- levels(dframe[, j])[-which(levfreq == 0)]
  } else {
    newlev <- levels(dframe[,j])
  }
  dframe[,j] <- factor(dframe[,j], levels = sort(newlev))
}

varcheck <- apply(dframe, 2, function(tl) length(table(tl)))	
if (any(varcheck == 1)) stop("Variable with only 1 value detected! Can't proceed with estimation!")
#-----------------------------parameter consistency-----------------------------

active <- checkPars(active,nvar)
rank <- checkPars(rank,nvar)
level <- checkPars(level,nvar)

if (length(sets) == 1) sets <- lapply(1:nvar,"c")
if (!all(sort(unlist(sets)) == (1:nvar))) {
	print(cat("sets union",sort(unlist(sets)),"\n"))
	stop("inappropriate set structure !")
	}
nset <- length(sets)

mis<-rep(0,nobj)
for (l in 1:nset) {
	lset<-sets[[l]]
	if (all(!active[lset])) next()
	jset<-lset[which(active[lset])]
	for (i in 1:nobj) {
		if (any(is.na(dframe[i,jset]))) 
			dframe[i,jset] <- NA
		else mis[i] <- mis[i] + 1
		}
	}
	
for (j in 1:nvar) {
	k<-length(levels(dframe[,j]))
	if (rank[j] > min(ndim,k-1)) rank[j]<-min(ndim,k-1)
	}
	
#----------------initialize scores and counters-----------------------------

x <- cbind(orthogonalPolynomials(mis,1:nobj,ndim))
x <- normX(centerX(x,mis),mis)$q
#x <- normX(centerX(x,mis),mis)$q*sqrt(nobj*nvar)           #norm to X'MX=nmI --> X are z-scores

y <- lapply(1:nvar, function(j) computeY(dframe[,j],x))
#y <- updateY(dframe,x,y,active,rank,level,sets)
sold <- totalLoss(dframe,x,y,active,rank,level,sets)
iter <- pops <- 0

#-----------------------------main computation--------------------------------

repeat {
	iter <- iter + 1
	y <- updateY(dframe,x,y,active,rank,level,sets,verbose=verbose)
	smid <- totalLoss(dframe, x, y, active, rank, level, sets)/(nobj * nvar * ndim)
	ssum <- totalSum(dframe, x, y, active, rank, level, sets)
        qv <- normX(centerX((1/mis)*ssum,mis),mis)
	z <- qv$q
        #z <- qv$q*sqrt(nobj*nvar)                                   #norm to var = 1 
  
	snew <- totalLoss(dframe, z, y, active, rank, level, sets)/(nobj * nvar * ndim)
	if (verbose > 0) cat("Iteration:",formatC(iter,digits=3,width=3),"Loss Value: ", formatC(c(smid),digits=6,width=6,format="f"),"\n")
	
        r <- abs(qv$r)/2                                             #eigenvalues
	ops <- sum(r)                                                #convergence criteria
        aps <- sum(La.svd(crossprod(x,mis*z),0,0)$d)/ndim

	if (iter == itermax) {
		stop("maximum number of iterations reached")
		}
	
	if (smid > sold) {
		warning(cat("Loss function increases in iteration ",iter,"\n"))
	}

	if ((ops - pops) < eps) break 
		else {x <- z; pops <- ops; sold <- smid}	
}

#-----------------------------store final version--------------------------------

#---------------------------- Cone restricted SVD ------------------------------
ylist<-alist<-clist<-ulist<-NULL
for (j in 1:nvar) {                        #final score computation based on SVD
  gg <- dframe[,j]
  c <- computeY(gg,z)                      #category centroids based on object scores z
  d <- as.vector(table(gg))
  lst <- restrictY(d,c,rank[j],level[j])   #based on category centroids
  y <- lst$y                               #category scores Y = ZA'
  a <- lst$a                               #category loadings (weights A)
  u <- lst$z                               #(single) category scores (low rank quantifications)
  ylist <- c(ylist,list(y))                #list of final category scores Y
  alist <- c(alist,list(a))                #list of category loadings A
  clist <- c(clist,list(c))                #list of category centroids C
  ulist <- c(ulist,list(u))                #list of low rank quantifications U (aka Z)
}

#--------------------------preparing/labeling output----------------------------

dimlab <- paste("D", 1:ndim, sep = "")
for (i in 1:nvar) {
  if (ndim == 1) {
    ylist[[i]] <- cbind(ylist[[i]])
    ulist[[i]] <- cbind(ulist[[i]])
    clist[[i]] <- cbind(clist[[i]])
    #alist[[i]] <- cbind(alist[[i]])
  }
  
  options(warn = -1)
  rnames <- sort(as.numeric((rownames(clist[[i]]))))             #convert row names into integers
  options(warn = 0)
  if ((any(is.na(rnames))) || (length(rnames) == 0)) rnames <- rownames(clist[[i]])
  if (!is.matrix(ulist[[i]])) ulist[[i]] <- as.matrix(ulist[[i]])
  
  rownames(ylist[[i]]) <- rownames(ulist[[i]]) <- rownames(clist[[i]]) <- rnames
  rownames(alist[[i]]) <- paste(1:dim(alist[[i]])[1])
  colnames(clist[[i]]) <- colnames(ylist[[i]]) <- colnames(alist[[i]]) <- dimlab
  colnames(ulist[[i]]) <- paste(1:dim(as.matrix(ulist[[i]]))[2])
}
names(ylist) <- names(ulist) <- names(clist) <- names(alist) <- colnames(dframe)
rownames(z) <- rownames(dframe)
colnames(z) <- dimlab
#alist.t <- lapply(alist,t)

#------ score and dummy matrix -------
dummymat <- as.matrix(expandFrame(dframe, zero = FALSE, clean = FALSE))         #indicator matrix
dummymat01 <- dummymat                           #final indicator matrix
dummymat[dummymat == 2] <- NA                    #missing observations
dummymat[dummymat == 0] <- Inf                   #irrelevant entries

scoremat <- array(NA, dim = c(dim(dframe), ndim), dimnames = list(rownames(dframe), colnames(dframe), paste("dim", 1:ndim, sep = "")))  #initialize array for score matrix (n*p*ndim)

for (i in 1:ndim) {
  catscores.d1 <-  do.call(rbind, ylist)[,i]       #category scores for dimension i
  dummy.scores <- t(t(dummymat) * catscores.d1)    #full score matrix (Inf, -Inf)

  freqlist <- apply(dframe, 2, function(dtab) as.list(table(dtab)))
  cat.ind <- sequence(sapply(freqlist, length))     #category indices (sequence)
  scoremat[,,i] <- t(apply(dummy.scores, 1, function(ds) {             #data matrix with category scores
                          ind.infel <- which(ds == Inf)                         #identify Inf entries
                          ind.minfel <- which(ds == -Inf)                       #identify -Inf entries
                          ind.nan <- which(is.nan(ds))
                          ind.nael <- which((is.na(ds) + (cat.ind != 1)) == 2)  #identify NA entries
                          ds[-c(ind.infel, ind.minfel, ind.nael, ind.nan)]               #return scored entries
                      } ))
}

#----------------------- discrimination measures ---------------------------
disc.mat <- apply(scoremat, 3, function(xx) {                                  #matrix with discrimination measures
                  apply(xx, 2, function(cols) {
                    (sum(cols^2, na.rm = TRUE))/nobj
                    })})

#---------------------- end discrimination measures ------------------------

result <- list(datname = name, catscores = ylist, scoremat = scoremat, objscores = z, 
               cat.centroids = clist, ind.mat = dummymat01, loadings = alist, 
               low.rank = ulist, discrim = disc.mat, ndim = ndim, niter = iter, level = level, 
               eigenvalues = r, loss = smid, rank.vec = rank, active = active, dframe = dframe, call = match.call())
class(result) <- "homals"
result
}

