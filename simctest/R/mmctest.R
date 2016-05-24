# source("simctest/R/simctest.r");

# mmctest -- Multiple Monte-Carlo Tests





# export: mmctest, hBH, hPC, hIndep
# some FDR controlling procedures

# FDR control by Benjamini-Hochberg
hBH <- function(p, threshold) {

    return(rank(p) <= max( c(which(sort(p)<=(1:length(p))*threshold/length(p)),-1) ));
}

# FDR control by Pounds&Cheng
hPC <- function(p, threshold) {

  m <- length(p);
  pi0_hat <- min(1, 2/m*sum(p));
  return( rank(p) <= max( c(which(sort(p)<=(1:length(p))*threshold/ (pi0_hat*m) ),-1) ));
}

# independent testing via Bonferroni
hBonferroni <- function(p, threshold) {

  m <- length(p);
  return(p<=threshold/m);
}

tau <- function(p,threshold) {
  m <- length(p);

  # step-up
  t <- threshold/rep(m,m);#Bonferroni
#   t <- (1:m)/m*threshold;#Simes
#   t <- threshold/(m+1-1:m);#Hochberg
#   t <- 1;#Rom
#   t <- (1:m)/m*threshold;#Benjamini-Hochberg
#   t <- (1:m)/m*(threshold/sum(1/(1:m)));#Benjamini-Yekutieli

  # step-down
#   t <- 1-(1-threshold)**(1/(m+1-(1:m)));#Sidak
#   t <- threshold/(m+1-(1:m));#Holm
#   t <- threshold/(m+1-(1:m));#Shaffer for a general scenario

  return(t);
}

hStepUp <- function(p, threshold) {
  return(rank(p) <= max(c(which(sort(p)<=tau(p,threshold)),-1)) );
}

hStepDown <- function(p, threshold) {
  return(rank(p) <= min(c(which(sort(p)>tau(p,threshold)), length(p))) );
}





# class mmctSamplerGeneric
# derive a class from mmctSamplerGeneric to implement the interface used to draw new samples
# or use class mmctSampler to directly pass a function and the number of hypotheses
# getSamples has to be implemented in such a way as to draw n[i] new samples for each
# hypothesis ind[i], where i=1...length(ind) and return the number of exceedances
# getNumber has to return the number of hypotheses
setClass("mmctSamplerGeneric")
setGeneric("getSamples", def=function(obj, ind, n){standardGeneric("getSamples")})
setGeneric("getNumber", def=function(obj){standardGeneric("getNumber")})





# class mmctSampler
# directly pass a function "f" and the number of hypotheses "num", class has a slot "data" for additional data
setClass("mmctSampler", contains="mmctSamplerGeneric", representation=representation(f="function", num="numeric", data="numeric"))

# implemented method
# for each hypotheses in vector ind[i], request n[i] new samples
setMethod("getSamples", signature(obj="mmctSampler"), function(obj, ind, n) {
    return(obj@f(ind, n, obj@data));
  }
)

# implemented method
setMethod("getNumber", signature(obj="mmctSampler"), function(obj) {
    return(obj@num);
  }
)

# pseudo constructor
mmctSampler <- function(f, num, data=NULL) {
  obj <- new("mmctSampler");
  obj@f <- f;
  obj@num <- num;
  obj@data <- data;
  return(obj);
}





# class mmctestres
# exportMethods: cont, show, pEstimate
setClass("mmctestres", contains="sampalgPrecomp", representation=representation(internal="environment", epsilon="numeric", threshold="numeric", r="numeric", R="numeric",
h="function", gensample="mmctSamplerGeneric", g="numeric", num="numeric", A="numeric", B="numeric", C="numeric", thompson="logical", rejprob="numeric"))

setGeneric("mainalg", def=function(obj, stopcrit){standardGeneric("mainalg")})
setMethod("mainalg", signature(obj="mmctestres"), function(obj, stopcrit) {

    g <- obj@g;
    num <- obj@num;
    batchN <- obj@internal$batchN;
    pu <- obj@internal$pu;
    pl <- obj@internal$pl;
    it <- obj@internal$it;

    A <- obj@A;
    B <- obj@B;
    C <- obj@C;

    m <- getNumber(obj@gensample);
    currentBatch <- max(batchN);
    timer <- proc.time()[[3]];
    copying <- 0;

    prior_alpha <- rep(1,m);
    prior_beta <- rep(1,m);

    tryCatch({
      while(length(B)>stopcrit$undecided) {

	# sample all \hat{p}_i in B
	currentBatch <- floor(1.25*currentBatch);
	if(currentBatch > 100000) { currentBatch <- 100000; }

	# Thompson Sampling: in each iteration, allocate samples to each undecided hypothesis
	# proportional to its probability of being randomly classified, i.e. proportional to
	# min(r_i,1-r_i), where r_i is the empirical probability of H_{0i} being rejected in R repetitions.
	# The parameter R can be set in the constructor (default value 1000).
	# Use Residual Sampling to allocate samples with weights min(r_i,1-r_i).
	if(obj@thompson==T) {
	  if((stopcrit$maxnum>0) && (stopcrit$maxit>0)) { constantBatch <- floor(stopcrit$maxnum/stopcrit$maxit); }
	  else { constantBatch <- currentBatch*m; }
	  if(it>0) {
	    prior_alpha <- 1+g;
	    prior_beta <- 1+num-g;
	    rej <- rowSums( replicate(obj@R, obj@h(rbeta(m,prior_alpha,prior_beta),obj@threshold)) );
	    sampleprob <- pmin(rej/obj@R,1-rej/obj@R);
	    
	    # correction
	    if(sum(sampleprob)==0) { sampleprob <- rep(1,m); }

	    # residual sampling
	    wiDelta <- constantBatch*sampleprob/sum(sampleprob);
	    batchN <- floor(wiDelta);
	    sampleprob <- wiDelta-batchN;
	    if(constantBatch-sum(batchN)>0) {
	      s <- sample(1:m,size=constantBatch-sum(batchN),replace=T,prob=sampleprob);
	      batchN <- batchN + hist(s, breaks=0:m+0.5, plot=F)$counts;
	    }
	  }
	  else { batchN[B] <- floor(constantBatch/m); }
	}
	else { batchN[B] <- currentBatch; }

	g[B] <- g[B] + getSamples(obj@gensample, B, batchN[B]);
	num[B] <- num[B] + batchN[B];

	# compute upper "exact" confidence level (Clopper-Pearson)
	a <- (num/(num+obj@r) - (num-batchN)/(num-batchN+obj@r)) * (obj@epsilon/m);
	pu_ <- pu;
	pl_ <- pl;
	qindex <- (g>0) & (g<num);
	pu[qindex] <- 1 - qbeta(a[qindex]/2, num[qindex]-g[qindex], g[qindex]+1);
	pu[g==0] <- 1-(a[g==0]/2)**(1/num[g==0]);
	pu[g==num] <- 1;

	# ...and lower "exact" confidence level (Clopper-Pearson)
	pl[qindex] <- 1 - qbeta(1-a[qindex]/2, num[qindex]+1-g[qindex], g[qindex]);
	pl[g==0] <- 0;
	pl[g==num] <- (a[g==num]/2)**(1/num[g==num]);

	# intersect confidence intervals and set confidence intervals of unconsidered hypotheses back to last value
	pu <- pmin(pu_, pu);
	pl <- pmax(pl_, pl);
	pu[setdiff(1:m,B)] <- pu_[setdiff(1:m,B)];
	pl[setdiff(1:m,B)] <- pl_[setdiff(1:m,B)];

	# compute classifications
	A <- which(obj@h(pu, obj@threshold));
	C <- which(obj@h(pl, obj@threshold));
	B <- setdiff(C, A);
	
	it <- it+1;

	copying <- 1;
	obj@g <- g;
	obj@num <- num;
	obj@internal$batchN <- batchN;
	obj@internal$pu <- pu;
	obj@internal$pl <- pl;
	obj@internal$it <- it;
	obj@A <- A;
	obj@B <- B;
	obj@C <- C;
	copying <- 0;

	if((stopcrit$maxnum>0) && (sum(num)+length(B)*floor(1.25*currentBatch)>=stopcrit$maxnum)) { break; }
	if((stopcrit$maxit>0) && (it>=stopcrit$maxit)) { break; }
	if((stopcrit$elapsedsec>0) && (proc.time()[[3]]-timer>=stopcrit$elapsedsec)) { break; }
      } # end while
    }, interrupt = function(interrupt){
      # finish copying if appropriate
      if(copying) {
	obj@g <- g;
	obj@num <- num;
	obj@internal$batchN <- batchN;
	obj@internal$pu <- pu;
	obj@internal$pl <- pl;
	obj@internal$it <- it;
	obj@A <- A;
	obj@B <- B;
	obj@C <- C;
	obj@rejprob <- NULL;
	return(obj);
      }
    }) # end tryCatch

    # update weights/rejection probabilities before termination
    prior_alpha <- 1+g;
    prior_beta <- 1+num-g;
    rej <- rowSums( replicate(obj@R, obj@h(rbeta(m,prior_alpha,prior_beta),obj@threshold)) );
    obj@rejprob <- rej/obj@R;
    return(obj);
  }
)

# cont offers four stopping criteria given as list "steps"
# steps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)
# corresonding to a maximal number of iterations, maximal number of samples drawn, number of yet undecided hypotheses, elapsed time
setMethod("cont", signature(data="mmctestres"), function(data, steps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)) {

    if(length(setdiff(ls(steps),c("maxit","maxnum","undecided","elapsedsec")))>0) { stop("parameter steps contains invalid data"); }

    if(length(steps)==0) { stopcrit <- list(maxit=0, maxnum=0, undecided=0, elapsedsec=0); }
    else {
      stopcrit <- list();
      if("maxit" %in% ls(steps)) { stopcrit<-c(stopcrit,maxit=steps$maxit); } else { stopcrit<-c(stopcrit,maxit=0); }
      if("maxnum" %in% ls(steps)) { stopcrit<-c(stopcrit,maxnum=steps$maxnum); } else { stopcrit<-c(stopcrit,maxnum=0); }
      if("undecided" %in% ls(steps)) { stopcrit<-c(stopcrit,undecided=steps$undecided); } else { stopcrit<-c(stopcrit,undecided=0); }
      if("elapsedsec" %in% ls(steps)) { stopcrit<-c(stopcrit,elapsedsec=steps$elapsedsec); } else { stopcrit<-c(stopcrit,elapsedsec=0); }
    }

    return(mainalg(data, stopcrit));
  }
)

setMethod("show", signature(object="mmctestres"), function(object) {

    m <- getNumber(object@gensample);
    cat(paste("Number of rejected hypotheses: ",length(object@A),"\n",sep=""));
    cat(paste("Number of non-rejected hypotheses: ",length(setdiff(1:m, object@C)),"\n",sep=""));
    cat(paste("Number of undecided hypotheses: ",length(object@B),"\n",sep=""));
    cat(paste("Total number of samples: ",sum(object@num),"\n",sep=""));
  }
)

setGeneric("pEstimate", def=function(obj){standardGeneric("pEstimate")})
setMethod("pEstimate", signature(obj="mmctestres"), function(obj) {

    return((obj@g+1)/(obj@num+1));
  }
)

setGeneric("rejProb", def=function(obj){standardGeneric("rejProb")})
setMethod("rejProb", signature(obj="mmctestres"), function(obj) {

    return(obj@rejprob);
  }
)

setGeneric("confidenceLimits", def=function(obj){standardGeneric("confidenceLimits")})
setMethod("confidenceLimits", signature(obj="mmctestres"), function(obj) {

    l <- list();
    l$lowerLimits <- obj@internal$pl;
    l$upperLimits <- obj@internal$pu;
    return(l);
  }
)

setGeneric("testResult", def=function(obj){standardGeneric("testResult")})
setMethod("testResult", signature(obj="mmctestres"), function(obj) {

    m <- getNumber(obj@gensample);
    l <- list();
    l$rejected <- obj@A;
    l$nonrejected <- setdiff(1:m, obj@C);
    l$undecided <- obj@B;
    return(l);
  }
)

setGeneric("summary.mmctestres", def=function(object,...){standardGeneric("summary.mmctestres")})
setMethod("summary.mmctestres", signature(object="mmctestres"), function(object) {

    cat(paste("Number of hypotheses: ",getNumber(object@gensample),sep=""));
    cat(strwrap(paste("Indices of rejected hypotheses:", paste(object@A,collapse=" ")),prefix="\n"));
    cat(strwrap(paste("Indices of undecided hypotheses:", paste(object@B,collapse=" ")),prefix="\n"));
    cat("\nAll hypotheses not listed are classified as not rejected.\n");
  }
)





# class mmctest
# exportMethods: run
setClass("mmctest", contains="mmctestres", representation=representation(internal="environment"))

setMethod("run", signature(alg="mmctest", gensample="mmctSamplerGeneric"), function(alg, gensample, maxsteps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)) {

    obj <- new("mmctestres");

    obj@epsilon <- alg@internal$epsilon;
    obj@threshold <- alg@internal$threshold;
    obj@thompson <- alg@internal$thompson;
    obj@r <- alg@internal$r;
    obj@R <- alg@internal$R;
    obj@h <- alg@internal$h;
    obj@gensample <- gensample;

    m <- getNumber(obj@gensample);
    obj@g <- rep(0, m);
    obj@num <- rep(0, m);
    obj@internal$batchN <- rep(10, m);
    obj@internal$pu <- rep(1, m);
    obj@internal$pl <- rep(0, m);
    obj@internal$it <- 0;

    obj@A <- 0;
    obj@B <- 1:m;
    obj@C <- 0;

    return(cont(obj, steps=maxsteps));
  }
)

# pseudo constructor
mmctest <- function(epsilon=0.01, threshold=0.1, r=10000, h, thompson=F, R=1000) {

  obj <- new("mmctest");
  obj@internal$epsilon=epsilon;
  obj@internal$threshold=threshold;
  obj@internal$thompson=thompson;
  obj@internal$r=r;
  obj@internal$R=R;
  obj@internal$h=h;
  return(obj);
}
