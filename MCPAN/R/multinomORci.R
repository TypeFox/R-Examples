

Y2mat <- function(Ymat){
  if(is.table(Ymat)){
    if(length(dim(Ymat))!=2){
      stop(paste("Ymat should be a two-dimensional contingency table but input has", length(dim(Ymat)), "dimensions"))}else{YMAT <- as.matrix(Ymat)}}else{
                   if(is.data.frame(Ymat)){
                     wcn <- unlist(lapply(as.list(Ymat), function(x){is.numeric(x)|is.integer(x)}))
                     if(length(wcn)>sum(wcn)){warning("Some columns in Ymat do not contain numbers and will be omitted.")}
                     if(sum(wcn)>=2){YMAT <- data.matrix(Ymat[,wcn])}else{stop("Ymat should contain at least 2 columns with numeric/integer entries")}
                   }else{
                     if(is.matrix(Ymat)){
                       if(length(dim(Ymat))!=2){stop(paste("Ymat should be a two-dimensional contingency table but input has",
                                                                 length(dim(Ymat)), "dimensions"))}else{ YMAT <- Ymat}
                     }else{stop("Ymat should be a 2-dimensional matrix but input can not be coerced to a 2 dimensional matrix")}
                     
                   }
                 }
 
  if(is.numeric(YMAT)|is.integer(YMAT)){return(as.matrix(YMAT))}else{stop("Ymat must have integer (or numeric) entries")} 
  
}



datcmcheck <- function(Ymat, cmcat=NULL, cmgroup=NULL, ...)
{
  
  Ymat <- Y2mat(Ymat)
  if(!(is.numeric(Ymat)|is.integer(Ymat))){stop("Ymat must have integer (or numeric entries)")}
  if(any(Ymat-round(Ymat) >=  10 * .Machine$double.eps)){warning("Ymat should have whole numbers as entries")}

    ncat <- ncol(Ymat)
    ngrp <- nrow(Ymat)
    
    catni <- colSums(Ymat)
    groupni <- rowSums(Ymat)    
    
  if(is.null(cmcat)){
    CMCAT <- contrMat(n=catni, type="Dunnett", ...) 
  }else{
  if(is.character(cmcat)){
    ctype <- match.arg(cmcat, choices=c("Dunnett", "Tukey", "Sequen"))
    CMCAT <- contrMat(n=catni, type=ctype, ...)
  }else{
  if(!is.matrix(cmcat)){stop("cmcat must be a matrix or a single character string")}else{
  if(!(is.numeric(cmcat)|is.integer(cmcat))){stop("cmcat must have integer (or numeric entries)")}else{
  CMCAT <- cmcat  
  }}}}

    if(is.null(cmgroup)){
      CMGROUP <- contrMat(n=groupni, type="Dunnett", ...) 
    }else{
      if(is.character(cmgroup)){
        ctypegroup <- match.arg(cmgroup, choices=c("Dunnett", "Tukey", "Sequen"))
        CMGROUP <- contrMat(n=groupni, type=ctypegroup, ...)
      }else{
        if(!is.matrix(cmgroup)){stop("cmgroup must be a matrix or a single character string")}else{
          if(!(is.numeric(cmgroup)|is.integer(cmgroup))){stop("cmgroup must have integer (or numeric entries)")}else{
            CMGROUP <- cmgroup
          }}}}
    

    if(ncol(CMCAT) != ncat){stop("Number of columns in cmcat must be the same as number of columns (categories) in Ymat")}
    if(ncol(CMGROUP) != ngrp){stop("Number of columns in cmgroup must be the same as number of rows (groups) in Ymat")}
    
    RScat <- rowSums(CMCAT)
    RScatpos <- apply(CMCAT, 1, function(x){sum(x[sign(x)==+1])})
    
    RSgrp <- rowSums(CMGROUP)
    RSgrppos <- apply(CMGROUP, 1, function(x){sum(x[sign(x)==+1])})
    
    if(any(RScat>=10 * .Machine$double.eps)){warning("Some rows of 'contrast' matrix cmcat do not sum to 0.")}
    if(any(RSgrp>=10 * .Machine$double.eps)){warning("Some rows of 'contrast' matrix cmgrp do not sum to 0.")}
 
    if(any((RScatpos-1)>=10 * .Machine$double.eps)){cat("Note: The positive coefficients in some rows of 'contrast' matrix cmcat do not sum to 0.")}
    if(any((RSgrppos-1)>=10 * .Machine$double.eps)){cat("Note: The positive coefficients in some rows of 'contrast' matrix cmgrp do not sum to 0.")}
  return(list(Ymat=Ymat, cmgroup=CMGROUP, cmcat=CMCAT))  
}
  
  


##################################
# helper functions

addmat <- function(Ymat, alpha=1, type="group"){
  ngroup <- nrow(Ymat); ncat <- ncol(Ymat) 
  switch(type,
         "group"={priormat <- matrix(rep(alpha/ncat, times=ngroup*ncat),ncol=ncat)},
         "cell"={priormat <- matrix(rep(alpha, times=ngroup*ncat),ncol=ncat)},
         "total"={priormat <- matrix(rep(alpha/(ncat*ngroup), times=ngroup*ncat),ncol=ncat)})
  return(priormat)
}


makediriprior <- function(Ymat, prior){
if(is.null(prior)){PRIOR <- addmat(Ymat=Ymat, alpha=1, type="cell")}else{
  
  if(!(is.numeric(prior)|is.integer(prior))){stop("prior must be numeric")}

  ncat <- ncol(Ymat)
  ngrp <- nrow(Ymat)
  
  if(!is.matrix(prior)){
     if(length(prior)==1){ PRIOR <- addmat(Ymat=Ymat, alpha=prior, type="cell")}else{
     if(length(prior)==ncat){ PRIOR <- matrix(rep(prior, times=nrow(Ymat)), ncol=ncat)}else{
     if(length(prior)==(ncat*ngrp)){PRIOR <- matrix(prior, ncol=ncat)}else{stop("Prior is not a single number nor a vector, nor a matrix")
     }}}
  }else{
    if(ncol(prior)==ncol(Ymat) & nrow(prior)==nrow(Ymat)){
    PRIOR <- prior
    }else{stop("prior does not match dimensions of Ymat, and can not be coerced to match it")}}}
  if(is.null(rownames(PRIOR))){rownames(PRIOR) <- rownames(Ymat)}
  if(is.null(colnames(PRIOR))){colnames(PRIOR) <- colnames(Ymat)}
  return(PRIOR)
}


# paramat: a matrix with B=many samples 
# of parameters (columns) from a posterior distrribution
# or a bootstrapping process with B bootstrap replications

# CM a contrast matrix; ncol(CM) must be equal ncol(paramat)
# CM is %*%-multiplied with each row of paramat


# METHOD: sample from the dirichlet-posterior and
# compute SCS acc. to Besag et al. (1995)

# Ymat: data matrix
# prior: prior matrix
# cmcat: def. log-odds difference
# cmgroup def. between-group comparisons
# BSIM: No. of samples from dirichlet posterior

# bycomp: how to order parameters: bycomp=TRUE, parameters primarily ordered by betweengroup-comparisons, sec. by odds
#         bycomp=TRUE, parameters primarily ordered by betweengroup-comparisons, sec. by odds
#         bycomp=FALSE, parameters primarily ordered by odds, sec. by between-group-comparisons

# returns SCS at logit-scale!

############################################################

# SCI-method

simpostdiri <- function(Ymat, prior=NULL, cmcat=NULL, cmgroup=NULL, BSIM=10000,
 alternative="two.sided", conf.level=0.95, bycomp=FALSE, bychr=" btw ", ...){
 

YCG <- datcmcheck(Ymat=Ymat, cmcat=cmcat, cmgroup=cmgroup, ...)

YMAT <- YCG$Ymat
CMCAT <- YCG$cmcat
CMGROUP <- YCG$cmgroup

PRIOR <- makediriprior(Ymat=YMAT, prior=prior)

if(!is.numeric(BSIM)||length(BSIM)!=1){stop("BSIM must be a single, large integer")}
alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
if(!is.numeric(conf.level)||length(conf.level)!=1||(conf.level>=1 |conf.level<=0)){stop("conf.level must be a single numeric, e.g. 0.95")}


# apply CM%*%vector to each row of a 2-dimensional array (matrix)
CM1apply <- function(paramat, CM){apply(paramat, MARGIN=1, function(x){CM%*%matrix(x, ncol=1)})}

# result is a matrix with nrow=nrow(paramat) and ncol=nrow(CM)

# apply CM%*%vector to the third dimension vector of row of a 3-dimensional array
CM12apply <- function(paraarr, CM){apply(paraarr, MARGIN=1:2, function(x){CM%*%matrix(x, ncol=1)})}

# no. of groups in data
ngr <- nrow(YMAT)

# no. of categories in data
ncat <- ncol(YMAT)

# no. of odds (w.rt. categories)
nodds <- nrow(CMCAT)

# no. of between-group-comp.
ncomp <- nrow(CMGROUP)

compnames <- rownames(CMGROUP)
oddsnames <- rownames(CMCAT)

# use data + PRIOR as the dirichlet posterior
postpar <- YMAT+PRIOR

# define different 3-dimensional arrays

# sample of the dirichlet posterior
postarr <- array(dim=c(BSIM, ncat, ngr))

# sample of log-odds (for each treatment group separately)
postlogodds <- array(dim=c(BSIM, nodds, ngr))

# sample of log-odds differences (for each treatment group separately)
postdifflogit <- array(dim=c(BSIM, nodds, ncomp))

# sample of differences (comp. between treatments) for each log-odds
postcomplogit <- array(dim=c(BSIM, nodds*ncomp, ngr))

# sample from the dirichlet posterior, separately 
# 3 dim. array: BSIM x K x I
# (i.e. independently) for each treatment group (i=1,...,I)
for(i in 1:ngr){postarr[,,i] <- rdirichlet(n=BSIM, alpha=postpar[i,])}
for(i in 1:ngr){postlogodds[,,i] <- t(CM1apply(paramat=log(postarr[,,i]), CM=CMCAT))}

# compute samples of the between-group differences of the logits
for(d in 1:nodds){
postdifflogit[,d,] <- t(CM1apply(paramat=postlogodds[,d,], CM=CMGROUP))
}

# coerce the MxD matrix in each of the B samples into a vector of length=M*D

LNAM <- oddstrt2logitnames(cmcat=CMCAT, cmgroup=CMGROUP, bycomp=bycomp, bychr=bychr)

if(bycomp){
postlogitmat <- t(apply(X=postdifflogit, MARGIN=1, function(x){as.numeric(x)}))
colnames(postlogitmat) <- LNAM$logitnames
}else{
postlogitmat <- t(apply(X=postdifflogit, MARGIN=1, function(x){as.numeric(t(x))}))
colnames(postlogitmat) <- LNAM$logitnames}

SCS <- SCSrank(x=postlogitmat, conf.level=conf.level, alternative=alternative)

OUT<-list(SCS=SCS, Ymat=YMAT, prior=PRIOR, postarr=postarr, postlogodds=postlogodds, 
postdifflogit=postdifflogit, postlogitmat=postlogitmat)

return(c(OUT, LNAM))
}




##############################################

### Simultaneous Wald-type add-x-intervals: ##

##############################################

vcovA <- function(Yvecg, cma){
  ncat <- length(Yvecg)
  nyg <- sum(Yvecg)
  pi <- Yvecg/nyg
  dipi <- diag(1/pi)
  v1c <- matrix(rep(1, ncat), nrow=ncat)
 SIGMAg <- ((cma %*% dipi %*% t(cma))-(cma %*% v1c %*% t(v1c) %*% t(cma)))/nyg
return(SIGMAg)
}


vcovBLg <- function(Yvecg, base=1){
  nbl <- length(Yvecg)-1
  vcm <- matrix(1/Yvecg[base], nrow=nbl, ncol=nbl)
  diag(vcm) <- 1/Yvecg[base] + 1/Yvecg[-base]
  return(vcm)
}


vcovBL <- function(Ymat, base=1){
  # list with groupwise vcov matrices
  matlist <- alply(.data=Ymat, .margins=1, .fun=vcovBLg, base=base) 
  # blockdiagonal matrix with groupwise mats in the diag
  bdmat <- do.call("adiag", args=matlist)
  return(bdmat)
}

oddstrt2logitnames <- function(cmcat, cmgroup, bycomp=TRUE, bychr=" btw ")
{
  
  if(is.null(rownames(cmcat))){rownames(cmcat) <- paste("odds", 1:nrow(cmcat), sep="")}
  if(is.null(rownames(cmgroup))){rownames(cmgroup) <- paste("comp", 1:nrow(cmgroup), sep="")}
  
  compnames <- rownames(cmgroup)
  oddsnames <- rownames(cmcat)
  nodds <- nrow(cmcat)
  ncomp <- nrow(cmgroup)
  
  if(bycomp){
    cnamout <- rep(compnames, each=nodds)
    onamout <- rep(oddsnames, times=ncomp)
    logitnames <- paste(onamout,cnamout,  sep=bychr)
  }else{
    cnamout <- rep(compnames, times=nodds)
    onamout <- rep(oddsnames, each=ncomp)
    logitnames <- paste( onamout,cnamout, sep=bychr)
    
  }
  return(list(logitnames=logitnames, compnames=cnamout, oddsnames=onamout, cmcat=cmcat, cmgroup=cmgroup))
}


waldaddx <- function(Ymat, addx=0, addtype="group",
 cmcat, cmgroup, alternative="two.sided",
 conf.level=0.95, base=1, bycomp=TRUE, bychr=" btw ", ...){
  
# flip the order of a vector (that resulted from as.vector(matrix))
# as if matrix would have been transposed
vecflip <- function(x, nrow, ncol){as.vector(t(matrix(x, nrow=nrow, ncol=ncol)))}

YCG <- datcmcheck(Ymat=Ymat, cmcat=cmcat, cmgroup=cmgroup, ...)

YMAT <- YCG$Ymat
CMCAT <- YCG$cmcat
CMGROUP <- YCG$cmgroup


if(!(is.numeric(addx)|is.integer(addx))){stop("addx must be a single number")}
if(length(addx)>1){warning("Only first entry of addx will be used"); addx<-addx[1]}
if(addx<0 | addx>1){warning("Quantity to be added to cells (addx) is outside [0,1]")}

alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
if(!is.numeric(conf.level)||length(conf.level)!=1||(conf.level>=1 |conf.level<=0)){stop("conf.level must be a single numeric, e.g. 0.95")}

  ADDMAT <- addmat(Ymat=YMAT, alpha=addx, type=addtype)
  Ymataddx <- YMAT + ADDMAT

# replace 0  
  if(any(Ymataddx <=10 * .Machine$double.eps)){
    Ymataddx[Ymataddx <=10 * .Machine$double.eps] <- 0.5
    evcovBL <- vcovBL(Ymataddx, base=base)
    replace0 <- TRUE
  }else{replace0 <- FALSE}
  
  ncat <- ncol(YMAT); cni <- rep(3,ncat)
  if(!is.null(colnames(YMAT))){names(cni)<-colnames(YMAT)}
  
  if(length(base)>1|!(is.numeric(base)|is.integer(base))){stop("base must be a single interger")}
  if(base>ncat|base<=0){stop(paste("There are ", ncat, " categories; base must be in [1,", ncat, ",]"))}
  
  cmBL <- contrMat(cni, type="Dunnett", base=base)
  
  #PsimatBL <- log(Ymataddx) %*% t(cmBL)
  PsimatBL <- cmBL %*% t(log(Ymataddx))
  PsivecBL <- as.vector(PsimatBL)
  evcovBL <- vcovBL(Ymataddx, base=base)
  
  cmcatBL <- CMCAT[,-base]
  
  cmLOR <- kronecker(CMGROUP, cmcatBL) 
  cmLOR

  npara <- nrow(cmLOR)
  
  #round(cbind(PsivecBL,  evcovBL),3)
  # estimate, corresp. covariance mat, and corresp. correlation matrix
  LORest <- cmLOR %*% matrix(PsivecBL, ncol=1)
  LORevcov <- cmLOR %*% evcovBL %*% t(cmLOR)
  LORecorrmat <- cov2cor(LORevcov)
  LORvarest <- diag(LORevcov)
  
  switch(alternative,
         "two.sided"={
           quanti <- qmvnorm(p = conf.level, sigma = LORecorrmat, tail = "both.tails")$quantile
           stderr <- sqrt(LORvarest)
           lCI <- LORest - quanti * stderr
           uCI <- LORest + quanti * stderr
         },
         "less"={
           quanti <- qmvnorm(p = conf.level, sigma = LORecorrmat,  tail = "lower.tail")$quantile
           stderr <- sqrt(LORvarest)
           lCI <- rep(-Inf, npara)
           uCI <- LORest + quanti * stderr
         },
         "greater"={
           quanti <- qmvnorm(p = conf.level, sigma = LORecorrmat,  tail = "upper.tail")$quantile
           stderr <- sqrt(LORvarest)
           lCI <- LORest + quanti * stderr
           uCI <- rep(Inf, npara)
         }
  )
  
  LNAM <- oddstrt2logitnames(cmcat=CMCAT, cmgroup=CMGROUP, bycomp=bycomp, bychr=bychr)
  SCS <- list()
  
  if(bycomp){
    SCS$conf.int <- cbind(lCI, uCI)
    colnames(SCS$conf.int) <- c("lower", "upper")
    rownames(SCS$conf.int) <- LNAM$logitnames
    SCS$estimate <- LORest
    names(SCS$estimate) <- LNAM$logitnames
    SCS$quanti <- quanti
    OUT<-list(SCS=SCS, Ymataddx=Ymataddx, addmat=ADDMAT, LORevcov=LORevcov, LORecorrmat=LORecorrmat, replace0=replace0)
  }
  if(!bycomp){
    flipid <- vecflip(x=1:length(LORest), nrow=nrow(cmcatBL), ncol=nrow(CMGROUP))
    SCS <- list()  
    SCS$conf.int <- cbind(lCI[flipid], uCI[flipid])
    colnames(SCS$conf.int) <- c("lower", "upper")
    rownames(SCS$conf.int) <- LNAM$logitnames
    SCS$estimate <- LORest[flipid]
    names(SCS$estimate) <- LNAM$logitnames
    SCS$quanti <- quanti

    OUT<-list(SCS=SCS, Ymataddx=Ymataddx, addmat=ADDMAT, LORevcov=LORevcov[flipid,flipid], LORecorrmat=LORecorrmat[flipid,flipid], replace0=replace0)
  }
  
return(c(OUT, LNAM))
  
}

  
multinomORci <- function(Ymat, cmcat=NULL, cmgroup=NULL, cimethod="DP", alternative="two.sided",
                                conf.level=0.95, bycomp=FALSE, bychr=" btw ", ...){
args<-list(...)
argnames <- names(args)

cimethod <- match.arg(cimethod, choices=c("DP", "Wald"))

if(is.null(args$base)){BASE <- 1}else{BASE <- args$base}

switch(cimethod,
       "DP"={
         if(any(!argnames %in% c("prior", "BSIM"))){
           warg <- which(!argnames %in% c("prior", "BSIM", "base"))
           warning(paste("Argument(s)", paste(argnames[warg], collapse=", "), " not available for method DP, and thus ignored"))}
         if(is.null(args$prior)){prior<-1}else{prior<-args$prior}
         if(is.null(args$BSIM)){BSIM<-10000}else{BSIM<-args$BSIM}
         RES <- simpostdiri(Ymat=Ymat, prior=prior, cmcat=cmcat, cmgroup=cmgroup,
                            BSIM=BSIM, alternative=alternative, conf.level=conf.level, bycomp=bycomp, bychr=bychr, base=BASE)
         METHODTXT <- paste(signif(conf.level*100, digits=3), "% SCI, ", BSIM, "samples from Dirichlet posterior" )
       },
       "Wald"={
         if(any(!argnames %in% c("addx", "addtype", "base"))){
           warg <- which(!argnames %in% c("addx", "addtype", "base"))
           warning(paste("Arguments", paste(argnames[warg], collapse=", "), " not available for method Wald, and thus ignored"))}
         if(is.null(args$addx)){addx<-0}else{addx<-args$addx}
         if(is.null(args$addtype)){addtype<-"cell"}else{addtype<-args$addtype}        
         RES <- waldaddx(Ymat=Ymat, cmcat=cmcat, cmgroup=cmgroup, alternative=alternative,
                           conf.level=conf.level, bycomp=bycomp, bychr=bychr, addx=addx, addtype=addtype,base=BASE)
         if(addx==0){ADDTXT <- ""}else{ADDTXT <- paste("(after adding", addx, "per", addtype, ")")}
         METHODTXT <- paste(signif(conf.level*100, digits=3), "% Wald-type SCI", ADDTXT)
       })
  
SCI <- cbind(data.frame(odds=RES$oddsname, comp=RES$compname, oddsratio=RES$logitnames#, estimate=RES$SCS$estimate
                        ), RES$SCS$conf.int)

OUT <-list(SCI=SCI, 
           details=c(RES, list(
          cimethod = cimethod,
          METHODTXT = METHODTXT
           )))
class(OUT) <- "multinomORci"

return(OUT)

}


as.data.frame.multinomORci <- function(x, row.names = NULL, optional = FALSE, exp=TRUE, ...){
dargs <- list(...)
dargs$row.names <- row.names
dargs$optional <- optional
DAT <- x$SCI
if(is.null(dargs$row.names)){row.names(DAT) <- NULL}
if(exp){
  DAT$lower <- exp(DAT$lower)
  DAT$upper <- exp(DAT$upper)
}
dargs$x <- DAT
do.call("as.data.frame", args=dargs)
}


print.multinomORci <- function(x, exp=TRUE, ...){
  dargs <- list(...)
  DAT <- as.data.frame(x, exp=exp)
  dargs$x <- DAT

  cat(x$details$METHODTXT, "\n")
  if(exp){cat("Intervals transformed to scale of odds ratios: \n")}else{cat("Intervals on logit-scale: \n")}
  do.call("print", args=dargs)
  
  if(x$details$cimethod == "DP") { cat("with Drichlet prior: \n")
    dargs$x <- x$details$prior
    do.call("print", args=dargs)}
  return(invisible(x))
}

