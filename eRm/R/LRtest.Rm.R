`LRtest.Rm` <-
function(object, splitcr = "median", se = TRUE)
{
# performs Andersen LR-test
# object... object of class RM
# splitcr... splitting criterion for LR-groups. "all.r" corresponds to a complete
#            raw score split (r=1,...,k-1), "median" to a median raw score split,
#            "mean" corresponds to the mean raw score split.
#            optionally also a vector of length n for group split can be submitted.
# se...whether standard errors should be computed


call<-match.call()

spl.gr<-NULL

X.original<-object$X
if((length(splitcr) > 1) & is.character(splitcr)){    # if splitcr is character vector, treated as factor
  splitcr<-as.factor(splitcr)
}
if(is.factor(splitcr)){
   spl.nam<-deparse(substitute(splitcr))
   spl.lev<-levels(splitcr)
   spl.gr<-paste(spl.nam,spl.lev,sep=" ")
   splitcr<-unclass(splitcr)
}

numsplit<-is.numeric(splitcr)
if (any(is.na(object$X))) {
  if (!numsplit && splitcr=="mean") {                                   #mean split
    spl.gr<-c("Raw Scores < Mean", "Raw Scores >= Mean")
    X<-object$X
    # calculates index for NA groups
    # from person.parameter.eRm
      dichX <- ifelse(is.na(X),1,0)
      strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
      gmemb <- as.vector(data.matrix(data.frame(strdata)))
    gindx<-unique(gmemb)
    rsum.all<-rowSums(X,na.rm=T)
    grmeans<-tapply(rsum.all,gmemb,mean)      #sorted
    ngr<-table(gmemb)                         #sorted
    m.all<-rep(grmeans,ngr)                   #sorted,expanded
    rsum.all<-rsum.all[order(gmemb)]
    spl<-ifelse(rsum.all<m.all,1,2)
    splitcr<-spl
    object$X<-X[order(gmemb),]
  }
  if (!numsplit && splitcr=="median") {                                   #median split
    spl.gr<-c("Raw Scores <= Median", "Raw Scores > Median")
    # cat("Warning message: Persons with median raw scores are assigned to the lower raw score group!\n")
    X<-object$X
    # calculates index for NA groups
    # from person.parameter.eRm
      dichX <- ifelse(is.na(X),1,0)
      strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
      gmemb <- as.vector(data.matrix(data.frame(strdata)))
    gindx<-unique(gmemb)
    rsum.all<-rowSums(X,na.rm=T)
    grmed<-tapply(rsum.all,gmemb,median)      #sorted
    ngr<-table(gmemb)                         #sorted
    m.all<-rep(grmed,ngr)                     #sorted,expanded
    rsum.all<-rsum.all[order(gmemb)]
    spl<-ifelse(rsum.all<=m.all,1,2)
    splitcr<-spl
    object$X<-X[order(gmemb),]
  }
}

if (!is.numeric(splitcr)) {
  if (splitcr=="all.r") {                               #full raw score split   ### begin MjM 2012-03-18
    rvind <- rowSums(object$X, na.rm=TRUE)              #person raw scoobject
    excl_0_k <- (rvind > 0) & (rvind < sum(apply(object$X, 2, max, na.rm=T)))
    Xlist <- by(object$X[excl_0_k,], rvind[excl_0_k], function(x) x)
    names(Xlist) <- as.list(paste("Raw Score =", sort(unique(rvind[excl_0_k]))))
    spl.gr <- unlist(names(Xlist))
  }                                                                             ### end MjM 2012-03-18

  if (splitcr=="median") {                                   #median split
    spl.gr<-c("Raw Scores <= Median", "Raw Scores > Median")
    #removed rh 2010-12-17
    #cat("Warning message: Persons with median raw scores are assigned to the lower raw score group!\n")
    rv <- apply(object$X,1,sum,na.rm=TRUE)
    rvsplit <- median(rv)
    rvind <- rep(0,length(rv))
    rvind[rv > rvsplit] <- 1                                 #group with highraw scoobject
    Xlist <- by(object$X,rvind,function(x) x)
    names(Xlist) <- list("low","high")
  }

  if (splitcr=="mean") {                                     #mean split
    spl.gr<-c("Raw Scores < Mean", "Raw Scores >= Mean")
    rv <- apply(object$X,1,sum,na.rm=TRUE)
    rvsplit <- mean(rv)
    rvind <- rep(0,length(rv))
    rvind[rv > rvsplit] <- 1                                 #group with highraw scoobject
    Xlist <- by(object$X,rvind,function(x) x)
    names(Xlist) <- list("low","high")
    }
}

if (is.numeric(splitcr)) {                                 #manual raw score split
  spl.nam<-deparse(substitute(splitcr))
  if (length(splitcr)!=dim(object$X)[1]){
    stop("Mismatch between length of split vector and number of persons!")
  } else {
    rvind <- splitcr
    Xlist <- by(object$X,rvind, function(x) x)
    names(Xlist) <- as.list(sort(unique(splitcr)))
    if(is.null(spl.gr)){
      spl.lev<-names(Xlist)
      spl.gr<-paste(spl.nam,spl.lev,sep=" ")
    }
  }
}

#----------item to be deleted---------------
del.pos.l <- lapply(Xlist, function(x) {
                    it.sub <- datcheck.LRtest(x,object$X,object$model)  #items to be removed within subgroup
                    })

del.pos <- unique(unlist(del.pos.l))
if (length(del.pos) >= (ncol(object$X)-1)) {
  stop("\nNo items with appropriate response patterns left to perform LR-test!\n")
}

if(length(del.pos) > 0){                                                        ### begin MjM 2013-01-27
  warning(paste0(
    "\n", 
    prettyPaste("The following items were excluded due to inappropriate response patterns within subgroups:"),
    "\n",
    paste(colnames(object$X)[del.pos], collapse=" "),
    "\n\n",
    prettyPaste("Full and subgroup models are estimated without these items!")
  ), immediate.=TRUE)
}                                                                               ### end MjM 2013-01-27


if (length(del.pos) > 0) {
  X.el <- object$X[,-(del.pos)]
} else {
  X.el <- object$X
}

if(ifelse(length(splitcr) == 1, splitcr != "all.r", TRUE)){   ### begin MjM 2012-03-18   # for all cases except "all.r"
  Xlist.n <- by(X.el, rvind, function(y) y)
  names(Xlist.n) <- names(Xlist)
  if (length(del.pos) > 0) Xlist.n <- c(Xlist.n,list(X.el)) # X.el added since we must refit whole group without del.pos items
} else {
  Xlist.n <- by(X.el[excl_0_k,], rvind[excl_0_k], function(y) y)
  names(Xlist.n) <- names(Xlist)
  Xlist.n <- c(Xlist.n,list(X.el[excl_0_k,])) # X.el added since we must refit whole group without del.pos items
}                         ### end MjM 2012-03-18

if (object$model=="RM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                              # betalab <- colnames(objectg$X)
                               list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta,outobj=objectg)   # rh outobj added
                               ###list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta)   # rh outobj added
                               })
       }
if (object$model=="PCM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- PCM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta,outobj=objectg)   # rh outobj added
                               ###list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta)   # rh outobj added
                               })
       }
if (object$model=="RSM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RSM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta,outobj=objectg)   # rh outobj added
                               ###list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta)   # rh outobj added
                               })
       }

## extract fitted splitgroup models  # rh 02-05-2010
if(ifelse(length(splitcr) == 1, splitcr != "all.r", TRUE)){   ### begin MjM 2012-03-18
  fitobj <- likpar[6, 1:length(unique(rvind))]
} else {
  fitobj <- likpar[6, 1:length(unique(rvind[excl_0_k]))]
}                         ### end MjM 2012-03-18
likpar <- likpar[-6,]

if((length(del.pos) > 0) | ifelse(length(splitcr) == 1, splitcr == "all.r", FALSE)) {                  #re-estimate full model   ### MjM 2012-03-18
  pos <- length(Xlist.n)                    #position of the full model
  loglik.all <- likpar[1,pos][[1]]          #loglik full model
  etapar.all <- rep(0,likpar[2,pos])         #etapar full model (filled with 0 for df computation)
  likpar <- likpar[,-pos]
  Xlist.n <- Xlist.n[-pos]
} else {
  loglik.all <- object$loglik
  etapar.all <- object$etapar
}

loglikg <- sum(unlist(likpar[1,]))                    #sum of likelihood value for subgroups
LR <- 2*(abs(loglikg-loglik.all))                  #LR value
df = sum(unlist(likpar[2,]))-(length(etapar.all))  #final degrees of freedom
pvalue <- 1 - pchisq(LR, df)                             #pvalue

betalist <- likpar[3,]                                #organizing betalist


result <- list(X=X.original, X.list=Xlist.n, model=object$model,LR=LR,
               df=df, pvalue=pvalue, likgroup=unlist(likpar[1,],use.names=FALSE),
               betalist=betalist, etalist=likpar[4,],selist=likpar[5,], spl.gr=spl.gr, call=call, fitobj=fitobj)  ## rh fitobj added
class(result) <- "LR"

return(result)

}
