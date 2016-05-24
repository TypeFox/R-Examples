Waldtest.Rm <- function(object, splitcr = "median"){
  # performs item-based Wald test (Fischer & Molenaar, p.90)
  # object... object of class RM
  # splitcr... splitting criterion for LR-groups. "median" to a median raw score split,
  #            "mean" corobjectponds to the mean raw score split.
  #            optionally also a vector of length n for group split can be submitted.

  call<-match.call()

  spl.gr<-NULL

  X.original<-object$X
  if (length(splitcr)>1 && is.character(splitcr)){    # if splitcr is character vector, treated as factor
     splitcr<-as.factor(splitcr)
  }
  if (is.factor(splitcr)){
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
      rsum.all<-rowSums(X,na.rm=TRUE)
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
      #removed rh 2010-12-17
      #cat("Warning message: Persons with median raw scores are assigned to the lower raw score group!\n")
      X<-object$X
      # calculates index for NA groups
      # from person.parameter.eRm
        dichX <- ifelse(is.na(X),1,0)
        strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
        gmemb <- as.vector(data.matrix(data.frame(strdata)))
      gindx<-unique(gmemb)
      rsum.all<-rowSums(X,na.rm=TRUE)
      grmed<-tapply(rsum.all,gmemb,median)      #sorted
      ngr<-table(gmemb)                         #sorted
      m.all<-rep(grmed,ngr)                     #sorted,expanded
      rsum.all<-rsum.all[order(gmemb)]
      spl<-ifelse(rsum.all<=m.all,1,2)
      splitcr<-spl
      object$X<-X[order(gmemb),]
    }
  }


  if (is.numeric(splitcr)){
    spl.nam<-deparse(substitute(splitcr))
    if (length(table(splitcr)) > 2) stop("Dichotomous person split required!")
    if (length(splitcr) != dim(object$X)[1]) {
      stop("Mismatch between length of split vector and number of persons!")
    } else {
      rvind <- splitcr
      Xlist <- by(object$X,rvind, function(x) x)
      names(Xlist) <- as.list(sort(unique(splitcr)))
      if(is.null(spl.gr)){
        spl.lev<-names(Xlist)
        spl.gr<-paste(spl.nam,spl.lev,sep=" ")
      }
    }}

  if (!is.numeric(splitcr)) {
    if (splitcr=="median") {                                   #median split
      rv <- apply(object$X,1,sum,na.rm=TRUE)
      rvsplit <- median(rv)
      rvind <- rep(0,length(rv))
      rvind[rv > rvsplit] <- 1                                 #group with high raw score object
      Xlist <- by(object$X,rvind,function(x) x)
      names(Xlist) <- list("low","high")
      }

    if (splitcr=="mean") {                                     #mean split
      rv <- apply(object$X,1,sum,na.rm=TRUE)
      rvsplit <- mean(rv)
      rvind <- rep(0,length(rv))
      rvind[rv > rvsplit] <- 1                                 #group with highraw scoobject
      Xlist <- by(object$X,rvind,function(x) x)
      names(Xlist) <- list("low","high")
      }

  }

  del.pos.l <- lapply(Xlist, function(x) {
                      it.sub <- datcheck.LRtest(x,object$X,object$model)  #items to be removed within subgroup
                      })

  del.pos <- unique(unlist(del.pos.l))
  if ((length(del.pos)) >= (dim(object$X)[2]-1)) {
    stop("\nNo items with appropriate response patterns left to perform Wald-test!\n")
  }

  if(length(del.pos) > 0){
    warning(paste0(
      "\n", 
      prettyPaste("The following items were excluded due to inappropriate response patterns within subgroups:"),
      "\n",
      paste(colnames(object$X)[del.pos], collapse=" "),
      "\n\n",
      prettyPaste("Subgroup models are estimated without these items!")
    ), immediate.=TRUE)
  }

  if(length(del.pos) > 0){
    X.el <- object$X[,-(del.pos)]
  } else {
    X.el <- object$X
  }
  Xlist.n <- by(X.el,rvind,function(y) y)
  names(Xlist.n) <- names(Xlist)

  if (object$model=="RM") {
         likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                                 objectg <- RM(x)
                                 parg <- objectg$etapar
                                 seg <- objectg$se.eta
                                 list(parg,seg,objectg$betapar,objectg$se.beta)
                                 })
         }
  if (object$model=="PCM") {
         likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                                 objectg <- PCM(x)
                                 parg <- objectg$etapar
                                 seg <- objectg$se.eta
                                 list(parg,seg,objectg$betapar,objectg$se.beta)
                                 })
         }
  if (object$model=="RSM") {
         likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                                 objectg <- RSM(x)
                                 parg <- objectg$etapar
                                 seg <- objectg$se.eta
                                 list(parg,seg,objectg$betapar,objectg$se.beta)
                                 })
         }


  betapar1 <- likpar[3,][[1]]
  beta1.se <- likpar[4,][[1]]
  betapar2 <- likpar[3,][[2]]
  beta2.se <- likpar[4,][[2]]
  num <- (betapar1-betapar2)
  denom <- sqrt(beta1.se^2 + beta2.se^2)
  W.i <- num/denom
  pvalues <- (1-pnorm(abs(W.i)))*2

  coef.table <- cbind(W.i,pvalues)
  dimnames(coef.table) <- list(names(betapar1),c("z-statistic","p-value"))

  result <- list(coef.table=coef.table,betapar1=betapar1,se.beta1=beta1.se,betapar2=betapar2,
  se.beta2=beta2.se, spl.gr=spl.gr, call=call, it.ex = del.pos)
  class(result) <- "wald"
  result

}
