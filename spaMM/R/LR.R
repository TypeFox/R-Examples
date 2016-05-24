compare.model.structures <- function(object,object2) {
  if (inherits(object,"HLfitlist") || inherits(object2,"HLfitlist")) {
    stop("compare.model.structures does not yet work on HLfitlist objects")
  }
  X1 <- colnames(object$`X.pv`)
  X2 <- colnames(object2$`X.pv`)
  if (length(X1)==0L) {
    REML1 <- NULL ## compatible with both ML or REML tests
  } else REML1 <- (object$APHLs$p_v != object$APHLs$p_bv)
  if (length(X2)==0L) {
    REML2 <- NULL ## idem
  } else REML2 <- (object2$APHLs$p_v != object2$APHLs$p_bv)
  REML <- unique(c(REML1,REML2))
  meth1 <- object$HL
  meth2 <- object2$HL
  if (! identical(meth1,meth2) || length(REML)>1 ) {
    stop("object fitted by different methods cannot be compared")
  }
  if ( ! is.null(X1)) X1 <- sapply(strsplit(X1,':'),function(x) paste(sort(x),collapse=':')) ## JBF 2015/02/23: sort variables in interaction terms before comparison
  if ( ! is.null(X2)) X2 <- sapply(strsplit(X2,':'),function(x) paste(sort(x),collapse=':'))
  dX12 <- setdiff(X1,X2)
  dX21 <- setdiff(X2,X1)
  if (length(dX12)>0 && length(dX21)>0) {
    stop("Non-nested fixed-effect models")
  } else if (length(dX12)>0) {
    Xnest <- "2in1"
  } else if (length(dX21)>0) {
    Xnest <- "1in2"
  } else Xnest <- NULL
  ranterms1 <- attr(object$ZAlist,"ranefs")
  ranterms2 <- attr(object2$ZAlist,"ranefs")
  randist1 <- lapply(object$rand.families,function(v) paste(paste(v)[1:2],collapse="")) ## makes a string from each $family and $link 
  randist2 <- lapply(object2$rand.families,function(v) paste(paste(v)[1:2],collapse="")) ## makes a string from each $family and $link 
  ranterms1 <- paste(ranterms1,randist1) ## joins each term and its distrib
  ranterms2 <- paste(ranterms2,randist2) ## joins each term and its distrib
  dR12 <- setdiff(ranterms1,ranterms2)
  dR21 <- setdiff(ranterms2,ranterms1)
  if (length(dR12)>0 && length(dR21)>0) { 
    stop("Non-nested random-effect models")
  } else if (length(dR12)>0) {
    Rnest <- "2in1"
  } else if (length(dR21)>0) {
    Rnest <- "1in2"
  } else Rnest <- NULL
  nest <- c(Xnest,Rnest)
  unest <- unique(nest)
  if (length(nest)==0L) { ## NULL,NULL
    stop("models to not appear different from each other") 
  } else if (length(unest)==2) {
    stop("Models not nested (opposite nestings for fixed and random terms). ")
  } else {
    df1 <- length(X1)
    df2 <- length(X2)
    if (!is.null(Rnest)) {
      lambda.object <- object$lambda.object
      if (!is.null(lambda.object)) df1 <- df1+nrow(lambda.object$coefficients_lambda)
      cov.mats <- object$cov.mats
      if ( ! is.null(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,nrow))
        df1 <- df1+sum(nrows*(nrows-1)/2)
      }
      lambda.object <- object2$lambda.object
      if (!is.null(lambda.object)) df2 <- df2+nrow(lambda.object$coefficients_lambda)
      cov.mats <- object2$cov.mats
      if ( ! is.null(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,nrow))
        df2 <- df2+sum(nrows*(nrows-1)/2)
      }
    }
    if (unest=="1in2") {
      fullm <- object2
      nullm <- object
      df <- df2-df1
    } else {
      fullm <- object
      nullm <- object2
      df <- df1-df2
    }
    if (length(nest)==2) {
      message("Nested models differing both by in their fixed and in their random terms. ")
      message("Tentatively using marginal likelihood to compare them... ")
      testlik <- "p_v" 
    } else {
      if (is.null(Rnest)) { ## fixed effect test 
        if (REML) {
          ## checking the comparability of REML fits
          if ( ! is.null(fullm$X.Re) ) {
            df.f.Re <-ncol(fullm$X.Re)
          } else df.f.Re <-ncol(fullm$`X.pv`)
          if ( ! is.null(nullm$X.Re) ) {
            df.n.Re <-ncol(nullm$X.Re)
          } else df.n.Re <-ncol(nullm$`X.pv`)
          if ( df.f.Re !=  df.n.Re ) {
            warning("LRT comparing REML fits with different designs is highly suspect")
          }
        }
        testlik <- "p_v"
      } else { ## random effect test
        if ( ! REML) warning("ML fits used to compare different random-effects models...")
        testlik <- "p_bv" ## used in both case, identical to p_v in the non-standard case
        stop("The two models have identical fixed-effect formulas\n and cannot yet be compared properly by this function.")
        ## need to take into account correlations in random slope models for example
      }
    } 
  }
  return(list(fullm=fullm,nullm=nullm,testlik=testlik,df=df))
}

LRT <- function(object,object2,boot.repl=0) { ## compare two HM objects
  info <- compare.model.structures(object,object2)
  nullm <- info$nullm; fullm <- info$fullm; testlik <- info$testlik;df <- info$df
  LRTori <- 2*(fullm$APHLs[[testlik]]-nullm$APHLs[[testlik]])
  pvalue <- 1-pchisq(LRTori,df=df) ## but not valid for testing null components of variance
  resu <- list(nullfit=nullm,fullfit=fullm,basicLRT = data.frame(LR2=LRTori,df=df,pvalue=pvalue)) ## format appropriate for more tests  
  if (boot.repl>0) {
    simbData <- nullm$data
    if (tolower(nullm$family$family)=="binomial") {
      form <- attr(nullm$predictor,"oriFormula") ## this must exists...  
      if (is.null(form)) {
        mess <- pastefrom("a 'predictor' object must have an 'oriFormula' member.",prefix="(!) From ")
        stop(mess)
      }
      exprL <- as.character(form[[2]][[2]]) 
      exprR <- as.character(form[[2]][[3]]) 
    }
    if (boot.repl<100) print("It is recommended to set boot.repl>=100 for Bartlett correction",quote=FALSE)
    aslistfull <- as.list(getCallHL(fullm)) 
    ## problem is for corrHLfit etc this is the call of the final HLfit call with $processed and a lot of missing original arguments  
    aslistfull$processed <- NULL ## may capture bugs 
    aslistnull <- as.list(getCallHL(nullm))
    aslistnull$processed <- NULL ## may capture bugs
    computeBootRepl <- function() {
      ## draw sample
      newy <- simulate(nullm)  ## only a vector of response value ## cannot simulate all samples in one block since some may not be analyzable  
      if (tolower(nullm$family$family)=="binomial") {
        ## c'est bouseux: soit j'ai (pos, neg) et le remplacement est possible
        ##    soit j'ai (pos,ntot -pos) et le 2e remplacment n'est pas poss (et pas necess)
        ##    aussi (ntot - pos, pos) ...
        ## would be simple if always ntot-pos, but how to control this ? 
        if (length(exprL)==1L) simbData[[exprL]] <- newy 
        if (length(exprR)==1L) simbData[[exprR]] <- nullm$weights - newy                    
        ## if (length(exprR)! =1) exprRdoes not correspond to a column in the data;frmae so there is no column to replace                     
      } else {simbData[[as.character(nullm$predictor[[2]])]] <- newy}
      ## analyze under both models
      aslistfull$data <- simbData
      aslistnull$data <- simbData
      fullfit <- (eval(as.call(aslistfull)))
      if (inherits(fullfit,"try-error")) return(fullfit) ## eg separation in binomial models
      nullfit <- try(eval(as.call(aslistnull)))
      if (inherits(nullfit,"try-error")) return(nullfit) 
      ## return pair of likelihoods
      if (inherits(fullfit,"HLfitlist")) {
        return(c(attr(fullfit,"APHLs")[[testlik]],attr(nullfit,"APHLs")[[testlik]]))
      } else return(c(fullfit$APHLs[[testlik]],nullfit$APHLs[[testlik]]))
    }
    bootLs<-matrix(,nrow=boot.repl,ncol=2) 
    colnames(bootLs) <- paste(c("full.","null."),testlik,sep="")
    msg <- "bootstrap replicates: "
    msglength <- nchar(msg)
    cat(msg)
    t0 <- proc.time()["user.self"]
    for (ii in seq_len(boot.repl)) {
      locitError <- 0
      repeat { ## for each ii!
        bootrep <- (computeBootRepl())
        if ( ! inherits(bootrep,"try-error")) { 
          bootLs[ii,] <- bootrep
          break ## replicate performed, breaks the repeat
        } else { ## there was one error
          locitError <- locitError + 1
          if (locitError>10) { ## to avoid an infinite loop
            stop("Analysis of bootstrap samples fails repeatedly. Maybe no statistical information in them ?")
          } ## otherwise repeat!
        }
      } 
      tused <- proc.time()["user.self"]-t0
      ttotal <- tused* boot.repl/ii
      if (interactive()) {
        for (bidon in 1:msglength) cat("\b")
        msg <- paste("Estimated time remaining for bootstrap: ",signif(ttotal-tused,2)," s.",sep="")
        msglength <- nchar(msg)
        cat(msg)
      } else {
        cat(ii);cat(" ")
        if ((ii %% 40)==0L) cat("\n")
      }
    }
    cat("\n")
    bootdL <- bootLs[,1]-bootLs[,2]
    meanbootLRT <- 2*mean(bootdL)
    resu <- c(resu,list(rawBootLRT = data.frame(LR2=LRTori,df=df,pvalue=(1+sum(bootdL>=LRTori/2))/(boot.repl+1)))) ## format appropriate for more tests  
    LRTcorr <- LRTori*df/meanbootLRT
    resu <- c(resu,list(BartBootLRT = data.frame(LR2=LRTcorr,df=df,pvalue=1-pchisq(LRTcorr,df=df)))) ## format appropriate for more tests  
    bootInfo <- list(meanbootLRT = meanbootLRT,bootreps = bootLs)
    resu <- c(resu,list(bootInfo=bootInfo)) ## keeps the sublist structure, which is not compatible with hglmjob.R...  
  }
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}

## anova treated as alias for LRT
anova.HLfit <- function(object,object2,...) {
  LRT(object,object2,...)
}

