# wrapper to apply genomic prediction models to an object of class gpData
# models used: regress/asreml and BLR

# date: 2011 - 01 - 16
# author: Valentin Wimmer
# date: 2011 - 05 - 03
# changes: Hans-Juergen Auinger
# date: 2011 - 11 - 21
# changes: return argument by Valentin Wimmer
# date: 2011 - 11 - 30
# changes: Hans-Juergen Auinger
# date: 2012 - 10 - 26
# changes: including R-asreml for BLUP by Hans-Juergen Auinger

gpMod <- function(gpData,model=c("BLUP","BL","BRR"),kin=NULL,predict=FALSE,trait=1,repl=NULL,markerEffects=FALSE,fixed=NULL,random=NULL,...){

  #ASReml <- "package:asreml" %in% search()
  ans <- list()
  model <- match.arg(model)
  m <- NULL
  prediction <- NULL
  if(length(trait) > 1) cat("\n\tNOTE:   The return object will be a list of gpMod-objects!\n")
  if(is.null(fixed)) fixed <- ~1
  if(is.null(random)) {
    random <- " ~ "
    randomFormula <- ~1
  } else {
    randomFormula <- random
    random <- paste(paste(random, collapse=" "), "+ ")
  }
  if(model == "BLUP"){
    if (is.null(kin)){
      if(!gpData$info$codeGeno) stop("Missing object 'kin', or use function codeGeno first!")
      # transposed crossproduct of the genotype matrix is used as relationship to obtain the variance components and mean of RR-BLUP
      if(markerEffects) kin <- tcrossprod(gpData$geno) else kin <- kin(gpData, ret="realized")
    } else {
      if(markerEffects) kin <- tcrossprod(gpData$geno)
    }
  }
  for(i in trait){
    df.trait <- gpData2data.frame(gpData, i, onlyPheno=TRUE, repl=repl)
    # take data from gpData object
    vec.bool <- colnames(df.trait) == "ID" | (colnames(df.trait) %in% dimnames(attr(terms(fixed), "factors"))[[1]]) | colnames(df.trait) %in% dimnames(attr(terms(randomFormula), "factors"))[[1]]
    cnt <- dim(gpData$pheno)[2]
    if(!is.null(gpData$phenoCovars)) cnt <- cnt + dim(gpData$phenoCovars)[2]
    if(i %in% 1:cnt) {
      yName <- dimnames(gpData$pheno)[[2]][as.numeric(i)]
      vec.bool[colnames(df.trait) %in% yName] <- TRUE
    } else {
      vec.bool <- vec.bool | colnames(df.trait) == i
      yName <- i
    }
    df.trait <- df.trait[, vec.bool]
    df.trait <- df.trait[!apply(is.na(df.trait), 1, sum), ]
    kinNames <- unique(df.trait$ID[!df.trait$ID %in% rownames(kin)])
    if(length(kinNames) != 0 & model=="BLUP"){
      df.trait <- df.trait[!df.trait$ID %in% kinNames,]
      warning("Some phenotyped IDs are not in the kinship matrix!\nThese are removed from the analysis")
    }
    kinNew <- kin[unique(df.trait$ID), unique(df.trait$ID)]
    kinTS <- kin[df.trait$ID, df.trait$ID]# expand the matrix to what is needed
    if(model == "BLUP"){
      #if(!ASReml){
        res <- regress(as.formula(paste(yName, paste(fixed, collapse=" "))), Vformula=as.formula(paste(paste(random, collapse=" "), "kinTS")),data=df.trait, identity=TRUE, tol=1e-8,...)
        us <- BLUP(res)$Mean
        genVal <- us[grep("kinTS", names(us))]
        genVal <- genVal[!duplicated(names(genVal))]
        names(genVal) <-  unlist(strsplit(names(genVal), "kinTS."))[(1:length(genVal))*2]
      #} else {
      #  kinTS <- kinNew[dimnames(gpData$pheno)[[1]], dimnames(gpData$pheno)[[1]]]
      #  covM.I <- try(solve(kinTS),TRUE)
        # adding constant to diagonal, if covM is singular
      #  if(class(covM.I)=="try-error"){
      #    warning("Covariance matrix is computationally singular: constant 1e-6 is added to the diagonal elements of the covariance matrix")
      #    covM.I <- solve(kinTS + diag(1e-6,ncol(kinTS)))
       # }
        # print warning in case of numerical problems
       # if(any(covM.I>1e8)) warning("Large >1e8 entries in the inverse covariance matrix")
       # kinGinv <- write.relationshipMatrix(covM.I,file=NULL,type="none",sorting="ASReml",digits=10)
       # attr(kinGinv, "rowNames") <- rownames(kinTS)

       # if("ginverse" %in% names(list(...)))
        #  asrObj <- asreml(as.formula(paste(yName, paste(fixed, collapse=""))), random=as.formula(paste(random, "giv(ID, var=TRUE)")), ginverse=c(list(ID=kinGinv), ginverse), data=df.trait)
        #else
        #  asrObj <- asreml(as.formula(paste(yName, paste(fixed, collapse=""))), random=as.formula(paste(random, "giv(ID, var=TRUE)")), ginverse=list(ID=kinGinv), data=df.trait)

        #genVal <-  asrObj$coefficients$random[substr(names(asrObj$coefficients$random), 1, 8) == "giv(ID)_"]
        #names(genVal) <- substr(names(genVal), 9, nchar(names(genVal)))
      #}
      # genVal <- NULL
      if(markerEffects){
        m <- as.numeric(t(gpData$geno[rownames(kinNew), ]) %*% ginv(kinNew) %*% genVal[rownames(kinNew)])
        names(m) <- colnames(gpData$geno)
      }
      if(predict){
        if(markerEffects){
          prediction <- as.numeric(gpData$geno[!rownames(gpData$geno) %in% rownames(kinNew), ] %*% m)
          names(prediction) <- rownames(gpData$geno)[!rownames(gpData$geno) %in% rownames(kinNew)]
        } else {
          prediction <- gpData$geno %*% t(gpData$geno[rownames(kinNew), ]) %*% ginv(kinNew) %*% genVal[rownames(kinNew)]
          names(prediction) <- rownames(gpData$geno)
          prediction <- prediction[!dimnames(prediction)[[1]] %in% names(genVal)] / mean(prediction[names(genVal)]/genVal)
        }
      }
    }

    if(model=="BL"){
      if(random != " ~ ") stop("Random terms are not supported in model 'BL'!")
      if(dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis have to be done fist.")
      X <- gpData$geno[rownames(gpData$geno) %in% df.trait$ID,]
      y <- df.trait[df.trait$ID %in% rownames(gpData$geno), yName]
      if(fixed != " ~ 1"){
        XF <- as.formula(paste(yName, paste(fixed, collapse=" ")))
        XF <- model.matrix(XF,data=df.trait[df.trait$ID %in% rownames(gpData$geno),],na.action=na.pass)
      } else XF <- NULL
      if(is.null(kin)){
        if(is.null(XF))
          capture.output(res <- BLR(y=y,XL=X,...),file="BLRout.txt")
        else
          capture.output(res <- BLR(y=y,XL=X,XF=XF,...),file="BLRout.txt")
      } else {
        if(is.null(XF))
          capture.output(res <- BLR(y=y,XL=X,GF=list(ID=df.trait$ID,A=kinTS),...),file="BLRout.txt")
        else
          capture.output(res <- BLR(y=y,XL=X,GF=list(ID=df.trait$ID,A=kinTS),XF=XF,...),file="BLRout.txt")
      }
      genVal <- res$yHat - mean(res$yHat)
      names(genVal) <- rownames(X)
      m <- res$bL
      if(predict){
        prediction <- as.numeric(gpData$geno[!rownames(gpData$geno) %in% names(genVal), ] %*% m)
        names(prediction) <- rownames(gpData$geno)[!rownames(gpData$geno) %in% names(genVal)]
      }
    }
    if(model=="BRR"){
      if(random != " ~ ") stop("Random terms are not supported in model 'BRR'!")
      if(dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis has to be done fist.")
      X <-  gpData$geno[rownames(gpData$geno) %in% df.trait$ID,]
      y <- df.trait[df.trait$ID %in% rownames(gpData$geno), yName]
      if(fixed != " ~ 1"){
        XF <- as.formula(paste(yName, paste(fixed, collapse=" ")))
        XF <- model.matrix(XF,data=df.trait[df.trait$ID %in% rownames(gpData$geno),],na.action=na.pass)
      } else XF <- NULL
      if(is.null(kin)){
        if(is.null(XF))
          capture.output(res <- BLR(y=y,XR=X,...),file="BLRout.txt")
        else
          capture.output(res <- BLR(y=y,XR=X,XF=XF,...),file="BLRout.txt")
      } else {
        if(is.null(XF))
          capture.output(res <- BLR(y=y,XR=X,GF=list(ID=df.trait$ID,A=kinTS),...),file="BLRout.txt")
        else
          capture.output(res <- BLR(y=y,XR=X,GF=list(ID=df.trait$ID,A=kinTS),XF=XF,...),file="BLRout.txt")
      }
      genVal <- res$yHat - mean(res$yHat)
      names(genVal) <- rownames(X)
      m <- res$bR
      if(predict){
        prediction <- as.numeric(gpData$geno[!rownames(gpData$geno) %in% names(genVal), ] %*% m)
        names(prediction) <- rownames(gpData$geno)[!rownames(gpData$geno) %in% names(genVal)]
      }
    }

    y <- df.trait[,yName]
    names(y) <- df.trait[,"ID"]
    ret <- list(fit=res,model=model,y=y,g=genVal,prediction=prediction, markerEffects=m,kin=kin)
    class(ret) = "gpMod"
    ans[[i]] <- ret
    names(ans)[length(ans)] <- yName
  }
  if(length(trait) > 1){
   class(ans) <- "gpModList"
   return(ans)
  } else {
    return(ret)
  }
}

summary.gpMod <- function(object,...){
    ans <- list()
    ans$model <- object$model
    if(object$model %in% c("BLUP")) ans$summaryFit <- summary(object$fit)
    if(object$model=="BL") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, varU=object$fit$varU, lambda=object$fit$lambda, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    if(object$model=="BRR") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, varBr=object$fit$varBr, varU=object$fit$varU, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    ans$n <- sum(!is.na(object$y))
    ans$sumNA <- sum(is.na(object$y))
    ans$summaryG <- summary(as.numeric(object$g))
    if(is.null(object$prediction)) ans$summaryP <- NULL else ans$summaryP <- summary(as.numeric(object$prediction))
    class(ans) <- "summary.gpMod"
    ans
}

summary.gpModList <- function(object,...){
  ret <- list()
  for(i in 1:length(object)){
    ret[[i]] <- summary(object[[i]])
  }
  class(ret) <- "summary.gpModList"
  names(ret) <- names(object)
  ret
}

print.summary.gpMod <- function(x,...){
    cat("Object of class 'gpMod' \n")
    cat("Model used:",x$model,"\n")
    cat("Nr. observations ",x$n," \n",sep="")
    cat("Genetic performances: \n")
    cat("  Min.    1st Qu. Median  Mean    3rd Qu. Max    \n")
    cat(format(x$summaryG,width=7,trim=TRUE), "\n",sep=" ")
    if(!is.null(x$summaryP)){
      cat("\nGenetic performances of predicted individuals: \n", sep=" ")
      cat("  Min.    1st Qu. Median  Mean    3rd Qu. Max    \n")
      cat(format(x$summaryP,width=7,trim=TRUE), "\n",sep=" ")
    }
    cat("--\n")
    cat("Model fit \n")
    if(x$model %in% c("BLUP"))
      cat(print(x$summaryFit),"\n")
    else {
      cat("MCMC options: nIter = ",x$summaryFit$nIter,", burnIn = ",x$summaryFit$burnIn,", thin = ",x$summaryFit$thin,"\n",sep="")
      cat("             Posterior mean \n")
      cat("(Intercept) ",x$summaryFit$mu,"\n")
      cat("VarE        ",x$summaryFit$varE,"\n")
      if(!is.null(x$summaryFit$varU)) cat("VarU        ",x$summaryFit$varU,"\n")
      if(x$model=="BL"){
        cat("lambda      ",x$summaryFit$lambda,"\n")
      }
      if(x$model=="BRR"){
        cat("varBr       ",x$summaryFit$varBr,"\n")
      }
    }
}

print.summary.gpModList <- function(x,...) {
  for(i in 1:length(x)){
    cat(paste("\n\tTrait ", names(x)[i], "\n\n\n"))
    print(x[[i]])
    if(i != length(x)) cat("-------------------------------------------------\n\n")
  }
}
