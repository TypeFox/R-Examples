contrastmeans <- function(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", Df, permlist) {
  options(scipen=6)
  
  if (is.null(ctrnames) || ctrnames%in%c("NULL", "")) ctrnames <- NULL
  K <- Kmatrix(model, modelterm)$K
  termsLabel <- rownames(K)
  if (missing(ctrmatrix)) {
    ques <- paste("\nHow many contrasts you want set up for ", sQuote(modelterm), "? ", sep="")
    nctr <- as.integer(readline(ques))
    nrK <- length(termsLabel)
    xnew <- matrix(rep(0, nctr*nrK), ncol=nrK)
    colnames(xnew) <- termsLabel
    xnew <- edit(xnew) 
    ctrmatrix <- xnew[1:nctr, 1:nrK, drop=FALSE]
    rownames(ctrmatrix) <- ctrnames
    if (!(all(rowSums(ctrmatrix)==0)))  {
      cat("\n", "The contrast matrix is:\n\n") 
      print(ctrmatrix)
      cat("\n\n")
      stop("\n", "Please check the row",  sQuote(which(rowSums(ctrmatrix)!=0)), "of the contrast matrix!\n\n")
    }
  }else{
    colnames(ctrmatrix) <- termsLabel
    rownames(ctrmatrix) <- ctrnames
    nctr <- nrow(ctrmatrix)
  }
  if (all((is.null(dim(ctrmatrix)) | dim(ctrmatrix)[1]==1), missing(permlist)))  adj <- "none"
  rK <- ctrmatrix%*%K
  
  mp <- mymodelparm(model)  
  if (all(missing(permlist), missing(Df))){
	  if (mp$df!=0) {
		Df <- mp$df
	  }else{   
		  if (class(model)[1]=="lme") {
			vars <- c(unlist(strsplit(modelterm, "\\:")), modelterm)
			Df <- min(terms(model$fixDF)[vars], na.rm=TRUE)
			mDf <- max(terms(model$fixDF)[vars], na.rm=TRUE)
			if (length(vars) > 2) cat("\n", "Denominator degree of freedom for", 
			  sQuote(modelterm), "and its marginal terms vary between", sQuote(Df), "and", 
			  sQuote(mDf), ".\n","Probabilities will be calculated using", sQuote(Df), "Df.",  "\n")
		  }else if (class(model)[1] == "lmerMod") {
		    vars <- unlist(strsplit(modelterm, "\\:"))
			termlabel <- attr(terms(model),"term.labels")
			for (i in vars) termlabel <- termlabel[grep(i, termlabel)]
			termlabel <- paste(termlabel, collapse="-")
			model.b <- update( model, as.formula(paste(".~. -", termlabel)))
			Df <- getKR(KRmodcomp(model, model.b), "ddf")
      }else stop("You need provide Df for the model!")
	  }
  }

  cm <- rK%*%mp$coef
  vcovm <- mp$vcov
  vcov.contr <- rK %*% tcrossprod(vcovm, rK)
  ses <- sqrt(diag(vcov.contr))
  t.v <- cm/ses
  nr <- nrow(rK)
  dv <- t(1/ses)
  cor.contr <- as.matrix(vcov.contr * (t(dv) %*% dv))  
  
  if (missing(permlist)) {
    t.p.value <- 2*pt(-abs(t.v), Df)
    t.p.value <- p.adjust(t.p.value, adj)
    out.put <- cbind(cm, ses, t.v, Df, t.p.value)
    colnames(out.put) <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
    rownames(out.put) <- ctrnames
    attr(out.put,"Note") <- paste("The p-value is adjusted by", sQuote(adj), "method, if p-value = 0 means p-value < 0.0001.")
    }else{
  	nsim <- length(permlist[[1]])
    tValue <- function(x, rK){
    cm <- rK%*%x$coef
    vcovm <- x$vcov
    vcov.contr <- rK %*% tcrossprod(vcovm, rK)
    ses <- sqrt(diag(vcov.contr))
    t.v <- cm/ses
    return(t.v)	
  }
  
  if (nctr==1) {
    per.p <- (sum(sapply(permlist[[1]], function(x) abs(tValue(x, rK)) > abs(t.v)))+1)/(nsim + 1)
  }else{
    per.p <- (rowSums(sapply(permlist[[1]], function(x) abs(tValue(x, rK)) > abs(t.v)))+1)/(nsim + 1)
  }
  out.put <- cbind(cm, ses, t.v, per.p)
  colnames(out.put) <- c("Estimate", "Std. Error", "t value", "Permuted Pr(>|t|)")
  rownames(out.put) <- ctrnames
  attr(out.put,"Note") <- paste("The permuted p-value is obtained using", sQuote(nsim), "permutations.") 
  }
  return(list("The t tests of the specified contrasts"=round(out.put, 4), "K"=ctrmatrix))
} 
