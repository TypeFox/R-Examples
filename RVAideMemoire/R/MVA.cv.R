# pls : plsr,cppls
# MASS : lda,qda

print.MVA.cv <- function(x,...) {
  cat("\n        Cross validation\n\n")
  cat(paste0("Model:"),x$model,"\n")
  cat(paste0(x$k,"-fold validation\n"))
  if (x$repet>1) {cat(paste0("Validation repeated ",x$repet," times\n"))}
  ncomp <- if (x$type!="qual2") {
    x$ncomp
  } else {
    paste(x$ncomp,"+",length(x$models2.list[[1]]$prior)-1)
  }
  cat(paste0(x$repet*x$k," ",ifelse(x$type=="qual2","couples of ",""),"submodels generated (",
    ncomp," components)\n"))
  if (x$model %in% c("LDA","QDA","PLS-DA/LDA","PLS-DA/QDA","PPLS-DA/LDA","PPLS-DA/QDA")) {
    cat(paste0("\nClassification criterion:"),x$crit.DA,"\n")
  }
  cat("\n")
  if (x$type=="quant") {
    if (ncol(x$RMSEP)==1) {
	cat(paste0("Mean (standard error) RMSEP: ",signif(mean(x$RMSEP),3)," (",signif(se(x$RMSEP),2),")\n"))
	cat(paste0("Mean (standard error) Q2: ",signif(mean(x$Q2),3)," (",signif(se(x$Q2),2),")\n"))
    } else {
	cat("Mean (standard error) RMSEP:\n")
	to.print1 <- paste0(signif(colMeans(x$RMSEP),3),"(",signif(apply(x$RMSEP,2,se),2),")")
	names(to.print1) <- colnames(x$RMSEP)
	print(to.print1,quote=FALSE)
	cat("Mean (standard error) Q2:\n")
	to.print2 <- paste0(signif(colMeans(x$Q2),3),"(",signif(apply(x$Q2,2,se),2),")")
	names(to.print2) <- colnames(x$Q2)
	print(to.print2,quote=FALSE)
    }
  } else {
    cat(paste0("Mean (standard error) number of misclassifications (%): ",signif(100*mean(x$NMC),3)," (",
	signif(100*se(x$NMC),2),")\n"))
  }
  cat("\n")
}

MVA.cv <- function(X,Y,repet=10,k=7,ncomp=5,model=c("PLSR","CPPLS","PLS-DA","PPLS-DA","LDA","QDA","PLS-DA/LDA",
  "PLS-DA/QDA","PPLS-DA/LDA","PPLS-DA/QDA"),lower=0.5,upper=0.5,Y.add=NULL,weights=rep(1,nrow(X)),set.prior=FALSE,
  crit.DA=c("plug-in","predictive","debiased"),...) {
  model <- match.arg(model)
  crit.DA <- match.arg(crit.DA)
  type <- if (model %in% c("PLSR","CPPLS")) {"quant"} else {
    if (model %in% c("PLS-DA","PPLS-DA","LDA","QDA")) {"qual1"} else {"qual2"}
  }
  if (type=="quant" & is.factor(Y)) {
    if (model=="PLSR") {
	warning("'model' re-set to 'PLS-DA'")
    } else {
	warning("'model' re-set to 'PPLS-DA'")
    }
    model <- ifelse(model=="PLSR","PLS-DA","PPLS-DA")
    type <- "qual1"
  }
  if (is.factor(Y)) {
    Yfac <- Y
    lev <- paste0("Y.",levels(Y))
    Y <- I(model.matrix(~Y-1))
  }
  X <- as.matrix(as.data.frame(X))
  Y <- as.matrix(as.data.frame(Y))
  if (is.factor(Y)) {colnames(Y) <- lev}
  prior <- NULL
  if (set.prior) {
    mweights <- tapply(weights,Yfac,mean)
    prior <- as.vector(mweights/sum(mweights))
  }
  fun <- switch(type,quant=MVA.cv.quant,qual1=MVA.cv.qual1,qual2=MVA.cv.qual2)
  res <- if (type=="quant") {fun(X,Y,repet,k,ncomp,model,lower,upper,Y.add,weights,...)} else
    if (type=="qual1") {fun(X,Y,Yfac,groups=levels(Yfac),repet,k,ncomp,model,lower,upper,Y.add,weights,prior,
    crit.DA=crit.DA,...)} else
   {fun(X,Y,Yfac,groups=levels(Yfac),repet,k,ncomp,model,lower,upper,Y.add,weights,prior,crit.DA=crit.DA,...)}
  class(res) <- c("list","MVA.cv")
  return(res)
}

MVA.cv.quant <- function(X,Y,repet,k,ncomp,model,lower,upper,Y.add,weights,...) {
  whole.set <- as.data.frame(cbind(weights,Y.add,Y,X))
  rownames(whole.set) <- 1:nrow(whole.set)
  col.Yadd <- if (!is.null(Y.add)) {2:(2+ncol(Y.add)-1)} else {NULL}
  col.Y <- if (!is.null(Y.add)) {(2+ncol(Y.add)):(2+ncol(Y.add)+ncol(Y)-1)} else {2:(2+ncol(Y)-1)}
  col.X <- if (!is.null(Y.add)) {(2+ncol(Y.add)+ncol(Y)):ncol(whole.set)} else {(2+ncol(Y)):ncol(whole.set)}
  test.sets.list.repet <- list()
  length(test.sets.list.repet) <- repet
  test.sets.list.repet <- lapply(test.sets.list.repet,function(x) {split(whole.set,sample(gl(k,1,nrow(whole.set))))})
  models.list <- list()
  length(models.list) <- repet*k
  names(models.list) <- paste(rep(1:repet,each=k),rep(1:k,repet),sep=":")
  N <- nrow(Y)
  TSS <- apply(Y,2,function(x) {sum((x-mean(x))^2)})
  RMSEP <- matrix(0,nrow=repet,ncol=ncol(Y),dimnames=list(1:repet,colnames(Y)))
  Q2 <- matrix(0,nrow=repet,ncol=ncol(Y),dimnames=list(1:repet,colnames(Y)))
  for (i in 1:repet) {
    test.sets.list.k <- test.sets.list.repet[[i]]
    pred <- matrix(0,nrow=nrow(whole.set),ncol=ncol(Y),dimnames=list(1:nrow(whole.set),colnames(Y)))
    for (j in 1:k) {
	test.set <- test.sets.list.k[[j]]
	test.set.X <- as.matrix(as.data.frame(test.set[,col.X]))
	train.set <- whole.set[-as.numeric(rownames(test.set)),]
	train.set.X <- as.matrix(as.data.frame(train.set[,col.X]))
	train.set.weights <- train.set$weights
	train.set.Yadd <- train.set[,col.Yadd]
	train.set.Y <- as.matrix(as.data.frame(train.set[,col.Y]))
	if (j==1) {
	  nmax <- min(c(nrow(train.set)-max(unlist(lapply(test.sets.list.k,nrow))),ncol(X)+1))
	  if (ncomp>=nmax) {
	    ncomp <- nmax-1
	    warning(paste0("'ncomp' re-set to ",ncomp))
	  }
	}
	model.k <- if (model=="PLSR") {
	  pls::plsr(train.set.Y~train.set.X,ncomp=ncomp,...)
	} else {
	  if (!is.null(Y.add)) {
	    pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,Y.add=train.set.Yadd,
		weights=train.set.weights,...)
	  } else {
	    pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,weights=train.set.weights,...)
	  }
	}
	models.list[[i*k-(k-j)]] <- model.k
	pred[as.numeric(rownames(test.set)),] <- predict(model.k,newdata=test.set.X,ncomp=ncomp)
    }
    PRESS <- colSums((Y-pred)^2)
    RMSEP[i,] <- sqrt(PRESS/N)
    Q2[i,] <- 1-PRESS/TSS
  }
  return(list(model=model,type="quant",repet=repet,k=k,ncomp=ncomp,models.list=models.list,RMSEP=RMSEP,Q2=Q2))
}

MVA.cv.qual1 <- function(X,Y,Yfac,groups,repet,k,ncomp,model,lower,upper,Y.add,weights,prior,crit.DA,...) {
  whole.set <- as.data.frame(cbind(weights,Y.add,Y,X))
  rownames(whole.set) <- 1:nrow(whole.set)
  trueclass <- apply(Y,1,function(x) {colnames(Y)[which(x==1)]})
  levels(Yfac) <- unique(trueclass)
  col.Yadd <- if (!is.null(Y.add)) {2:(2+ncol(Y.add)-1)} else {NULL}
  col.Y <- if (!is.null(Y.add)) {(2+ncol(Y.add)):(2+ncol(Y.add)+ncol(Y)-1)} else {2:(2+ncol(Y)-1)}
  col.X <- if (!is.null(Y.add)) {(2+ncol(Y.add)+ncol(Y)):ncol(whole.set)} else {(2+ncol(Y)):ncol(whole.set)}
  test.sets.list.repet <- list()
  length(test.sets.list.repet) <- repet
  test.sets.list.repet <- lapply(test.sets.list.repet,function(x) {splitf(whole.set,factor(trueclass),k)})
  test.sets.repet.trueclass <- lapply(test.sets.list.repet,function(x) {lapply(x,function(y) {trueclass[as.numeric(rownames(y))]})})
  models.list <- list()
  length(models.list) <- repet*k
  names(models.list) <- paste(rep(1:repet,each=k),rep(1:k,repet),sep=":")
  NMC <- numeric(repet)
  for (i in 1:repet) {
    test.sets.list.k <- test.sets.list.repet[[i]]
    pred <- character(nrow(Y))
    for (j in 1:k) {
	test.set <- test.sets.list.k[[j]]
	test.set.X <- as.matrix(as.data.frame(test.set[,col.X]))
	train.set <- whole.set[-as.numeric(rownames(test.set)),]
	train.set.weights <- train.set$weights
	train.set.Yadd <- train.set[,col.Yadd]
	train.set.Y <- as.matrix(as.data.frame(train.set[,col.Y]))
	train.set.Yfac <- Yfac[-as.numeric(rownames(test.set))]
	train.set.X <- as.matrix(as.data.frame(train.set[,col.X]))
	train.set.trueclass <- trueclass[-as.numeric(rownames(test.set))]
	if (j==1) {
	  if (model %in% c("PLS-DA","PPLS-DA")) {
	    nmax <- min(c(nrow(train.set)-max(unlist(lapply(test.sets.list.k,nrow))),ncol(X)+1))
	    if (ncomp>=nmax) {
		ncomp <- nmax-1
		warning(paste0("'ncomp' re-set to ",ncomp))
	    }
	  } else {
	    ncomp <- nlevels(Yfac)-1
	  }
	}
	model.k <- if (model=="LDA") {
	  if (!is.null(prior)) {
	    MASS::lda(as.matrix(train.set.X),train.set.Yfac,prior=prior,tol=1.0e-8)
	  } else {
	    MASS::lda(as.matrix(train.set.X),train.set.Yfac,tol=1.0e-8)
	  }
	} else if (model=="QDA") {
	  if (!is.null(prior)) {
	    MASS::qda(as.matrix(train.set.X),train.set.Yfac,prior=prior,tol=1.0e-8)
	  } else {
	    MASS::qda(as.matrix(train.set.X),train.set.Yfac,tol=1.0e-8)
	  }
	} else if (model=="PLS-DA") {
	  pls::plsr(train.set.Y~train.set.X,ncomp=ncomp,...)
	} else {
	  if (!is.null(Y.add)) {
	    pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,Y.add=train.set.Yadd,
		weights=train.set.weights,...)
	  } else {
	    pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,weights=train.set.weights,...)
	  }
	}
	models.list[[i*k-(k-j)]] <- model.k
	if (model %in% c("LDA","QDA")) {
	  pred.lev <- as.character(predict(model.k,newdata=test.set.X,method=crit.DA)$class)
	} else {
	  pred.dummy <- predict(model.k,newdata=test.set.X,ncomp=ncomp)
	  pred.lev <- apply(pred.dummy,1,function(x) {colnames(Y)[which.max(x)]})
	}
	pred[as.numeric(rownames(test.set))] <- pred.lev
    }
    pred.correct <- pred==trueclass
    rate <- 1-sum(pred.correct)/length(pred.correct)
    NMC[i] <- rate
  }
  return(list(model=model,type="qual1",repet=repet,k=k,ncomp=ncomp,crit.DA=crit.DA,groups=groups,
    models.list=models.list,NMC=NMC))
}

MVA.cv.qual2 <- function(X,Y,Yfac,groups,repet,k,ncomp,model,lower,upper,Y.add,weights,prior,crit.DA,...) {
  whole.set <- as.data.frame(cbind(weights,Y.add,Y,X))
  rownames(whole.set) <- 1:nrow(whole.set)
  trueclass <- as.character(Yfac)
  col.Yadd <- if (!is.null(Y.add)) {2:(2+ncol(Y.add)-1)} else {NULL}
  col.Y <- if (!is.null(Y.add)) {(2+ncol(Y.add)):(2+ncol(Y.add)+ncol(Y)-1)} else {2:(2+ncol(Y)-1)}
  col.X <- if (!is.null(Y.add)) {(2+ncol(Y.add)+ncol(Y)):ncol(whole.set)} else {(2+ncol(Y)):ncol(whole.set)}
  test.sets.list.repet <- list()
  length(test.sets.list.repet) <- repet
  test.sets.list.repet <- lapply(test.sets.list.repet,function(x) {splitf(whole.set,Yfac,k)})
  test.sets.repet.trueclass <- lapply(test.sets.list.repet,function(x) {lapply(x,function(y) {trueclass[as.numeric(rownames(y))]})})
  models.list1 <- models.list2 <- list()
  length(models.list1) <- length(models.list2) <- repet*k
  names(models.list1) <- names(models.list2) <- paste(rep(1:repet,each=k),rep(1:k,repet),sep=":")
  NMC <- numeric(repet)
  for (i in 1:repet) {
    test.sets.list.k <- test.sets.list.repet[[i]]
    pred <- character(nrow(Y))
    for (j in 1:k) {
	test.set <- test.sets.list.k[[j]]
	test.set.X <- as.matrix(as.data.frame(test.set[,col.X]))
	train.set <- whole.set[-as.numeric(rownames(test.set)),]
	train.set.weights <- train.set$weights
	train.set.Yadd <- train.set[,col.Yadd]
	train.set.Y <- as.matrix(as.data.frame(train.set[,col.Y]))
	train.set.X <- as.matrix(as.data.frame(train.set[,col.X]))
	train.set.trueclass <- trueclass[-as.numeric(rownames(test.set))]
	train.set.Yfac <- Yfac[-as.numeric(rownames(test.set))]
	if (j==1) {
	  nmax <- min(c(nrow(train.set)-max(unlist(lapply(test.sets.list.k,nrow))),ncol(X)+1))
	  if (ncomp>=nmax) {
	    ncomp <- nmax-1
	    warning(paste0("'ncomp' re-set to ",ncomp))
	  }
	}
	model.k.temp <- if (model %in% c("PLS-DA/LDA","PLS-DA/QDA")) {
	  pls::plsr(train.set.Y~train.set.X,ncomp=ncomp,...)
	} else {
	  if (!is.null(Y.add)) {
	    pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,Y.add=train.set.Yadd,
		weights=train.set.weights,...)
	  } else {
	    pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,weights=train.set.weights,...)
	  }
	}
	models.list1[[i*k-(k-j)]] <- model.k.temp
	model.k <- if (model %in% c("PLS-DA/LDA","PPLS-DA/LDA")) {
	  if (!is.null(prior)) {
	    MASS::lda(as.matrix(model.k.temp$scores),train.set.Yfac,prior=prior,tol=1.0e-8)
	  } else {
	    MASS::lda(as.matrix(model.k.temp$scores),train.set.Yfac,tol=1.0e-8)
	  }
	} else {
	  if (!is.null(prior)) {
	    MASS::qda(as.matrix(model.k.temp$scores),train.set.Yfac,prior=prior,tol=1.0e-8)
	  } else {
	    MASS::qda(as.matrix(model.k.temp$scores),train.set.Yfac,tol=1.0e-8)
	  }
	}
	models.list2[[i*k-(k-j)]] <- model.k
	pred[as.numeric(rownames(test.set))] <- as.character(predict(model.k,predict(model.k.temp,test.set.X,
	  type="scores"),method=crit.DA)$class)
    }
    pred.correct <- pred==trueclass
    rate <- 1-sum(pred.correct)/length(pred.correct)
    NMC[i] <- rate
  }
  return(list(model=model,type="qual2",repet=repet,k=k,ncomp=ncomp,crit.DA=crit.DA,groups=groups,
    models1.list=models.list1,models2.list=models.list2,NMC=NMC))
}

