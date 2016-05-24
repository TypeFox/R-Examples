# pls : plsr,cppls
# MASS : lda,qda

print.MVA.cmv <- function(x,...) {
  cat("\n        Cross model validation (2CV)\n\n")
  cat(paste0("Model:"),x$model,"\n")
  cat(paste0("Inner loop: ",x$kinn,"-fold validation\n"))
  cat(paste0("Outer loop: ",x$kout,"-fold validation\n"))
  if (x$repet>1) {cat(paste0("Validation repeated ",x$repet," times\n"))}
  comp <- if (x$type!="qual2") {
    range(unlist(lapply(x$models.list,function(y) y$ncomp)))
  } else {
    range(unlist(lapply(x$models1.list,function(y) y$ncomp)))
  }
  cat(paste0(x$repet*x$kout," ",ifelse(x$type=="qual2","couples of ",""),"submodels generated (",
    comp[1]," to ",comp[2]," components)\n\n"))
  cat(paste0("Inner loop criterion:"),ifelse(x$crit.inn=="NMC","number of misclassifications",x$crit.inn))
  cat(ifelse(x$crit.inn=="Q2",paste0(" (",100*x$Q2diff," % threshold)\n"),"\n"))
  if (x$type=="qual2") {cat(paste0("Classification criterion (inner and outer loops):"),x$crit.DA,"\n")}
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

MVA.cmv <- function(X,Y,repet=10,kout=7,kinn=8,ncomp=8,model=c("PLSR","CPPLS","PLS-DA","PPLS-DA","PLS-DA/LDA",
  "PLS-DA/QDA","PPLS-DA/LDA","PPLS-DA/QDA"),crit.inn=c("RMSEP","Q2","NMC"),Q2diff=0.05,lower=0.5,upper=0.5,
  Y.add=NULL,weights=rep(1,nrow(X)),set.prior=FALSE,crit.DA=c("plug-in","predictive","debiased"),...) {
  model <- match.arg(model)
  crit.inn <- match.arg(crit.inn)
  crit.DA <- match.arg(crit.DA)
  type <- if (model %in% c("PLSR","CPPLS")) {"quant"} else {
    if (model %in% c("PLS-DA","PPLS-DA")) {"qual1"} else {"qual2"}
  }
  if (type=="quant" & is.factor(Y)) {
    if (model=="PLSR") {
	warning("'model' re-set to 'PLS-DA' and 'crit.inn' to 'NMC'")
    } else {
	warning("'model' re-set to 'PPLS-DA' and 'crit.inn' to 'NMC'")
    }
    model <- ifelse(model=="PLSR","PLS-DA","PPLS-DA")
    type <- "qual1"
    crit.inn <- "NMC"
  }
  if (type=="quant") {
    if (!crit.inn %in% c("RMSEP","Q2")) {stop("'crit.inn' must be 'RMSEP' or 'Q2'")}
  } else {
    if (!crit.inn %in% c("NMC")) {warning("'crit.inn' re-set to 'NMC'")}
    crit.inn <- "NMC"
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
  fun <- switch(type,quant=MVA.cmv.quant,qual1=MVA.cmv.qual1,qual2=MVA.cmv.qual2)
  res <- if (type=="quant") {fun(X,Y,repet,kout,kinn,ncomp,model,crit.inn,Q2diff,lower,upper,Y.add,
    weights,...)} else if (type=="qual1") {fun(X,Y,groups=levels(Yfac),repet,kout,kinn,ncomp,model,crit.inn,lower,
    upper,Y.add,weights,...)} else {fun(X,Y,Yfac,groups=levels(Yfac),repet,kout,kinn,ncomp,model,crit.inn,lower,upper,
    Y.add,weights,prior,crit.DA=crit.DA,...)}
  class(res) <- c("list","MVA.cmv")
  return(res)
}

MVA.cmv.quant <- function(X,Y,repet,kout,kinn,ncomp,model,crit.inn,Q2diff,lower,upper,Y.add,weights,...) {
  whole.set <- as.data.frame(cbind(weights,Y.add,Y,X))
  rownames(whole.set) <- 1:nrow(whole.set)
  col.Yadd <- if (!is.null(Y.add)) {2:(2+ncol(Y.add)-1)} else {NULL}
  col.Y <- if (!is.null(Y.add)) {(2+ncol(Y.add)):(2+ncol(Y.add)+ncol(Y)-1)} else {2:(2+ncol(Y)-1)}
  col.X <- if (!is.null(Y.add)) {(2+ncol(Y.add)+ncol(Y)):ncol(whole.set)} else {(2+ncol(Y)):ncol(whole.set)}
  test.sets.list.repet <- list()
  length(test.sets.list.repet) <- repet
  test.sets.list.repet <- lapply(test.sets.list.repet,function(x) {split(whole.set,sample(gl(kout,1,nrow(whole.set))))})
  models.list <- list()
  length(models.list) <- repet*kout
  names(models.list) <- paste(rep(1:repet,each=kout),rep(1:kout,repet),sep=":")
  N <- nrow(Y)
  TSS <- apply(Y,2,function(x) {sum((x-mean(x))^2)})
  RMSEP <- matrix(0,nrow=repet,ncol=ncol(Y),dimnames=list(1:repet,colnames(Y)))
  Q2 <- matrix(0,nrow=repet,ncol=ncol(Y),dimnames=list(1:repet,colnames(Y)))
  for (i in 1:repet) {
    test.sets.list.out <- test.sets.list.repet[[i]]
    pred.out <- matrix(0,nrow=nrow(whole.set),ncol=ncol(Y),dimnames=list(1:nrow(whole.set),colnames(Y)))
    for (j in 1:kout) {
	test.set <- test.sets.list.out[[j]]
	test.set.weights <- test.set$weights
	test.set.Yadd <- test.set[,col.Yadd]
	test.set.Y <- as.matrix(as.data.frame(test.set[,col.Y]))
	test.set.X <- as.matrix(as.data.frame(test.set[,col.X]))
	rest.set <- whole.set[-as.numeric(rownames(test.set)),]
	rest.set.weights <- rest.set$weights
	rest.set.Yadd <- rest.set[,col.Yadd]
	rest.set.Y <- as.matrix(as.data.frame(rest.set[,col.Y]))
	rest.set.X <- as.matrix(as.data.frame(rest.set[,col.X]))
	rownames(rest.set) <- 1:nrow(rest.set)
	val.sets.list <- split(rest.set,sample(gl(kinn,1,nrow(rest.set))))
	if (j==1) {
	  nmax <- min(c(nrow(rest.set)-max(unlist(lapply(val.sets.list,nrow))),ncol(X)+1))
	  if (ncomp>=nmax) {
	    ncomp <- nmax-1
	    warning(paste0("'ncomp' re-set to ",ncomp))
	  }
	}
	pred.inn <- list()
	length(pred.inn) <- ncomp
	for (n in 1:ncomp) {pred.inn[[n]] <-matrix(0,nrow=nrow(rest.set),ncol=ncol(Y))}
	for (k in 1:kinn) {
	  val.set <- val.sets.list[[k]]
	  val.set.X <- as.matrix(as.data.frame(val.set[,col.X]))
	  train.set <- rest.set[-as.numeric(rownames(val.set)),]
	  train.set.weights <- train.set$weights
	  train.set.Yadd <- train.set[,col.Yadd]
	  train.set.Y <- as.matrix(as.data.frame(train.set[,col.Y]))
	  train.set.X <- as.matrix(as.data.frame(train.set[,col.X]))
	  model.kinn <- if (model=="PLSR") {
	    pls::plsr(train.set.Y~train.set.X,ncomp=ncomp,...)
	  } else {
	    if (!is.null(Y.add)) {
		pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,Y.add=train.set.Yadd,
		  weights=train.set.weights,...)
	    } else {
		pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,weights=train.set.weights,...)
	    }
	  }
	  pred <- predict(model.kinn,newdata=val.set.X)
	  for (l in 1:ncomp) {pred.inn[[l]][as.numeric(rownames(val.set)),] <- pred[,,l]}
	}
	PRESS.inn <- do.call("rbind",lapply(pred.inn,function(x) {colSums((rest.set.Y-as.data.frame(x))^2)}))
	if (crit.inn=="RMSEP") {
	  N.inn <- nrow(rest.set)
	  RMSEP.inn <- sqrt(PRESS.inn/N.inn)
	  RMSEP.inn.sum <- rowSums(RMSEP.inn)
	  ncomp.kept <- which.min(RMSEP.inn.sum)[1]
	} else  if (crit.inn=="Q2") {
	  TSS.inn <- apply(rest.set.Y,2,function(x) {sum((x-mean(x))^2)})
	  Q2.inn <- 1-t(apply(PRESS.inn,1,function(x) {x/TSS.inn}))
	  if (ncol(Y)==1) {Q2.inn <- t(Q2.inn)}
	  Q2.inn.mean <- rowMeans(Q2.inn)
	  Q2d.inn <- numeric(length(Q2.inn.mean)-1)
	  for (q in 2:length(Q2.inn.mean)) {
	    Q2d.inn[q-1] <- abs((Q2.inn.mean[q]-Q2.inn.mean[q-1])/Q2.inn.mean[q-1])
	    Q2d.inn[q-1] <- Q2d.inn[q-1]*sign(Q2.inn.mean[q])
	  }
	  if (all(Q2d.inn<Q2diff)) {
	    ncomp.kept <- 1
	  } else if (all(Q2d.inn>=Q2diff)) {
	    ncomp.kept <- length(Q2.inn)
	  } else {
	    ncomp.kept <- which(Q2d.inn<Q2diff)[1]
	  }
	}
	model.kout <- if (model=="PLSR") {
	  pls::plsr(rest.set.Y~rest.set.X,ncomp=ncomp.kept,...)
	} else {
	  if (!is.null(Y.add)) {
	    pls::cppls(rest.set.Y~rest.set.X,ncomp=ncomp.kept,lower=lower,upper=upper,Y.add=rest.set.Yadd,
		weights=rest.set.weights,...)
	  } else {
	    pls::cppls(rest.set.Y~rest.set.X,ncomp=ncomp.kept,lower=lower,upper=upper,weights=rest.set.weights,...)
	  }
	}
	models.list[[i*kout-(kout-j)]] <- model.kout
	pred.out[as.numeric(rownames(test.set)),] <- predict(model.kout,newdata=test.set.X,ncomp=ncomp.kept)
    }
    PRESS <- colSums((Y-pred.out)^2)
    RMSEP[i,] <- sqrt(PRESS/N)
    Q2[i,] <- 1-PRESS/TSS
  }
  return(list(model=model,type="quant",repet=repet,kout=kout,kinn=kinn,crit.inn=crit.inn,Q2diff=Q2diff,
    models.list=models.list,RMSEP=RMSEP,Q2=Q2))
}

MVA.cmv.qual1 <- function(X,Y,groups,repet,kout,kinn,ncomp,model,crit.inn,lower,upper,Y.add,weights,...) {
  whole.set <- as.data.frame(cbind(weights,Y.add,Y,X))
  rownames(whole.set) <- 1:nrow(whole.set)
  trueclass <- apply(Y,1,function(x) {colnames(Y)[which(x==1)]})
  col.Yadd <- if (!is.null(Y.add)) {2:(2+ncol(Y.add)-1)} else {NULL}
  col.Y <- if (!is.null(Y.add)) {(2+ncol(Y.add)):(2+ncol(Y.add)+ncol(Y)-1)} else {2:(2+ncol(Y)-1)}
  col.X <- if (!is.null(Y.add)) {(2+ncol(Y.add)+ncol(Y)):ncol(whole.set)} else {(2+ncol(Y)):ncol(whole.set)}
  test.sets.list.repet <- list()
  length(test.sets.list.repet) <- repet
  test.sets.list.repet <- lapply(test.sets.list.repet,function(x) {splitf(whole.set,factor(trueclass),kout)})
  test.sets.repet.trueclass <- lapply(test.sets.list.repet,function(x) {lapply(x,function(y) {trueclass[as.numeric(rownames(y))]})})
  models.list <- list()
  length(models.list) <- repet*kout
  names(models.list) <- paste(rep(1:repet,each=kout),rep(1:kout,repet),sep=":")
  NMC <- numeric(repet)
  for (i in 1:repet) {
    test.sets.list.out <- test.sets.list.repet[[i]]
    pred.out <- character(nrow(Y))
    for (j in 1:kout) {
	test.set <- test.sets.list.out[[j]]
	test.set.weights <- test.set$weights
	test.set.Yadd <- test.set[,col.Yadd]
	test.set.Y <- as.matrix(as.data.frame(test.set[,col.Y]))
	test.set.X <- as.matrix(as.data.frame(test.set[,col.X]))
	rest.set <- whole.set[-as.numeric(rownames(test.set)),]
	rest.set.weights <- rest.set$weights
	rest.set.Yadd <- rest.set[,col.Yadd]
	rest.set.Y <- as.matrix(as.data.frame(rest.set[,col.Y]))
	rest.set.X <- as.matrix(as.data.frame(rest.set[,col.X]))
	rest.set.trueclass <- trueclass[-as.numeric(rownames(test.set))]
	rownames(rest.set) <- names(rest.set.trueclass) <- 1:nrow(rest.set)
	val.sets.list <- splitf(rest.set,factor(rest.set.trueclass),kinn)
	if (j==1) {
	  nmax <- min(c(nrow(rest.set)-max(unlist(lapply(val.sets.list,nrow))),ncol(X)+1))
	  if (ncomp>=nmax) {
	    ncomp <- nmax-1
	    warning(paste0("'ncomp' re-set to ",ncomp))
	  }
	}
	pred.inn <- matrix("",nrow=nrow(rest.set),ncol=ncomp,dimnames=list(rownames(rest.set),1:ncomp))
	for (k in 1:kinn) {
	  val.set <- val.sets.list[[k]]
	  val.set.X <- as.matrix(as.data.frame(val.set[,col.X]))
	  train.set <- rest.set[-as.numeric(rownames(val.set)),]
	  train.set.weights <- train.set$weights
	  train.set.Yadd <- train.set[,col.Yadd]
	  train.set.Y <- as.matrix(as.data.frame(train.set[,col.Y]))
	  train.set.X <- as.matrix(as.data.frame(train.set[,col.X]))
	  model.kinn <- if (model=="PLS-DA") {
	    pls::plsr(train.set.Y~train.set.X,ncomp=ncomp,...)
	  } else {
	    if (!is.null(Y.add)) {
		pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,Y.add=train.set.Yadd,
		  weights=train.set.weights,...)
	    } else {
		pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,weights=train.set.weights,...)
	    }
	  }
	  pred.inn.dummy <- predict(model.kinn,newdata=val.set.X)
	  pred.inn.lev <- apply(pred.inn.dummy,3,function(x) {apply(x,1,function(y) {colnames(Y)[which.max(y)]})})
	  pred.inn[as.numeric(rownames(val.set)),] <- pred.inn.lev
	}
	pred.inn.correct <- apply(pred.inn,2,function(x) {x==rest.set.trueclass})
	rate.inn <- 1-apply(pred.inn.correct,2,function(x) sum(x)/length(x))
	ncomp.kept <- which.min(rate.inn)[1]
	model.kout <- if (model=="PLS-DA") {
	  pls::plsr(rest.set.Y~rest.set.X,ncomp=ncomp.kept,...)
	} else {
	  if (!is.null(Y.add)) {
	    pls::cppls(rest.set.Y~rest.set.X,ncomp=ncomp.kept,lower=lower,upper=upper,Y.add=rest.set.Yadd,
		weights=rest.set.weights,...)
	  } else {
	    pls::cppls(rest.set.Y~rest.set.X,ncomp=ncomp.kept,lower=lower,upper=upper,weights=rest.set.weights,...)
	  }
	}
	models.list[[i*kout-(kout-j)]] <- model.kout
	pred.out.dummy <- predict(model.kout,newdata=test.set.X,ncomp=ncomp.kept)
	pred.out.lev <- apply(pred.out.dummy,1,function(x) {colnames(Y)[which.max(x)]})
	pred.out[as.numeric(rownames(test.set))] <- pred.out.lev
    }
    pred.out.correct <- pred.out==trueclass
    rate.out <- 1-sum(pred.out.correct)/length(pred.out.correct)
    NMC[i] <- rate.out
  }
  return(list(model=model,type="qual1",repet=repet,kout=kout,kinn=kinn,crit.inn=crit.inn,groups=groups,
    models.list=models.list,NMC=NMC))
}

MVA.cmv.qual2 <- function(X,Y,Yfac,groups,repet,kout,kinn,ncomp,model,crit.inn,lower,upper,Y.add,weights,prior,crit.DA,...) {
  whole.set <- as.data.frame(cbind(weights,Y.add,Y,X))
  rownames(whole.set) <- 1:nrow(whole.set)
  trueclass <- as.character(Yfac)
  col.Yadd <- if (!is.null(Y.add)) {2:(2+ncol(Y.add)-1)} else {NULL}
  col.Y <- if (!is.null(Y.add)) {(2+ncol(Y.add)):(2+ncol(Y.add)+ncol(Y)-1)} else {2:(2+ncol(Y)-1)}
  col.X <- if (!is.null(Y.add)) {(2+ncol(Y.add)+ncol(Y)):ncol(whole.set)} else {(2+ncol(Y)):ncol(whole.set)}
  test.sets.list.repet <- list()
  length(test.sets.list.repet) <- repet
  test.sets.list.repet <- lapply(test.sets.list.repet,function(x) {splitf(whole.set,Yfac,kout)})
  test.sets.repet.trueclass <- lapply(test.sets.list.repet,function(x) {lapply(x,function(y) {trueclass[as.numeric(rownames(y))]})})
  models.list1 <- models.list2 <- list()
  length(models.list1) <- length(models.list2) <- repet*kout
  names(models.list1) <- names(models.list2) <- paste(rep(1:repet,each=kout),rep(1:kout,repet),sep=":")
  NMC <- numeric(repet)
  for (i in 1:repet) {
    test.sets.list.out <- test.sets.list.repet[[i]]
    pred.out <- character(nrow(Y))
    for (j in 1:kout) {
	test.set <- test.sets.list.out[[j]]
	test.set.weights <- test.set$weights
	test.set.Yadd <- test.set[,col.Yadd]
	test.set.Y <- as.matrix(as.data.frame(test.set[,col.Y]))
	test.set.X <- as.matrix(as.data.frame(test.set[,col.X]))
	rest.set <- whole.set[-as.numeric(rownames(test.set)),]
	rest.set.weights <- rest.set$weights
	rest.set.Yadd <- rest.set[,col.Yadd]
	rest.set.Y <- as.matrix(as.data.frame(rest.set[,col.Y]))
	rest.set.X <- as.matrix(as.data.frame(rest.set[,col.X]))
	rest.set.trueclass <- trueclass[-as.numeric(rownames(test.set))]
	rest.set.Yfac <- Yfac[-as.numeric(rownames(test.set))]
	rownames(rest.set) <- names(rest.set.trueclass) <- 1:nrow(rest.set)
	val.sets.list <- splitf(rest.set,factor(rest.set.trueclass),kinn)
	if (j==1) {
	  nmax <- min(c(nrow(rest.set)-max(unlist(lapply(val.sets.list,nrow))),ncol(X)+1))
	  if (ncomp>=nmax) {
	    ncomp <- nmax-1
	    warning(paste0("'ncomp' re-set to ",ncomp))
	  }
	}
	pred.inn <- matrix("",nrow=nrow(rest.set),ncol=ncomp,dimnames=list(rownames(rest.set),1:ncomp))
	for (k in 1:kinn) {
	  val.set <- val.sets.list[[k]]
	  val.set.X <- as.matrix(as.data.frame(val.set[,col.X]))
	  train.set <- rest.set[-as.numeric(rownames(val.set)),]
	  train.set.weights <- train.set$weights
	  train.set.Yadd <- train.set[,col.Yadd]
	  train.set.Y <- as.matrix(as.data.frame(train.set[,col.Y]))
	  train.set.X <- as.matrix(as.data.frame(train.set[,col.X]))
	  train.set.Yfac <- rest.set.Yfac[-as.numeric(rownames(val.set))]
	  model.kinn.temp <- if (model %in% c("PLS-DA/LDA","PLS-DA/QDA")) {
	    pls::plsr(train.set.Y~train.set.X,ncomp=ncomp,...)
	  } else {
	    if (!is.null(Y.add)) {
		pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,Y.add=train.set.Yadd,
		  weights=train.set.weights,...)
	    } else {
		pls::cppls(train.set.Y~train.set.X,ncomp=ncomp,lower=lower,upper=upper,weights=train.set.weights,...)
	    }
	  }
	  models.kinn <- list()
	  length(models.kinn) <- ncomp
	  for (l in 1:ncomp) {
	    models.kinn[[l]] <- if (model %in% c("PLS-DA/LDA","PPLS-DA/LDA")) {
		if (!is.null(prior)) {
		  MASS::lda(as.matrix(model.kinn.temp$scores[,1:l]),train.set.Yfac,prior=prior,tol=1.0e-8)
		} else {
		  MASS::lda(as.matrix(model.kinn.temp$scores[,1:l]),train.set.Yfac,tol=1.0e-8)
		}
	    } else {
		if (!is.null(prior)) {
		  MASS::qda(as.matrix(model.kinn.temp$scores[,1:l]),train.set.Yfac,prior=prior,tol=1.0e-8)
		} else {
		  MASS::qda(as.matrix(model.kinn.temp$scores[,1:l]),train.set.Yfac,tol=1.0e-8)
		}
	    }
	    pred.inn[as.numeric(rownames(val.set)),l] <- as.character(predict(models.kinn[[l]],predict(model.kinn.temp,
		val.set.X,ncomp=1:l,type="scores"),method=crit.DA)$class)
	  }
	}
	pred.inn.correct <- apply(pred.inn,2,function(x) {x==rest.set.trueclass})
	rate.inn <- 1-apply(pred.inn.correct,2,function(x) sum(x)/length(x))
	ncomp.kept <- which.min(rate.inn)[1]
	model.kout.temp <- if (model %in% c("PLS-DA/LDA","PLS-DA/QDA")) {
	  pls::plsr(rest.set.Y~rest.set.X,ncomp=ncomp.kept,...)
	} else {
	  if (!is.null(Y.add)) {
	    pls::cppls(rest.set.Y~rest.set.X,ncomp=ncomp.kept,lower=lower,upper=upper,Y.add=rest.set.Yadd,
		weights=rest.set.weights,...)
	  } else {
	    pls::cppls(rest.set.Y~rest.set.X,ncomp=ncomp.kept,lower=lower,upper=upper,weights=rest.set.weights,...)
	  }
	}
	models.list1[[i*kout-(kout-j)]] <- model.kout.temp
	model.kout <- if (model %in% c("PLS-DA/LDA","PPLS-DA/LDA")) {
	  if (!is.null(prior)) {
	    MASS::lda(as.matrix(model.kout.temp$scores),rest.set.Yfac,prior=prior,tol=1.0e-8)
	  } else {
	    MASS::lda(as.matrix(model.kout.temp$scores),rest.set.Yfac,tol=1.0e-8)
	  }
	} else {
	  if (!is.null(prior)) {
	    MASS::qda(as.matrix(model.kout.temp$scores),rest.set.Yfac,prior=prior,tol=1.0e-8)
	  } else {
	    MASS::qda(as.matrix(model.kout.temp$scores),rest.set.Yfac,tol=1.0e-8)
	  }
	}
	models.list2[[i*kout-(kout-j)]] <- model.kout
	pred.out[as.numeric(rownames(test.set))] <- as.character(predict(model.kout,predict(model.kout.temp,test.set.X,
	  type="scores"),method=crit.DA)$class)
    }
    pred.out.correct <- pred.out==trueclass
    rate.out <- 1-sum(pred.out.correct)/length(pred.out.correct)
    NMC[i] <- rate.out
  }
  return(list(model=model,type="qual2",repet=repet,kout=kout,kinn=kinn,crit.inn=crit.inn,crit.DA=crit.DA,groups=groups,
    models1.list=models.list1,models2.list=models.list2,NMC=NMC))
}

