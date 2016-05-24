###-----------------------------------------------------------------
# 7 Sep 2015
# 
#-----------------------------------------------------------------

# CBPS parses the formula object and passes the result to CBPS.fit
CBPS <- function(formula, data, na.action, ATT=1, iterations=1000, standardize=TRUE, method="over", twostep=TRUE, ...) {
		if (missing(data)) 
			data <- environment(formula)
		call <- match.call()
		family <- binomial()

		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data", "na.action"), names(mf), 0L)
		mf <- mf[c(1L, m)]
		mf$drop.unused.levels <- TRUE
		mf[[1L]] <- as.name("model.frame")
		
			mf <- eval(mf, parent.frame())
			mt <- attr(mf, "terms")
			Y <- model.response(mf, "any")
			if (length(dim(Y)) == 1L) {
					nm <- rownames(Y)
					dim(Y) <- NULL
					if (!is.null(nm)) 
						names(Y) <- nm
			}
			  
			X <- if (!is.empty.model(mt)) model.matrix(mt, mf)#[,-2]
			else matrix(, NROW(Y), 0L)
				
			X<-cbind(1,X[,apply(X,2,sd)>0])

			
			fit <- eval(call("CBPS.fit", X = X, treat = Y, ATT=ATT, 
						           intercept = attr(mt, "intercept") > 0L, method=method, iterations=iterations, 
                       standardize = standardize, twostep = twostep))	
				
			fit$na.action <- attr(mf, "na.action")
			xlevels <- .getXlevels(mt, mf)
      fit$data<-data
			fit$call <- call
			fit$formula <- formula
			fit$terms<-mt
		fit
	}

# CBPS.fit determines the proper routine (what kind of treatment) and calls the
# approporiate function.  It also pre- and post-processes the data
CBPS.fit<-function(treat, X, ATT, method, iterations, standardize, twostep, ...){
    # Special clause interprets T = 1 or 0 as a binary treatment, even if it is numeric
    if ((levels(factor(treat))[1] %in% c("FALSE","0",0)) & (levels(factor(treat))[2] %in% c("TRUE","1",1))
        & (length(levels(factor(treat))) == 2))
    {
      treat<-factor(treat)
    }

      # Declare some constants and orthogonalize Xdf.
      
      k=0
      if(method=="over") bal.only=FALSE
      if(method=="exact") bal.only=TRUE
      X.bal<-X
      
      names.X<-colnames(X)
      names.X[apply(X,2,sd)==0]<-"(Intercept)"
      
      X.orig<-X
      x.sd<-apply(as.matrix(X[,-1]),2,sd)
      Dx.inv<-diag(c(1,x.sd))
      diag(Dx.inv)<-1
      x.mean<-apply(as.matrix(X[,-1]),2,mean)
      X[,-1]<-apply(as.matrix(X[,-1]),2,FUN=function(x) (x-mean(x))/sd(x))
      k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
      k<-floor(k+.1)
      if (k < ncol(X)) stop("X is not full rank")
      svd1<-svd(X)
      X<-svd1$u
      XprimeX.inv<-ginv(t(X)%*%X)
      
      # Determine the number of treatments
      if (is.factor(treat)) {
        no.treats<-length(levels(treat))
        if (no.treats > 4) stop("Parametric CBPS not defined for more than 4 treatment values.  Consider using a continuous value.")
        if (no.treats < 2) stop("Treatment must take more than one value")
        
        if (no.treats == 2)
        {
          output<-CBPS.2Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, ATT, standardize = standardize, twostep = twostep)
        }
        
        if (no.treats == 3)
        {
          output<-CBPS.3Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep)
        }
        
        if (no.treats == 4)
        {
          output<-CBPS.4Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep)
        }
        
        # Reverse the svd, centering and scaling
        d.inv<- svd1$d
        d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
        d.inv[d.inv<= 1e-5]<-0
        beta.opt<-svd1$v%*%diag(d.inv)%*%coef(output)
        beta.opt[-1,]<-beta.opt[-1,]/x.sd
        
        beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
        output$coefficients<-beta.opt
        rownames(output$coefficients)<-names.X
        output$x<-X.orig
        
        # Calculate the variance
        variance<-output$var

        if (no.treats == 2){
          colnames(output$coefficients)<-c("Treated")
          output$var<-ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%variance%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)
          colnames(output$var)<-names.X
          rownames(output$var)<-colnames(output$var)
        }
        
        if (no.treats == 3){
          colnames(output$coefficients)<-levels(as.factor(treat))[c(2,3)]
          var.1.1<-variance[1:k,1:k]
          var.1.2<-variance[1:k,(k+1):(2*k)]
          var.2.1<-variance[(k+1):(2*k),1:k]
          var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
          trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2),cbind(trans.var.2.1,trans.var.2.2))
          colnames(output$var)<-c(paste0(levels(as.factor(treat))[2],": ", names.X),paste0(levels(as.factor(treat))[3], ": ", names.X))
          rownames(output$var)<-colnames(output$var)
        }
        
        if (no.treats == 4)
        {
          colnames(output$coefficients)<-levels(as.factor(treat))[c(2,3,4)]
          var.1.1<-variance[1:k,1:k]
          var.1.2<-variance[1:k,(k+1):(2*k)]
          var.1.3<-variance[1:k,(2*k+1):(3*k)]
          var.2.1<-variance[(k+1):(2*k),1:k]
          var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
          var.2.3<-variance[(k+1):(2*k),(2*k+1):(3*k)]
          var.3.1<-variance[(2*k+1):(3*k),1:k]
          var.3.2<-variance[(2*k+1):(3*k),(k+1):(2*k)]
          var.3.3<-variance[(2*k+1):(3*k),(2*k+1):(3*k)]
          trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.1.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.2.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.3.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.3.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          trans.var.3.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
          output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2,trans.var.1.3),cbind(trans.var.2.1,trans.var.2.2,trans.var.2.3),cbind(trans.var.3.1,trans.var.3.2,trans.var.3.3))
          colnames(output$var)<-c(paste0(levels(as.factor(treat))[2],": ", names.X),paste0(levels(as.factor(treat))[3], ": ", names.X),paste0(levels(as.factor(treat))[4], ": ", names.X))
          rownames(output$var)<-colnames(output$var)
        }
      } else if (is.numeric(treat)) {
        # Warn if it seems like the user meant to input a categorical treatment
        if (length(unique(treat)) <= 4) warning("Treatment vector is numeric.  Interpreting as a continuous treatment.  To solve for a binary or multi-valued treatment, make treat a factor.")
        output<-CBPS.Continuous(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep)
        
        # Reverse svd, centering, and scaling
        d.inv<- svd1$d
        d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
        d.inv[d.inv<= 1e-5]<-0
        beta.opt<-svd1$v%*%diag(d.inv)%*%coef(output)
        beta.opt[-1,]<-beta.opt[-1,]/x.sd
        
        beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
        output$coefficients<-as.matrix(beta.opt)
        rownames(output$coefficients)<-c(names.X)
        output$x<-X.orig
        
        # Calculate variance
        var.1<-output$var
        output$var<-ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)
        rownames(output$var)<-names.X
        colnames(output$var)<-rownames(output$var)
      } else {
        stop("Treatment must be either a factor or numeric")
      }

	  output
	}

  # Print coefficients and model fit statistics
	print.CBPS <- function(x, digits = max(3, getOption("digits") - 3), ...) {
		cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		if (length(coef(x))) {
			cat("Coefficients:\n")
			print.default(format(x$coefficients, digits = digits), 
				print.gap = 2, quote = FALSE)
		}
		else cat("No coefficients\n\n")
    if (max(class(x) == "CBPScontinuous"))
      cat("\nSigma-Squared: ",x$sigmasq)
		if (nzchar(mess <- naprint(x$na.action))) 
			cat("  (", mess, ")\n", sep = "")
		cat("Residual Deviance:\t", format(signif(x$deviance, 
			digits)), "\n")
		cat("J-Statistic:\t		", format(signif(x$J)),"\n")
		cat("Log-Likelihood:\t ",-0.5*x$deviance, "\n")
		invisible(x)
	}

  # Expands on print by including uncertainty for coefficient estimates
	summary.CBPS<-function(object, ...){
	  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)	  
	  x<-NULL
	  names.X<-as.vector(names(object$coefficients))
	  sd.coef <- diag(object$var)^.5
	  coef.table<-(cbind(as.vector(object$coefficients),as.vector(sd.coef),as.vector(object$coefficients/sd.coef),as.vector(2-2*pnorm(abs(object$coefficients/sd.coef)))))
	  colnames(coef.table)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
	  if (ncol(coef(object)) == 1)
	  {
		rownames(coef.table)<-rownames(object$coefficients)#names.X
	  }
	  if (ncol(coef(object)) > 1)
	  {
		  rnames<-array()
		  for (i in 1:ncol(coef(object)))
		  {
			rnames[((i-1)*nrow(coef(object))+1):(i*nrow(coef(object)))]<-paste0(levels(as.factor(object$y))[i],": ",rownames(coef(object)))
		  }
		  rownames(coef.table)<-rnames
	  }
	  
	  pval <- coef.table[,4]
	  symp <- symnum(pval, corr=FALSE,
					 cutpoints = c(0,  .001,.01,.05, .1, 1),
					 symbols = c("***","**","*","."," "))
	  coef.print<-cbind(signif(coef.table,3),as.vector(symp))
	  coef.print[coef.print=="0"]<-"0.000"
		
	  cat("\nCall:	\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
		  "\n", sep = "")
	  
		
	  cat("\nCoefficients:\n")

	  print(noquote(coef.print))
	  cat("---\n")
	  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
		#cat("\n	Null J:	 ",object$J)
    if(max(class(object)=="CBPScontinuous")){
      cat("\nSigma-Squared: ",object$sigmasq)
    }
	  cat("\nJ - statistic:	 ",object$J)
	  cat("\nLog-Likelihood: ",-0.5*object$deviance, "\n")
		
	  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
	  invisible(out)
	}
  

	vcov.CBPS<-function(object,...){
		return(object$var)
	}

  # Plot binary and multi-valued CBPS.  Plots the standardized difference in means for each contrast
  # before and after weighting.  Defined for an arbitrary number of discrete treatments.
	plot.CBPS<-function(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...){ 
		bal.x<-balance(x)
    if(is.null(covars))
    {
      covars<-1:nrow(bal.x[["balanced"]])
    }
    
		no.treats<-length(levels(as.factor(x$y)))	
		balanced.std.mean<-bal.x[["balanced"]][covars,]
		original.std.mean<-bal.x[["original"]][covars,]
		no.contrasts<-ifelse(no.treats == 2, 1, ifelse(no.treats == 3, 3, 6))
        
		abs.mean.ori.contrasts<-matrix(rep(0,no.contrasts*length(covars)),length(covars),no.contrasts)
		abs.mean.bal.contrasts<-matrix(rep(0,no.contrasts*length(covars)),length(covars),no.contrasts)
		contrast.names<-array()
		true.contrast.names<-array()
		contrasts<-c()
		covarlist<-c()
		ctr<-1
		for (i in 1:(no.treats-1))
		{
			for (j in (i+1):no.treats)
			{
				abs.mean.ori.contrasts[,ctr]<-abs(original.std.mean[covars,i+no.treats]-original.std.mean[covars,j+no.treats])
				abs.mean.bal.contrasts[,ctr]<-abs(balanced.std.mean[covars,i+no.treats]-balanced.std.mean[covars,j+no.treats])
				contrast.names[ctr]<-paste0(i,":",j)
				true.contrast.names[ctr]<-paste0(levels(as.factor(x$y))[i],":",levels(as.factor(x$y))[j])

        contrasts<-c(contrasts, rep(true.contrast.names[ctr],length(covars)))
        covarlist<-c(covarlist, rownames(balanced.std.mean))
        ctr<-ctr+1
			}
		} 

    max.abs.contrast<-max(max(abs.mean.ori.contrasts),max(abs.mean.bal.contrasts))
		m <- matrix(c(1,1,1,2,2,2,3,3,3),nrow = 3,ncol = 3,byrow = TRUE)
		layout(mat = m,heights = c(0.4,0.4,0.3))
    
		par(mfrow=c(2,1))
    if (!boxplot){
  		plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n', ...)
  		axis(side=2, at=seq(1,no.contrasts),contrast.names)
  		mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
  		mtext("Contrasts",side=2,line=2)
  		mtext("Before Weighting",side=3,line=0.5,font=2)
  		for (i in 1:no.contrasts)
  		{
  			for (j in 1:length(covars))
  			{
  				points(abs.mean.ori.contrasts[covars[j],i],i, ...)
  			}
  		}
  		plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts), xlab="", ylab="", main="", yaxt='n', ...)
  		axis(side=2, at=seq(1,no.contrasts),contrast.names)
  		mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
  		mtext("Contrasts",side=2,line=2)
  		mtext("After Weighting",side=3,line=0.5,font=2)
      
  		for (i in 1:no.contrasts)
  		{
  			for (j in 1:length(covars))
  			{
  				points(abs.mean.bal.contrasts[covars[j],i],i, ...)
  			}
  		}
    }
    else{
      boxplot(abs.mean.ori.contrasts, horizontal = TRUE, ylim = c(0,max.abs.contrast), xlim=c(1-0.5,no.contrasts+0.5), 
              xlab="", ylab="", main="", yaxt='n', ...)
      axis(side=2, at=seq(1,no.contrasts),contrast.names)
      mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
      mtext("Contrasts",side=2,line=2)
      mtext("Before Weighting",side=3,line=0.5,font=2)
      
      boxplot(abs.mean.bal.contrasts, horizontal = TRUE, ylim = c(0,max.abs.contrast), xlim=c(1-0.5,no.contrasts+0.5), 
              xlab="", ylab="", main="", yaxt='n', ...)
      axis(side=2, at=seq(1,no.contrasts),contrast.names)
      mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
      mtext("Contrasts",side=2,line=2)
      mtext("After Weighting",side=3,line=0.5,font=2)
    }
    
    par(mfrow=c(1,1))

		if(is.null(rownames(balanced.std.mean))) rownames(balanced.std.mean)<-paste0("X",covars)
    if(!silent) return(data.frame("contrast" = contrasts, "covariate" = covarlist, 
                                  "balanced"=abs.mean.bal.contrasts,  
                                  "original"=abs.mean.ori.contrasts))
	}
	
  # Plot the pre-and-post weighting correlations between X and T
	plot.CBPSContinuous<-function(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...){     
		bal.x<-balance(x)
    if (is.null(covars))
    {
      covars<-1:nrow(bal.x[["balanced"]])
    }
		balanced.abs.cor<-abs(bal.x[["balanced"]][covars])
		original.abs.cor<-abs(bal.x[["unweighted"]][covars])

		max.abs.cor<-max(max(original.abs.cor),max(balanced.abs.cor))
    
    if(!boxplot){
      plot(1, type="n", xlim=c(0,max.abs.cor), ylim=c(1.5,3.5), xlab = "Absolute Pearson Correlation", ylab = "", yaxt = "n", ...)
      axis(side=2, at=seq(2,3),c("CBPS Weighted", "Unweighted"))
      points(x=original.abs.cor, y=rep(3, length(covars)), pch=19)
      points(x=balanced.abs.cor, y=rep(2, length(covars)), pch=19)
    }
    else{
      boxplot(balanced.abs.cor, original.abs.cor, horizontal = TRUE, yaxt = 'n', xlab = "Absolute Pearson Correlation", ...)
      axis(side=2, at=c(1,2),c("CBPS Weighted", "Unweighted"))
    }
  
    if(!silent) return(data.frame("covariate"=rownames(bal.x[["balanced"]]),"balanced"=balanced.abs.cor,
                                  "original"=original.abs.cor))
	}

	balance<-function(object, ...)
	{
		UseMethod("balance")
	}

  # Calculates the pre- and post-weighting difference in standardized means for covariate within each contrast
	balance.CBPS<-function(object, ...){
		treats<-as.factor(object$y)
		treat.names<-levels(treats)
		X<-object$x
		bal<-matrix(rep(0,(ncol(X)-1)*2*length(treat.names)),ncol(X)-1,2*length(treat.names))
		baseline<-matrix(rep(0,(ncol(X)-1)*2*length(treat.names)),ncol(X)-1,2*length(treat.names))
		w<-object$weights
		cnames<-array()
		
		jinit<-ifelse(class(object)[1] == "npCBPS", 1, 2)
		for (i in 1:length(treat.names))
		{
			for (j in jinit:ncol(X))
			{
				bal[j-1,i]<-sum((treats==treat.names[i])*X[,j]*w)/sum(w*(treats==treat.names[i]))
				bal[j-1,i+length(treat.names)]<-bal[j-1,i]/sd(X[,j])
				baseline[j-1,i]<-mean(X[which(treats==treat.names[i]),j])
				baseline[j-1,i+length(treat.names)]<-baseline[j-1,i]/sd(X[,j])
			}
			cnames[i]<-paste0(treat.names[i],".mean")
			cnames[length(treat.names)+i]<-paste0(treat.names[i],".std.mean")
		}
		colnames(bal)<-cnames
		rownames(bal)<-colnames(X)[-1]
		colnames(baseline)<-cnames
		rownames(baseline)<-colnames(X)[-1]
		out<-list(balanced=bal,original=baseline)
		out
	}
	
  # Calculates the pre- and post-weighting correlations between each covariate and the T
	balance.CBPSContinuous<-function(object, ...){
		treat<-object$y
		X<-object$x
		w<-object$weights
    n<-length(w)
		cnames<-array()
		
    if ("npCBPS" %in% class(object)){
      jinit<-1
      bal<-matrix(rep(0,ncol(X)),ncol(X),1)
      baseline<-matrix(rep(0,ncol(X)),ncol(X),1)
      for (j in 1:ncol(X))
      {
        bal[j,1]<-(mean(w*X[,j]*treat) - mean(w*X[,j])*mean(w*treat)*n/sum(w))/(sqrt(mean(w*X[,j]^2) - mean(w*X[,j])^2*n/sum(w))*sqrt(mean(w*treat^2) - mean(w*treat)^2*n/sum(w)))
        baseline[j,1]<-cor(treat, X[,j], method = "pearson")
      }
      rownames(bal)<-colnames(X)
      rownames(baseline)<-colnames(X)   
    }
    else{
      bal<-matrix(rep(0,(ncol(X)-1)),ncol(X)-1,1)
      baseline<-matrix(rep(0,(ncol(X)-1)),ncol(X)-1,1)
      for (j in 2:ncol(X))
      {
        bal[j-1,1]<-(mean(w*X[,j]*treat) - mean(w*X[,j])*mean(w*treat)*n/sum(w))/(sqrt(mean(w*X[,j]^2) - mean(w*X[,j])^2*n/sum(w))*sqrt(mean(w*treat^2) - mean(w*treat)^2*n/sum(w)))
        baseline[j-1,1]<-cor(treat, X[,j], method = "pearson")
      }
      rownames(bal)<-colnames(X)[-1]
      rownames(baseline)<-colnames(X)[-1]      
    }
		
		colnames(bal)<-"Pearson Correlation"
		colnames(baseline)<-"Pearson Correlation"
		out<-list(balanced=bal,unweighted=baseline)
		out
	}	
	#########################################
