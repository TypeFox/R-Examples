# print summary objects
# Author: Yi Wang
# 05-Jan-2010

print.summary.manylm <- function (x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"),  dig.tst = max(1, min(5, digits - 1)),eps.Pvalue = .Machine$double.eps, ... )
{  
    resid 	<- x$residuals
    df  	<- x$df
    rdf 	<- df[2]
    n   	<- NROW(resid)
    show.est    <- x$show.est
    test	<- x$test
    symbolic.cor <- x$symbolic.cor

    n.bootsdone <- x$n.bootsdone
    if(all(n.bootsdone==n.bootsdone[1]))  n.bootsdone <- n.bootsdone[1] else
      n.bootsdone <- paste(n.bootsdone, collapse = ", ")
    ##################### BEGIN print Table of Coefficients ####################
    if (x$cor.type=="R")  
         corname <- "unconstrained correlation"
    else if (x$cor.type=="I")  
         corname <- "response assumed to be uncorrelated"
    else if (x$cor.type=="shrink") 
         corname <- paste("correlation matrix shrunk by parameter",round(x$shrink.param, digits = 2) )
    else if (x$cor.type=="blockdiag")
      corname <- paste("blockdiagonal correlation matrix with", x$shrink.param,"variables in each block")
    else if (x$cor.type=="augvar")
      corname <- paste("correlation matrix augmented with parameter",round(x$shrink.param, digits = 2))
    else corname <- ""
	
    if(x$resamp == "perm.resid")
       x$resamp <- "residual (without replacement)"
	
    # print significance tests
    if (length(x$aliased) == 0) 
	cat("\nNo Coefficients\n") 
    else {
        if (nsingular <- df[3] - df[1])
           cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep="")
        else if (!is.null(test)) {	
	   cat("\nTest statistics:\n")
	   coefs <- x$coefficients
	}
	est <- x$est
		
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
	    cn <- names(aliased)
            if (show.est) {
		est <- matrix(NA, length(aliased), NCOL(est) , dimnames = list(cn,               colnames(est)))
		est[!aliased, ] <- x$est 
	}    } 
	if (show.est) {
	    acs <- abs(coef.gv <- est[,1, drop = FALSE])
	    if (any(is.finite(acs)))  {
		# format the Generalised Variance Values
	       digmin <- 1 + floor(log10(range(acs[acs != 0], na.rm = TRUE)))
	       acs <- format(round(coef.gv, max(1, digits-digmin)), digits=digits)
	    }
	    # format the estimates values
	    estimates <- format(round(est[, 2:ncol(est)], digits = dig.tst),digits = digits)
	}
	
	if(!is.null(test)) {
  	    # format the test statistic
	    coefs <- x$coefficients
	    testvalue <- format(round(coefs[,1], digits = dig.tst), digits = digits)
            if (!is.logical(signif.stars) || is.na(signif.stars)) {
                warning("option \"show.signif.stars\" is invalid: assuming TRUE")
                signif.stars <- TRUE
            }
            zap.i <- 2
            pval <- coefs[,zap.i]
            ok <- !(is.na(pval))
            # formating the p values
            pval[ok]<- format.pval(pval[ok], digits = dig.tst, eps = eps.Pvalue)
            if(x$resamp=="none") pval[] <- ""

            signif.stars <- signif.stars && any(coefs[ok, 2] < 0.1)
	    tests <- cbind(testvalue, pval)
	    colnames(tests)<-colnames(coefs)[1:2]

            if (signif.stars) {
                Signif <- symnum(coefs[, 2], corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
               print.default(cbind(tests, format(Signif)), quote = FALSE, right = TRUE, na.print = "NA",...)
               cat("--- \nSignif. codes: ", attr(Signif, "legend"), "\n") 
	       if(x$p.uni == "none"){
                   if(x$resamp!="none")
		      cat("Arguments: with", n.bootsdone, "resampling iterations using",    x$resamp, "resampling and",corname, "\n")  
            }    } 
            else {
               if(x$resamp!="none")
	          print.default(tests, quote = FALSE, right = TRUE, na.print = "NA",...)
	       else 
                  print.default(tests[,-zap.i, drop=FALSE], quote = FALSE, right = TRUE, na.print = "NA",...)
               if (x$p.uni == "none" & x$resamp!="none")
                  cat("Arguments: with", n.bootsdone, "resampling iterations using",        x$resamp, "resampling and",corname, "\n")
	    }

  	   if(x$p.uni == "none" & x$resamp=="case" & sum(x$n.iter.sing)>0) {
		cat("\nNumber of iterations with adjusted tests (including skipped tests)      because of singularities in X due to the case resampling\n")
		print.default(x$n.iter.sing, quote = FALSE, right = TRUE, na.print = "", ...)	
		if(sum(x$nBoot-x$n.bootsdone)>0){
                   cat("\nNumber of iterations with skipped test statistic as the respective      variable to test became linear dependent during the case resampling step\n")
		   print.default(x$nBoot-x$n.bootsdone, quote = FALSE, right = TRUE,        na.print = "", ...) 
       }     }    } 
 
      if (show.est) {
	  cat("\nCoefficients:\n")
          print.default( cbind(acs, estimates), quote = FALSE, right = TRUE,
          na.print = "NA",...)
      } 
     if (x$p.uni != "none" & !is.null(test)) { 
	 if (!is.null(x$uni.p)){
            uni.p <- x$uni.p
            cols <- ncol(uni.p)   # = rank x
            coeffs.j <- matrix(NA, ncol=2*cols, nrow=nrow(x$uni.p))
            ok <- !(is.na(uni.p))
		uni.p[ok]<- format.pval(uni.p[ok], digits = dig.tst, eps = eps.Pvalue)
	    if(x$resamp=="none")  uni.p[] <- ""
            coeffs.j[,2*(1:cols)-1]  <- format(round(x$uni.test, digits = dig.tst), digits = digits)
            coeffs.j[,2*(1:cols)]    <- uni.p
            colna <- rep.int("", times=2*cols)
            colna[2*(1:cols)-1] <- colnames(x$uni.p)
            dimnames(coeffs.j) <- list(rownames(x$uni.p), colna )
            testname <- paste(x$test,"value")
            pname  <- paste("Pr(>",x$test,")", sep="") 
            coeffs.j[x$n.bootsdone==0,2*(1:cols)]  <-  NA 
            cona <- rep.int(c(testname, pname), times=cols) 
	    coeffs.j <- rbind(cona, coeffs.j) 
	    rownames(coeffs.j)[1] <- ""
	}
        cat(paste("\nUnivariate test statistic: \n" ))
        if(x$resamp!="none"){
            print.default(coeffs.j, quote =FALSE, right = TRUE, na.print = "NA",...)
            cat("\nArguments: with", n.bootsdone, "resampling iterations using",             x$resamp, "resampling and",corname, "\n")
        } 
        else 
            print.default(coeffs.j[,2*(1:cols)-1, drop=FALSE], quote =FALSE, right = TRUE, na.print = "NA",...)
	    if(x$resamp=="case" & sum(x$n.iter.sing)>0) {
		cat("\nNumber of iterations with adjusted tests (including skipped tests)      because of singularities in X due to the case resampling\n")
		print.default(x$n.iter.sing, quote = FALSE, right = TRUE, na.print = "", ...)	
		if (sum(x$nBoot-x$n.bootsdone)>0) {
		    cat("\nNumber of iterations with skipped test statistic as the respective      variable to test became linear dependent during the case resampling step\n")
		    print.default(x$nBoot-x$n.bootsdone, quote = FALSE, right = TRUE, na.print = "", ...) 
   }     }    }   } 
    ###################### END print Table of Coefficients #####################

    ##################### BEGIN print overall test statistics ##################
    if (!is.null(x$statistic)) {
       # ie if (is.null(test)) or if X only consists of the Intercept of 0
       if (!is.null(x$R2) & (length(x$r.squared)==1))
	  cat(paste("\n", x$R2, ":",sep=""), formatC(x$r.squared, digits = digits),"\n")
       else if (!is.null(x$R2)){
	  cat("\n", x$R2, ":\n",sep="" )
	  print.default(formatC(x$r.squared, digits = digits),quote = FALSE, right = TRUE, na.print = "NA",...)
	  cat("\n")
	}
       if (x$test=="LR")
          statname <- "\nLikelihood Ratio statistic: "
       else if (x$test=="F") statname <- "\nLawley-Hotelling trace statistic: "
       else statname <- "\nTest statistic: "

       if(x$resamp!="none"){
          cat(statname, paste(formatC(x$statistic[1], digits = digits),",",sep=""),         "p-value:", format.pval(x$statistic[2], digits = dig.tst, eps = eps.Pvalue),"\n")
          if (x$p.uni == "none") 
             cat("Arguments: with", n.bootsdone, "resampling iterations using", x$resamp, "resampling and",corname, "\n")    
       } 
       else cat(statname, paste(formatC(x$statistic[1], digits = digits),"\n",sep=""))
    
       if (x$p.uni != "none") {
          zap.ij <- 2
	  pvalj <- x$statistic.j[,zap.ij]
	  ok <- !(is.na(pvalj))
	  pvalj[ok]<- format.pval(pvalj[ok], digits = dig.tst, eps = eps.Pvalue)
	  if(x$resamp=="none") pvalj[] <- ""
	  x$statistic.j[,1] <- format(round(x$statistic.j[,1], digits = dig.tst), digits = digits)
	  x$statistic.j[,2] <- pvalj
          cat("\nUnivariate test statistic: \n")
          if(x$resamp!="none"){
             print.default(t(x$statistic.j), quote =FALSE, right = TRUE, na.print = "NA",...)
             cat("\nArguments: with", n.bootsdone, "resampling iterations using", x$resamp, "resampling and",corname, "\n")
          } else {
             uni.stat <- x$statistic.j[,-zap.ij, drop=FALSE]
             rownames(x$statistic.j)
             print.default(t(uni.stat),quote=FALSE,right=TRUE,na.print="NA",...)
    }   }   }
    ###################### END print overall test statistics ###################

    ######################## BEGIN print residual summary ######################
    if (x$show.residuals)
    {
      cat("\n\n",if (!is.null(x$w) && diff(range(x$w))) "Weighted ", "Residuals:\n", sep = "")
      if (rdf > 5| df[1]>4)  {
	  nam <- c("Min", "1Q", "Median", "3Q", "Max")
          rq <- if (length(dim(resid)) == 2)
                    structure(apply(t(resid), 1, quantile, na.rm = TRUE), dimnames = list(nam,     dimnames(resid)[[2]]))
	        else structure(quantile(resid, na.rm = TRUE), names = nam)
          print(rq, digits = digits, ...)
       }
      else if (rdf > 0) 
          print(resid, digits = digits, ...)
      else 
          cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")  

      cat("\nGeneralised Variance of the Residuals:", format(signif(x$genVar, digits)), "on", rdf, "degrees of freedom\n") 
    }
    ####################### END print residual summary #########################

    ####################### BEGIN print correlation ############################
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
           cat("\nCorrelation of Coefficients:\n")
           if (is.logical(symbolic.cor) && symbolic.cor)
              print(symnum(correl, abbr.colnames = NULL))
           else {
               correl <- format(round(correl, 2), nsmall = 2, digits = digits)
               correl[!lower.tri(correl)] <- ""
               print(correl[-1, -p, drop = FALSE], quote = FALSE)
    }   }   }
    ####################### ENDE print correlation #############################
    cat("\n")
    invisible(x)
}

# setMethod("print", "summary.manylm", print.summary.manylm )


