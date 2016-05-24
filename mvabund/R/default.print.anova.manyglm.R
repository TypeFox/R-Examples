# Pring anova objects
# Author: Yi Wang 
# 05-Jan-2010

# last modifed: David Warton, 06-May-2015

default.print.anova.manyglm <- function( x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"),  dig.tst = max(1, min(5, digits - 1)), eps.Pvalue = .Machine$double.eps, ...) 
{
   allargs <- match.call(expand.dots=FALSE)
   dots <- allargs$...
        
    anova   <- x
    x       <- anova$table
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
        warning("option \"show.signif.stars\" is invalid: assuming TRUE")
        signif.stars <- TRUE
    }
    test <- anova$test
    n.bootsdone <- anova$n.bootsdone
    if(all(n.bootsdone==n.bootsdone[1]))  n.bootsdone <- n.bootsdone[1] 
    else n.bootsdone <- paste(n.bootsdone, collapse = ", ")
 
    nModels <- nrow(x)

    modelnames <- dimnames(anova$table)[[1]]
    if(anova$resamp == "perm.resid")
       anova$resamp <- "residual permutation (without replacement)"
   if(anova$resamp == "pit.trap")
     anova$resamp <- "PIT-trap"
   if(anova$resamp == "monte.carlo")
     anova$resamp <- "parametric bootstrap"
   
    if (anova$cor.type=="R")  corname <- "unconstrained correlation response"
    else if (anova$cor.type=="I")  corname <- "uncorrelated response (for faster computation)"
    else if (anova$cor.type=="shrink") 
        corname <- "correlated response via ridge regularization"	  
    else corname <- ""

    if(is.null(anova$block))
      block.text=""
    else
     block.text= " block"
    ############## Anova Table for the simultaneous tests x ####################

    if (!is.null(heading <- attr(x, "heading"))) 
        cat(heading, sep = "\n")
    if (!is.null(title <- attr(x, "title")))    
        cat(title)   else cat("\n")

    nc <- dim(x)[2]
    if (is.null(cn <- colnames(x))) 
        stop("anova table must have colnames")

    # the p value is supposed to be in the last column of the table
    has.P <- substr(cn[nc], 1, 3) == "Pr("
    zap.i <- 1:(if (has.P) nc - 1 else nc)
    # Get columns with teststat.
    i <- which(substr(cn, 2, 7) == " value" | substr(cn, 3, 8) == " value")
    i <- c(i, which(!is.na(match(cn, c("Wald", "Score", "LR"))))) 
    # ie i is columns with teststat with name  ... " value"  or  one of "F", "LR"
    if (length(i)) zap.i <- zap.i[!(zap.i %in% i)]
    tst.i <- i
    if (length(i <- grep("Df$", cn)))     # df s not shown as zap.i
        zap.i <- zap.i[!(zap.i %in% i)]

    if(substr(anova$resamp,1,1)=="n") colnames(x)[nc[has.P]]  <- ""
    # "no p-values calculated as 'resample=none'
    printCoefmat(round(x, digits=dig.tst), digits = digits, signif.stars = signif.stars, has.Pvalue = has.P, P.values = has.P, cs.ind = NULL, zap.ind = zap.i, tst.ind = tst.i, na.print = "", ...)
    
    if(!is.null(test) & substr(anova$resamp,1,1)!="n"){
       if(substr(anova$p.uni,1,1)=="n") {
         if(dim(anova$uni.p)[2]>1)
         {   
           cat("Arguments:\n", "Test statistics calculated assuming", corname, 
               "\n P-value calculated using", n.bootsdone, "resampling iterations via",       paste(anova$resamp,block.text,sep=""), "resampling (to account for correlation in testing).\n")
         }
         if(dim(anova$uni.p)[2]==1)
         {   
           cat("Arguments: P-value calculated using", n.bootsdone, "resampling iterations via",       paste(anova$resamp,block.text,sep=""), "resampling (to account for correlation in testing).\n")
         }
#         cat("Arguments:\n", "Test statistics calculated assuming", corname, 
#              "\n P-value calculated using", n.bootsdone, "resampling iterations via",       paste(anova$resamp,block.text,sep=""), "resampling (to account for correlation in testing).\n")
                if(sum(anova$nBoot - anova$n.bootsdone)>1){
                    cat("\nNumber of iterations with skipped test statistic as the respective variable/variable-group to test became linear dependent during the case resampling step\n")
                    print.default(anova$nBoot - anova$n.bootsdone - 1, quote = FALSE, right = TRUE, na.print = "", ...) 
		}
        }
    }
  ########### END Anova Table for the simultaneous tests #################

  ###################### Anova Table for the univariate tests ################
    if (substr(anova$p.uni,1,1)!="n" & !is.null(test) ) {
    # no significance stars for the univariate table! 
        dimnam.ab <- colnames(anova$uni.p)
        col.dimnab <- rep.int("", times= 2*length(dimnam.ab))
        col.dimnab[2*(1:length(dimnam.ab))-1] <- dimnam.ab
        pmabund  <- ncol(anova$uni.p)
        testname <- anova$test
        pname    <- paste("Pr(>",testname,")", sep="")
        colna    <- c(rep.int(c(testname, pname), times=pmabund))

        uni.table <- matrix(NA, nrow(anova$uni.p), pmabund*2)

        uni.table[,2*(1:pmabund)-1]  <- round(anova$uni.test, digits=dig.tst)
        uni.table[,2*(1:pmabund)]    <- round(anova$uni.p, digits=digits)

        # rbind( colna, uni.table)
	dimnames(uni.table) <- list(c( rownames(anova$uni.p)), col.dimnab)
        if (!is.null(heading.uni <- attr(uni.table, "heading"))) 
            cat(heading.uni, sep = "\n")
        if (!is.null(title.uni <- "\nUnivariate Tests:\n"))
          cat(title.uni)         
         if (is.null(col.names <- colna))
            stop("uni.table must have attribute columnames")
#         col.names  <- substr( col.names , 1,8)     
        
#        if (!anova$one){
#             first.line <- uni.table[1,]
#             first.line[is.na(first.line)]<-""
#         }
 
        # If test = NULL, or rank=0, there is no test and no test statistics
#        i.uni <- which(substr(col.names, 2, 7) == " value" | substr(col.names, 3, 8) == " value")
        # columns with teststat
#        i.uni <- c(i.uni, which(!is.na(match(col.names, c("wald", "score", "LR"))))) 
        # which ... has only TRUE value/s, if there is a column with name
        # "F" or "LR"   --> should be changed to sth. more general
        
         zap.iuni <- which(substr(col.names, 1, 3) == "Pr(" ) # pvalues 
         pvalj <- uni.table[, zap.iuni, drop = FALSE]
         ok <- !(is.na(pvalj))
         pvalj[ok]<- format.pval(round(pvalj[ok], digits=dig.tst), digits = dig.tst, eps=eps.Pvalue)   
#         if(!anova$one) uni.table[1,] <- first.line
         uni.table[,zap.iuni]   <- pvalj # [ok]  

         uni.table <- rbind(col.names, uni.table)
         rownames(uni.table)[1] <- ""

         if(substr(anova$resamp,1,1)=="n") {
            print.default(uni.table[,-zap.iuni, drop=FALSE], quote = FALSE,
              right = TRUE, na.print = "", ...)
         } else print.default(uni.table, quote = FALSE, right = TRUE,
              na.print = "", ...)

        if( substr(anova$resamp,1,1)!="n"){
              cat("Arguments:\n", "Test statistics calculated assuming", corname, 
              "\nP-value calculated using", n.bootsdone, "resampling iterations via",       paste(anova$resamp,block.text,sep=""), "resampling (to account for correlation in testing.\n")
           if(sum(anova$nBoot - anova$n.bootsdone)>1){
              cat("\nNumber of iterations with skipped test statistic as the respective              variable/variable-group to test became linear dependent during the case resampling step\n")
              print.default(anova$nBoot - anova$n.bootsdone - 1, quote = FALSE, right = TRUE, na.print = "", ...)
            }
        }
    }

    ###################### END Anova Table for the univariate tests ############
    cat("\n")
    invisible(anova)
}


