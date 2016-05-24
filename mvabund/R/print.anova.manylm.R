# Pring anova objects
# Author: Yi Wang 
# 05-Jan-2010

print.anova.manylm <- function( x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"),  dig.tst = max(1, min(5, digits - 1)), eps.Pvalue = .Machine$double.eps, ...) 
{
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

    if(anova$resamp == "perm.resid")
      anova$resamp <- "residual (without replacement)"

    if (anova$cor.type=="R")  corname <- "unconstrained correlation"
    else if (anova$cor.type=="I")  corname <- "response assumed to be uncorrelated"
    else if (anova$cor.type=="shrink")
        corname <- paste("correlation matrix shrunk by parameter",
          round(anova$shrink.param, digits = dig.tst))
    else if (anova$cor.type=="blockdiag")
      corname <- paste("blockdiagonal correlation matrix with", anova$shrink.param,
       "variables in each block")
    else if (anova$cor.type=="augvar")
      corname <- paste("correlation matrix augmented with parameter",
          round(anova$shrink.param, digits = dig.tst))
    else corname <- ""

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
    i <- c(i, which(!is.na(match(cn, c("F", "LR"))))) 
    # ie i is columns with teststat with name  ... " value"  or  one of "F", "LR"
    if (length(i)) 
        zap.i <- zap.i[!(zap.i %in% i)]
    tst.i <- i
    if (length(i <- grep("Df$", cn)))     # df s not shown as zap.i
        zap.i <- zap.i[!(zap.i %in% i)]

    if(substr(anova$resamp,1,1)=="n") colnames(x)[nc[has.P]]  <- ""
    # "no p-values calculated as 'resample=none'

    printCoefmat(x, digits = digits, signif.stars = signif.stars, has.Pvalue = has.P, P.values = has.P, cs.ind = NULL, zap.ind = zap.i, tst.ind = tst.i, na.print = "", ...)
    
    if(!is.null(test) & substr(anova$resamp,1,1)!="n"){
        if(anova$p.uni=="none") {
            if(inherits(anova, "anova.manyglm") )
                cat("Arguments: with", n.bootsdone, "resampling iterations using",       anova$resamp, "resampling,", anova$teststat, "and",corname, "\n")
            else
                cat("Arguments: with", n.bootsdone, "resampling iterations using",        anova$resamp, "resampling and",corname, "\n") 
            if(anova$resamp=="case" & sum(anova$n.iter.sing)>0) {
                cat("\nNumber of iterations with adjusted tests (including skipped tests)      because of singularities in X due to the case resampling\n")
                print.default(anova$n.iter.sing, quote = FALSE, right = TRUE, na.print = "", ...)
                if(sum(anova$nBoot - anova$n.bootsdone)>0){
                    cat("\nNumber of iterations with skipped test statistic as the respective      variable/variable-group to test became linear dependent during the case resampling step\n")
                    print.default(anova$nBoot - anova$n.bootsdone, quote = FALSE, right = TRUE, na.print = "", ...) }
            }
        }
    }
    ############# END Anova Table for the simultaneous tests ###################
     
    ###################### RSS Table ###########################################
    if(anova$calc.rss) {
        RSStable <- anova$RSS$RSS
        if (!is.null(titleRSS <- attr(RSStable, "title")))  
            cat(titleRSS)   
        if(nrow(RSStable)==2){
            RSStable <- rbind(RSStable, anova$RSS$Diff)
            dimnames(RSStable)[[1]][3] <- "Diff Sum Sq"
            RSStable <- format(signif(RSStable, digits = dig.tst), digits = digits)
            print.default(RSStable, quote = FALSE, right = TRUE, na.print = "", ...)
        } else {
            RSStable <- format(signif(RSStable, digits = dig.tst), digits = digits)
            print.default(RSStable, quote = FALSE, right = TRUE, na.print = "", ...)
            cat("\n")
            Difftable <- anova$RSS$Diff
            if (!is.null(titleDiff <- attr(Difftable, "title")))    
                cat(titleDiff)  
            Difftable <- format(signif(Difftable, digits = dig.tst), digits = digits)
            print.default(Difftable, quote = FALSE, right = TRUE, na.print = "", ...)   
        }
    }
    ###################### END RSS Table #######################################
     
    ###################### Anova Table for the univariate tests ################
    if(anova$p.uni!="none" & !is.null(test) ) {
    # no significance stars for the univariate table! 

        dimnam.ab <- colnames(anova$uni.p)
        col.dimnab <- rep.int("", times= 2*length(dimnam.ab))
        col.dimnab[2*(1:length(dimnam.ab))-1] <- dimnam.ab
        pmabund  <- ncol(anova$uni.p)
        testname <- paste(anova$test,"value")
        pname    <- paste("Pr(>",anova$test,")", sep="")
        colna    <- c(rep.int(c(testname, pname), times=pmabund))

        uni.table <- matrix(NA, nrow(anova$uni.p), pmabund*2)

        uni.table[,2*(1:pmabund)-1]  <- round(anova$uni.test, digits=dig.tst)
        uni.table[,2*(1:pmabund)]    <- anova$uni.p

        # rbind( colna, uni.table)
        dimnames(uni.table) <- list(c( rownames(anova$uni.p)), col.dimnab) 
        if (!is.null(heading.uni <- attr(uni.table, "heading"))) 
            cat(heading.uni, sep = "\n")
        if (!is.null(title.uni <- "\nUnivariate Tests\nTest statistics:\n"))
          cat(title.uni)         
         if (is.null(col.names <- colna))
            stop("uni.table must have attribute columnames")
         col.names  <- substr( col.names , 1,8)     
        
        if (!anova$one){
             first.line <- uni.table[1,]
             first.line[is.na(first.line)]<-""
         }
 
        # If test = NULL, or rank=0, there is no test and no test statistics
        i.uni <- which(substr(col.names, 2, 7) == " value" | substr(col.names, 3, 8) == " value")
        # columns with teststat
        i.uni <- c(i.uni, which(!is.na(match(col.names, c("F", "LR"))))) 
        # which ... has only TRUE value/s, if there is a column with name
        # "F" or "LR"   --> should be changed to sth. more general
        
        zap.iuni <- which(substr(col.names, 1, 3) == "Pr(" ) 
         pvalj <- uni.table[,zap.iuni, drop = FALSE]
         ok <- !(is.na(pvalj))
         pvalj[ok]<- format.pval(pvalj[ok], digits = dig.tst, eps=eps.Pvalue)   
         if(!anova$one) uni.table[1,] <- first.line
         uni.table[,zap.iuni]   <- pvalj # [ok]  

         uni.table <- rbind(col.names, uni.table)
         rownames(uni.table)[1] <- ""

         if(substr(anova$resamp,1,1)=="n") {
            print.default(uni.table[,-zap.iuni, drop=FALSE], quote = FALSE,
              right = TRUE, na.print = "", ...)
         } else print.default(uni.table, quote = FALSE, right = TRUE,
              na.print = "", ...)


        if( substr(anova$resamp,1,1)!="n"){
           if(inherits(anova, "anova.manyglm") )
              cat("\nArguments: with", n.bootsdone, "resampling iterations using",              anova$resamp, "resampling,", anova$teststat, "and",corname, "\n")
           else 
              cat("\nArguments: with", n.bootsdone, "resampling iterations using", anova$resamp, "resampling and",corname, "\n")
           if(anova$resamp=="case" & sum(anova$n.iter.sing)>0) {
              cat("\nNumber of iterations with adjusted tests (including skipped tests)              because of singularities in X due to the case resampling\n")
              print.default(anova$n.iter.sing, quote = FALSE, right = TRUE, na.print = "", ...)
            }
           if(sum(anova$nBoot - anova$n.bootsdone)>0){
              cat("\nNumber of iterations with skipped test statistic as the respective              variable/variable-group to test became linear dependent during the case resampling step\n")
              print.default(anova$nBoot - anova$n.bootsdone, quote = FALSE, right = TRUE, na.print = "", ...)
            }
        }
    }

    ###################### END Anova Table for the univariate tests ############
    cat("\n")
    invisible(anova)
}


