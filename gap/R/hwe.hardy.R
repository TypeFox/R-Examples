hwe.hardy<-function(a, alleles=3, seed=3000, sample=c(1000, 1000, 5000)) {
#   require(genetics)
#    if (!missing(x)) {
#        if (!is.genotype(x)) {
#            stop("'x' must be of class 'genotype' or 'haplotype'")
#        } else {
#            # Get genotype counts
#            tab <- table(factor(allele(x, 1), levels=allele.names(x)),
#                         factor(allele(x, 2), levels=allele.names(x)))
#            a <- as.integer(t(tab)[lower.tri(t(tab), diag=T)])
#            a <- a[order(a)]
#            # Get number of alleles
#            alleles <- length(allele.names(x))
#        }        
#    }
    if (alleles<3) stop("number of alleles should be at least 3")
    p <- 1.0
    se <- 0.0
    swp <- rep(0,3)
    z <- .C("hwe_hardy", a=as.integer(a), alleles=as.integer(alleles),
            seed=as.integer(seed), gss=as.integer(sample),
            p=as.double(p), se=as.double(se), swp=as.double(swp),
            PACKAGE="gap")
    # Printout (partly taken from htest)
    z$method <- "Hardy-Weinberg equilibrium test using MCMC"
    z$data.name <- deparse(substitute(x))
    cat("\n")
    writeLines(strwrap(z$method, prefix = "\t"))
    cat("\n")
    cat("data: ", z$data.name, "\n")
    out <- character()
    fp <- format.pval(z$p, digits=4)
    out <- c(out, paste("p-value",
                  if(substr(fp,1,1) == "<") fp else paste("=",fp)))    
    out <- c(out, paste("p-value.se", "=", format(round(z$se, 4))))
    writeLines(strwrap(paste(out, collapse = ", ")))    
    cat("percentage of switches:\n")    
    cat(paste(" - partial", "=", format(round(z$swp[1], 2))), "\n")
    cat(paste(" - full", "=", format(round(z$swp[2], 2))), "\n")
    cat(paste(" - altogether", "=", format(round(z$swp[3], 2))), "\n\n")
    RVAL <- list(method=z$method, data.name=z$data.name, p.value=z$p, 
                 p.value.se=z$se, switches=z$swp)
    names(RVAL$switches) <- c("partial", "full", "altogether")
    return(invisible(RVAL))
}
