#############
## process ##
#############
process.bmmix <- function(x, what="post", burnin=1e4,
                         ggplot=TRUE){
    ## CHECKS ##
    if(!any(inherits(x, "bmmix"))) warning("x is not a bmmix object")
    x <- x[x$step>burnin,]
    toKeep <- grep(what, names(x))
    if(length(toKeep)==0) stop(paste(what, "not found in x"))


    ## EXTRACT RELEVANT COLUMNS ##
    out <- list(data.frame(step=x$step, x[,toKeep,drop=FALSE]))
    names(out) <- what

    ## GENERATE PLOTS ##
    if(ggplot){
       ## if(!require(ggplot2) | !require(reshape2)) stop("ggplot2 and reshape2 are needed.")

       ## MORE THAN ONE VARIABLE TO PLOT
       if(ncol(out[[1]])>2){
           ## reshape data if needed
           fig.dat <- melt(out[[1]], id=1)
           names(fig.dat) <- c("step", "group", what)
           fig.dat$group <- factor(sub(paste(what,".", sep=""), "", fig.dat$group))

           ## base
           base <- ggplot(data=fig.dat, aes_string(colour="group", fill="group"))

           ## trace
           out$traces <- base + geom_line(aes_string(x="step", y=what)) + labs(title=paste(title="trace of", what))

           ## densities
           out$densities <- base + aes_string(x=what) + geom_density(alpha=.3) +
               geom_point(aes_string(x=what, y=0), shape="|", alpha=.3) +
                   labs(title=paste(title="Density of", what))

           ## histograms
           out$histograms <- base + geom_histogram(aes_string(x=what), alpha=.3) +
               geom_point(aes_string(x=what, y=0), shape="|", alpha=.3) +
               labs(title=paste(title="Density of", what))

           ## violinplots
           out$violins <- base + geom_violin(aes_string(x="group", y=what), alpha=.3) +
               labs(title=paste(title="Violinplot of", what))
       } else {
           ## JUST ON VARIABLE TO PLOT
           fig.dat <- out[[1]]

           ## base
           base <- ggplot(data=fig.dat)

           ## trace
           out$traces <- base + geom_line(aes_string(x="step", y=what)) +
               labs(title=paste(title="trace of", what))

           ## densities
           out$densities <- base + geom_density(aes_string(x=what), alpha=.3) +
               geom_point(aes_string(x=what, y=0), shape="|", alpha=.3) +
                   labs(title=paste(title="Density of", what))

           ## histograms
           out$histograms <- base + geom_histogram(aes_string(x=what), alpha=.3) +
               geom_point(aes_string(x=what, y=0), shape="|", alpha=.3) +
                   labs(title=paste(title="Density of", what))
       }
   }

    ## RETURN ##
    return(out)
} # end process.bmmix
