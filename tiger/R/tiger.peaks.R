tiger.peaks <- function(result, synthetic.errors){

    if(identical(result$synthetic.errors , synthetic.errors)){
       warning("Synthetic errors are equivalent to existing synthetic errors. Not performing any calculations")
       return(result)
    }

    result$synthetic.errors <- synthetic.errors
    result$has.peaks <- TRUE

    cat("Calculating measures for synthetic peaks\n")
    #Performance for synthetic errors
    peaks.length<-dim(synthetic.errors)[4]
    n.errors <- dim(synthetic.errors)[2]
    n.levels <- dim(synthetic.errors)[3]
    

    synth_diag <- array(dim= c(n.errors,n.levels,53))

    for(error in 1:n.errors){
       for(level in 1:n.levels){
           if(!is.na(synthetic.errors[1,error,level,1])){
               bla <- diagnostic_dawson(modelled=synthetic.errors[1,error,level,],measured=synthetic.errors[2,error,level,], use_qualV=TRUE )
               synth_diag[error,level,] <- drop(as.matrix(bla))
           }
       }
    }
    diag.names <- names(bla)

    #Remove unnecessary measures
    synth_diag.red<-synth_diag[,,-which(diag.names %in% c("AIC","BIC","span"))]  
    diag.names.red <- diag.names[!diag.names %in% c("AIC","BIC","span")]

    #Error quantile for synthetic peaks according to time series data
    if(result$multi.model){
         modelled <- drop.dimension(x=result$modelled)
         measured <- rep(result$measured, result$count.model)
    } else {
         modelled <- result$modelled
         measured <- result$measured
    }
    t.diff <-  measured- modelled
    diff.ecdf <- ecdf(t.diff)

    for(i in 1:n.errors) {
       rank.diff <- diff.ecdf(synthetic.errors[2,i,,] - synthetic.errors[1,i,,])
       dim(rank.diff) <- c(n.levels,peaks.length)
       synth_diag.red[,,which(diag.names.red=="EQ")][i+(1:n.levels-1)*n.errors] <-  apply(rank.diff, MARGIN=1, FUN=mean)
    }

    #relative diff nicht Inf (symmetrisches Verhalten von 0-Werten in Steigung)
    check.inf <- which(diag.names.red %in% c("rel_diff", "krel"))
    synth_diag.red[,,check.inf][is.infinite(synth_diag.red[,,check.inf])] <- 0 

    result$measures.synthetic.peaks <- synth_diag.red

    synth_diag.rer <- synth_diag.red 
    for(i in 1:length(diag.names.red)){
       stopifnot(diag.names.red[i]==names(result$measures[i]))
       synth_diag.rer[,,i] <-to.uniform(ref=result$measures[,i], val=synth_diag.red[,,i])
    }

    result$measures.uniform.synthetic.peaks <- synth_diag.rer

    if(result$use.som){
        result <- tiger.som.peaks(result)
    }
    result$error.names <- paste(c("peak size", "shift", "recession", "lag",
"size./integr", "width", "false peak", "undeteced peak", "shift w/o peak"), " (",1:9,")", sep="")

    return(result)


}
