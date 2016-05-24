# function to trim weights
ebalance.trim <-
  function(
           ebalanceobj,
           max.weight=NULL,
           min.weight=0,
           max.trim.iterations = 200,
           max.weight.increment=.92,
           min.weight.increment= 1.08,
           print.level=0
           )
    {
    if( class(ebalanceobj)!= "ebalance" ){
     stop("ebalanceobj must be an ebalance object from a call to ebalance()")
    }
    minimization <- FALSE
    if(is.null(max.weight)){
     max.weight <- max(ebalanceobj$w/mean(ebalanceobj$w))
     minimization <- TRUE
    }
    if(length(max.weight) != 1){
      stop("length(max.weight) != 1")
    }
    if(length(min.weight) != 1){
      stop("length(min.weight) != 1")
    }

  # starting setup for trimming
   w.trimming <- 1
   coefs <- ebalanceobj$coefs

 ### Trimming to user supplied max weight or trim once starting from max weight obtained from distribution)
   for(iter.trim in 1:max.trim.iterations) {
    if( minimization == FALSE ) {
    cat("Trim iteration",format(iter.trim,digits=3),"\n")
    }
      eb.out <- eb(tr.total=ebalanceobj$target.margins,
               co.x=ebalanceobj$co.xdata,
               coefs=coefs,
               base.weight=ebalanceobj$w*w.trimming,
               max.iterations=ebalanceobj$max.iterations,
               constraint.tolerance=ebalanceobj$constraint.tolerance,
               print.level=print.level
               )
     weights.ratio = eb.out$Weights.ebal/mean(eb.out$Weights.ebal)
     coefs <- eb.out$coefs
     if(max(weights.ratio) <= max.weight && min(weights.ratio) >= min.weight) {
      if( minimization == FALSE ) { cat("Converged within tolerance \n") }
     break
     }
      w.trimming <- w.trimming*ifelse(weights.ratio>max.weight,w.trimming*((max.weight*max.weight.increment)/weights.ratio),1)
     if(min.weight) w.trimming <- w.trimming*
       ifelse(weights.ratio<min.weight,w.trimming*((min.weight*min.weight.increment)/weights.ratio),1)
    }

 ### automated trimming to minimize max weight
   if( minimization == TRUE ) {
   cat("Automated trimmig of max weight ratio \n")
       for(iter.max.weight in 1:max.trim.iterations) {
        max.weight.old     <- max.weight
        weights.ratio.old  <- max(eb.out$Weights.ebal/mean(eb.out$Weights.ebal))
        eb.out.old <- eb.out # store old weights
         if(print.level >= 0) {
           cat("Trim iteration",format(iter.max.weight,digits=3),"Max Weight Ratio:",format(weights.ratio.old,digits=4),"\n")
           }
        suppressWarnings(IsError <- try(
            for(iter.trim in 1:max.trim.iterations) {
              eb.out <- eb(tr.total=ebalanceobj$target.margins,
               co.x=ebalanceobj$co.xdata,
               coefs=coefs,
               base.weight=ebalanceobj$w*w.trimming,
               max.iterations=ebalanceobj$max.iterations,
               constraint.tolerance=ebalanceobj$constraint.tolerance,
               print.level=0
               )
             weights.ratio <- eb.out$Weights.ebal/mean(eb.out$Weights.ebal)
             coefs <- eb.out$coefs
            if(max(weights.ratio) <= max.weight && min(weights.ratio)>=min.weight) break
             w.trimming <- w.trimming*ifelse(weights.ratio>max.weight,w.trimming*((max.weight*max.weight.increment)/weights.ratio),1)
            if(min.weight) w.trimming = w.trimming*
            ifelse(weights.ratio<min.weight,w.trimming*((min.weight*min.weight.increment)/weights.ratio),1)
           }
    ,silent=TRUE)) # end try call
    
    if(class(IsError)=="try-error") {
       if(print.level >= 2) { cat("no further decrease in max weight ratio \n") }
     break
     }
    if(weights.ratio.old < max(weights.ratio)) {
     if(print.level >= 2) { cat("no further decrease in max weight ratio \n") }
     break
    }
      # exit loop if alog doesn converge
     max.weight <- max.weight*max.weight.increment   # otherwise lower max weight and try again
    # cat("Moving into next minimisation loop with old max weight ratio:",max(weights.ratio),"and new attempted max weight ratio",max.weight,"\n")

   } # end max weight loop
  # play back the results
     eb.out <- eb.out.old
     if(eb.out$converged == TRUE) {
        cat("Converged within tolerance \n")
     }
  } # end automated if

z <- list(
          target.margins = ebalanceobj$tr.total,
          co.xdata = ebalanceobj$co.xdata,
          w=eb.out$Weights.ebal,
          coefs=eb.out$coefs,
          maxdiff=eb.out$maxdiff,
          norm.constant = ebalanceobj$norm.constant,
          constraint.tolerance=ebalanceobj$constraint.tolerance,
          max.iterations=ebalanceobj$max.iterations,
          base.weight=ebalanceobj$base.weight,
          converged =  eb.out$converged
    )

class(z) <- "ebalance.trim"
return(z)
    
   }
    

