synth.tab <-
function(
           synth.res    = NA,
           dataprep.res = NA,
           round.digit  = 3
           )

  {
    
    if(sum(is.na(synth.res))>0)
     {
      stop(
            "\n\n#####################################################",
            "\n No synth.res provided\n\n)"
            )
     }
     
         if(sum(is.na(dataprep.res))>0)
     {
      stop(
            "\n\n#####################################################",
            "\n No dataprep.res provided\n\n)"
            )
     }
     
     
    tab.v           <- matrix(round(synth.res$solution.v,round.digit),
                              nrow=length(synth.res$solution.v))
    rownames(tab.v) <- rownames(dataprep.res$X1)
    colnames(tab.v) <- c("v.weights")
   
    treat.no <-  dataprep.res$tag$treatment.identifier
    nmat     <-  dataprep.res$names.and.numbers[-which(
                                                       dataprep.res$names.and.numbers[,2]==treat.no)
                                               ,]
 
 
    tab.w           <- data.frame(
                                  round(synth.res$solution.w,round.digit),
                                  nmat
                                  )
    colnames(tab.w) <- c(
                         "w.weights",
                         "unit.names",
                         "unit.numbers"
                         )
  
    tab.loss           <- cbind(
                                synth.res$loss.w,
                                synth.res$loss.v
                                )
    rownames(tab.loss) <- NULL
    colnames(tab.loss) <- c(
                            "Loss W",
                            "Loss V"
                            )
    
    tab.pred <-     round(
                          cbind(
                                dataprep.res$X1,
                                as.matrix(dataprep.res$X0)%*%as.matrix(synth.res$solution.w),
                                apply(as.matrix(dataprep.res$X0),1,mean)
                                )
                          ,round.digit)
                          
    colnames(tab.pred) <- c(
                            "Treated",
                            "Synthetic",
                            "Sample Mean"
                            )
                        
    output <- list(
                   tab.pred = tab.pred,
                   tab.v = tab.v,
                   tab.w = tab.w,
                   tab.loss = tab.loss
                   )

    return(invisible(output))



  }

