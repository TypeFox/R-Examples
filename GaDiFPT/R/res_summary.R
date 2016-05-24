res_summary <- 
  function(obj,Nspikes,fileout)
  {
    
    # ndig = total digits required 
    ndig <- 10
    
    M <- length(Nspikes)
    
    cat('#################################################################### \n')
    cat('#####               First Passage Time Simulation              ##### \n')
    cat('#####              for Gaussian Diffusions: Results            ##### \n')
    cat('#################################################################### \n \n \n')  
    
    
    
    # FPT density statistics
    
    deltat <- obj$time[2]-obj$time[1]
    
    finalmass <- obj$gg0[length(obj$gg0)]
    
    cat ('PROBABILITY MASS REACHED = ', finalmass, '\n \n')
    
    Meantime <- sum(obj$g0*obj$time)*deltat
    cat('Mean time from the density (ms) = ', format(Meantime,digits=ndig), '\n')
    Stdtime <- sqrt(sum(obj$g0*obj$time*obj$time)*deltat-Meantime*Meantime)
    cat('Standard deviation (ms) = ', format(Stdtime,digits=ndig), '\n')
    Medtime <- obj$time[which.min(abs(obj$gg0-0.5))]
    cat('Median time from the distribution (ms) = ', format(Medtime,digits=ndig), '\n \n')
    
    # output to file fileout
    
    results <-cbind(format(obj$time,digits=ndig),format(obj$g0, digits=ndig),format(obj$gg0,digits=ndig),format(obj$g0/(1-obj$gg0),digits=ndig))
    write.table(results, file=fileout,quote=FALSE,
                row.names=FALSE,col.names=c("time","g0","gg0","rate"),sep="  ")
    
    if (M > 0) {
      spikes <- get("spikes")  
      write(spikes,file=fileout,ncolumns=1,append=TRUE)  
    }
    
    cat('\n \n ********** SUMMARY OF RESULTS ********** \n \n', 
        file = fileout, append = TRUE)
    cat('Mean time from the density (ms) = ', format(Meantime,digits=ndig), '\n',file = fileout, append = TRUE)
    cat('Standard deviation (ms) = ', format(Stdtime,digits=ndig), '\n',file = fileout, append = TRUE)
    cat('Median time from the distribution (ms) = ', format(Medtime,digits=ndig), '\n \n \n',file = fileout, append = TRUE)
    
        
    
     if (M > 0) {
    
       # spikes statistics
       cat('Mean time for the spikes (ms) = ', format(mean(spikes),digits=ndig), '\n')
       cat('Standard deviation (ms) = ', format(sd(spikes),digits=ndig), '\n')
       cat('Median time for the spikes (ms) = ', format(median(spikes),digits=ndig), '\n \n')
      
       cat('Mean time for the spikes (ms)   = ', format(mean(spikes),digits=ndig), '\n',file = fileout, append = TRUE)
       cat('Standard deviation (ms)         = ', format(sd(spikes),digits=ndig), '\n',file = fileout, append = TRUE)
       cat('Median time for the spikes (ms) = ', format(median(spikes),digits=ndig), '\n \n',file = fileout, append = TRUE)
       
     }
    }
    
    