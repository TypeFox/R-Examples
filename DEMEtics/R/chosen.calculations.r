chosen.calculations <- function(bias,format.table,pm,statistics,bt=1000,x){

  out.a <- print("you chose tho evaluate:")


  
  if (bias=="correct"){out.b <- cat("bias corrected")
                     }else(if(bias=="uncorrected"){out.b <- cat("non-bias corrected")
                                                 }else{out.b <- cat("FAILURE in 'bias' argument!")})
                         

  out.c <- paste(x,"values")



  ifelse(format.table==TRUE,out.d <- cat("table format is transformed"),out.d <- cat("table format is not transformed"))

  if (pm=="pairwise"){out.e <- cat("pairwise between populations")
                    }else(if(pm=="overall"){out.e <- cat("averaged over all populations")
                                          }else{out.e <- cat("FAILURE in 'pm' argument!")})

  if (statistics=="all"){out.f <- cat("p-values and confidence intervals ","based on ",paste(bt," bootstrap resamplings"))
                       }else(if(statistics=="p"){out.f <- cat("p-values","based on ",paste(bt," bootstrap resamplings"))
                                               }else (if (statistics=="CI"){out.f <- cat("confidence intervals","based on ",paste(bt," bootstrap resamplings"))
                                                                          }else(if(statistics=="none"){out.f <- cat("no statistics")
                                                                                                     }else{out.f <- cat("FAILURE in 'statistics' argument!")})))

  
  out.all.1 <- paste(out.a,"\n",out.b,out.c,"\n",out.e)
       out.all.2 <- paste("statistical evaluation comprises:",out.f)
}
