gaps.plot <-
function(
                      synth.res = NA,
                      dataprep.res = NA,
                      Ylab = c("Title"),
                      Xlab = c("Time"),
                      Main = c("Gaps: Treated - Synthetic"),
                      tr.intake = NA,
                      Ylim = NA,
                      Z.plot = FALSE
                      
                     )

              { 
                    if(Z.plot==FALSE)
                      { 
                     
                     if(sum(is.na(dataprep.res$Y1plot)) > 0)
                       {
                         stop(
                              "\n\n#####################################################",
                              "\nYou have missing Y data for the treated!\n\n"
                             )  
                        }
                     
                     if(sum(is.na(dataprep.res$Y0plot)) > 0)
                       {
                         stop(
                              "\n\n#####################################################",
                              "\nYou have missing Y data for the controls!\n\n"
                             )  
                        } 
                    
                       
                     # Compute Y-Gap: Treated - Synthtic              
                     gap <- dataprep.res$Y1plot - (dataprep.res$Y0plot %*% synth.res$solution.w)
                     
                     # Get Ylim right 
                     if(sum(is.na(Ylim))>0)
                      {
                        Ylim <- c(
                                  -(.3*max(abs(gap)) + max(abs(gap))),
                                  (.3*max(abs(gap))  + max(abs(gap)))
                                  )
                      }
                                                        
                     plot(
                          dataprep.res$tag$time.plot,
                          gap,
                          t="l",
                          col="black",
                          lwd=2,
                          main=Main,
                          ylab=Ylab,
                          xlab=Xlab,
                          ylim=Ylim,
                          xaxs="i",
                          yaxs="i"
                          )
         
                     } else {
                     
                     # Compute Y-Gap: Treated - Synthtic              
                     gap <- dataprep.res$Z1 - (dataprep.res$Z0 %*% synth.res$solution.w)
                     
                     # Get Ylim right 
                     if(sum(is.na(Ylim))>0)
                      {
                        Ylim <- c(
                                  -(.3*max(abs(gap)) + max(abs(gap))),
                                  (.3*max(abs(gap))  + max(abs(gap)))
                                  )
                      }
                                                        
                     plot(
                          dataprep.res$tag$time.optimize.ssr,
                          gap,
                          t="l",
                          col="black",
                          lwd=2,
                          main=Main,
                          ylab=Ylab,
                          xlab=Xlab,
                          ylim=Ylim,
                          xaxs="i",
                          yaxs="i"
                          )           
                  
                   }           
                     
                     abline(h=0, col="black",lty="dashed",lwd=2)
                     abline(v=tr.intake,col="black",lty="dotted",lwd=2)
           
            }

