path.plot <-
function(
                      synth.res = NA,
                      dataprep.res = NA,
                      tr.intake = NA,
                      Ylab = c("Y Axis"),
                      Xlab = c("Time"),
                      Ylim = NA ,
                      Legend=c("Treated","Synthetic"),
                      Legend.position=c("topright"),
                      Main = NA,
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
                     
                     
                     y0plot1 <- dataprep.res$Y0plot %*% synth.res$solution.w                
                                       
                     # Get Ylim right 
                     if(sum(is.na(Ylim))>0)
                      {
                        Y.max <-  max(c(y0plot1,dataprep.res$Y1plot))
                        Y.min <-  min(c(y0plot1,dataprep.res$Y1plot))
                     
                        Ylim <- c(
                                   (Y.min - .3*Y.min ),
                                   (.3*Y.max + Y.max)
                                  )
                      }
                              
                      
                     plot(
                      dataprep.res$tag$time.plot
                     ,dataprep.res$Y1plot,
                     t="l",
                     col="black",
                     lwd=2,
                     main=Main,
                     ylab=Ylab,
                     xlab=Xlab,xaxs="i",yaxs="i",ylim=Ylim)
                     
                     lines(
                           dataprep.res$tag$time.plot,
                           y0plot1 ,
                           col="black",
                           lty="dashed",
                           lwd=2,
                           cex=4/5
                           )
                                  
                 } else {
                 
                 z0plot <- dataprep.res$Z0 %*% synth.res$solution.w                
                                       
                     # Get Ylim right 
                     if(sum(is.na(Ylim))>0)
                      {
                        Y.max <-  max(c(z0plot,dataprep.res$Z1))
                        Y.min <-  min(c(z0plot,dataprep.res$Z1))
                        
                        Ylim <- c(
                                   (Y.min - .3*Y.min),
                                   (.3*Y.max + Y.max)
                                  )
                      }
                              
                      
                     plot(
                     dataprep.res$tag$time.optimize.ssr
                     ,z0plot,
                     t="l",
                     col="black",
                     lwd=2,
                     main=Main,
                     ylab=Ylab,
                     xlab=Xlab,xaxs="i",yaxs="i",ylim=Ylim)
                     
                     lines(
                           dataprep.res$tag$time.optimize.ssr,
                           dataprep.res$Z1,
                           col="black",
                           lty="dashed",
                           lwd=2,
                           cex=4/5
                           ) 

            }

            abline(v=tr.intake,lty=3,col="black",lwd=2)
            
            if(sum(is.na(Legend))==0)
                     {
                     legend(Legend.position,legend=Legend,
                     lty=1:2,col=c("black","black"),lwd=c(2,2),cex=6/7)
                     }
            
     }

