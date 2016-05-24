# GUI for outlier detector.
# TODO: polish code - this code is quick'n'dirty to put it mildly.
# mvdl, 14.05.2010
# 03.11.2011 : Changed dependency on gWidgets and gWidgetstcltk.
# NOTE: This gui will be replaced by a more decent one based on gtk.

evGui <- function(y){
## setup
   distributions <- c(normal="normal",
                      pareto="pareto",
                      exponential="exponential",
                      lognormal="lognormal",
                      weibull="weibull")
   methods <- c(methodI="I",
                methodII="II")
   plotModes <- c(qq="qq", 
                  residual="residual")
   plotFat <- c(Off=FALSE, On=TRUE)


## plot updater
   updatePlot <- function(h, ...){
      if ( svalue(slFLeft) >= svalue(slFRight) ){
            gmessage("Left quantile limit >= Right quantile limit.\n Resetting default values",
               icon="error",parent=win)
            svalue(slFLeft) <- 10
            svalue(slFRight) <- 90
            return()
        }
      method <- methods[svalue(rdMethod)]
      if ( method == methods[1] ){
         enabled(tbAlpha) <- FALSE 
         enabled(tbRho) <- TRUE
         svalue(cbResidual) <- FALSE
         enabled(cbResidual) <- FALSE
         rhoOrAlpha <- c(svalue(rhoLeft),svalue(rhoRight))
     } else if ( method == methods[2] ){
         enabled(tbRho) <- FALSE
         enabled(tbAlpha) <- TRUE
         enabled(cbResidual) <- TRUE
         rhoOrAlpha <- c(svalue(alphaLeft),svalue(alphaRight))
     }
     # note: this line uses argument order. Could be more robust with string parsing method.
     L <- getOutliers(y, method=method,
               rhoOrAlpha, 
               distribution=distributions[svalue(distribution)],
               FLim=c(svalue(slFLeft)/100, svalue(slFRight)/100 )
               )
      mode<-"qq"
      if ( svalue(cbResidual) ){ mode<-"residual" } 
      outlierPlot(y, L, mode=mode, fat=plotFat[svalue(rdFat)])
      # update result block
      # empty fields for formatting purpose
        tblResult[2,1] <- paste("                                ")
        tblResult[3,1] <- paste("                                ")
        tblResult[2,2] <- paste("                                ")
        tblResult[3,2] <- paste("                                ")

        tblResult[2,2] <- paste("Left ", L$nOut[1])
        tblResult[3,2] <- paste("Right", L$nOut[2])
      switch(distributions[svalue(distribution)],
        normal = { tblResult[2,1] <-    paste("mu          ", round(L$mu,4))
                   tblResult[3,1] <-    paste("sigma       ", round(L$sigma,4))},
        pareto = { tblResult[2,1] <-    paste("ym          ", round(L$ym,4) )
                   tblResult[3,1] <-    paste("alpha       ", round(L$alpha,4))},
        exponential = {tblResult[2,1] <-paste("lambda      ",round(L$lambda,4))   
                       tblResult[3,1] <-paste("                              ") },
        lognormal = { tblResult[2,1] <- paste("mu(ln y)    ",round(L$mu,4))
                      tblResult[3,1] <- paste("sigma(ln y) ", round(L$sigma,4))},
        weibull = { tblResult[2,1] <-   paste("lambda      ",round(L$lambda,4))
                    tblResult[3,1] <-   paste("k           ",round(L$k,4))}
      )
       

   }

## Widgets
   win <- gwindow("Outlier detector")

   modelgp <- gframe("Model parameters",horizontal=FALSE, container=win)
   
   tbld <- glayout(container=modelgp)
   distribution <- gcombobox(names(distributions),
                     container=tbld,
                     handler=updatePlot)
   tbld[1,1] <- "Model distribution:"
   tbld[1,2] <- distribution
   
   tblp <- glayout(container=modelgp)
   slFLeft <- gslider(from=0, to=100, value=10, 
                     container=tblp,
                     handler=updatePlot)
   slFRight <- gslider(from=0, to=100, value=90, 
                     container=tblp,
                     handler=updatePlot)
   tblp[1,1] = "Left quantile limit"
   tblp[2,1] = slFLeft
   tblp[1,2] = "Right quantile limit"
   tblp[2,2] = slFRight

# detection parameters
   frDetect <- gframe("Detection parameters", horizontal=FALSE, container=win)

# method picker
   tbMethod <- glayout(container=frDetect)
   rdMethod <- gcombobox(names(methods),
                     container=tbMethod,
                     handler=updatePlot)
   tbMethod[1,1] <- "Detection method"
   tbMethod[2,1] <- rdMethod

# rho spinners
   tbRho <- glayout(container=frDetect)
   rhoLeft <- gspinbutton(from=0.1, to=length(y), by=0.1,
                      value=1,
                      container=tbRho, 
                      handler=updatePlot)
   
   rhoRight <- gspinbutton(from=0.1, to=length(y), by=0.1,
                       value=1,
                       container=tbRho, 
                       handler=updatePlot)
   tbRho[1,1] <- "rho- (left outliers)"
   tbRho[2,1] <- rhoLeft
   tbRho[1,2] <- "rho+ (right outliers)"
   tbRho[2,2] <- rhoRight

# alpha spinners
   tbAlpha <- glayout(container=frDetect)
   alphaLeft <- gspinbutton(from=0.1,to=1,by=0.01,
                        value=0.05,
                        container=tbAlpha, 
                        handler=updatePlot,)
   
   alphaRight <- gspinbutton(from=0.01,to=1,by=0.01,
                         value=0.05,
                         container=tbAlpha, 
                         handler=updatePlot)
   tbAlpha[1,1] <- "alpha- (left outliers)"
   tbAlpha[2,1] <- alphaLeft
   tbAlpha[1,2] <- "alpha+ (right outliers)"
   tbAlpha[2,2] <- alphaRight

# plot options
   frPlot <- gframe("Plot options", horizontal=TRUE, container=win)
   tbPlot <- glayout(container=frPlot)
   rdFat <- gcombobox(names(plotFat),
               handler=updatePlot,
               container=tbPlot)
   tbPlot[1,1] <- "Fat plot"
   tbPlot[1,2] <- rdFat
   
   cbResidual <- gcheckbox("Show residuals", 
                  handler=updatePlot, 
                  container=frPlot)

# Button with script generator
  # handler
  showScript <- function(h, ...){
      detParStr <- paste("rho=c(", svalue(rhoLeft), ",", svalue(rhoRight),")")
      if ( methods[svalue(rdMethod)] == methods[2] )
         detParStr <- paste("alpha=c(",svalue(alphaLeft),",",svalue(alphaRight),")")
      
      outlierString <- paste(
      "L <- getOutliers(<your variable>, method=","'",methods[svalue(rdMethod)],"'",  ", ",detParStr,
      ", distribution=","'",distributions[svalue(distribution)],"', ",
      "FLim=c(",svalue(slFLeft)/100,",",svalue(slFRight)/100,"))",sep="")

      mode<-"qq"
      if ( svalue(cbResidual) ){ mode<-"residual" } 

      plotString <- paste("outlierPlot(<your variable>, L, mode=","'",mode,"')",sep="")
      
      gtext(text=paste(outlierString,plotString,sep="\n\n"), container=gwindow("Code"))
      }

  
  # the actual button
  frResult <- gframe("Results", container=win)
  btScript <- gbutton("Code",
                      container=frResult,
                      handler=showScript)
  
  tblResult <- glayout(homgeneous=TRUE,container=frResult)
  tblResult[1,1] <- "Parameters"
  tblResult[1,2] <- "Outliers"

  # initialize plot.
   updatePlot()
}














