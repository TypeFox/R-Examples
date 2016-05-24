# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"find.right.table" <-
  function(object,
           inclZeroWRES,
           onlyfirst,
           samp,
           PI.subset,
           subscripts,
           PI.bin.table,
           panel.number,
           ...
           ){

    tmp.table <- NULL

    ## choose the right conditioning variable
    if(!is.null(samp)) {
      data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,samp=samp,subset=PI.subset)
    } else {
      data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=PI.subset)
    }
    tmp.data <- data[subscripts,] 
    stratas <- PI.bin.table[[length(PI.bin.table)]]
    num.stratas <- length(stratas)
    
    ## first check if the panel.number is the same as the VPC strata number
    tmp.strata=stratas[panel.number]

    dim.sub.data <- dim(subset(tmp.data,eval(parse(text=tmp.strata))))

    if (dim.sub.data[1] == dim(tmp.data)[1]){
      tmp.table <- PI.bin.table[[panel.number]]
    } else {
      cat(paste("The conditioning variable for the plot\n"))
      cat(paste("  and the conditioning variable from the VPC file\n"))
      cat(paste("  are not in the same order. Searching for the right \n"))
      cat(paste("  prediction interval values to use \n"))
      cat(paste("\n"))

      for(i in 1:num.stratas){
        if (is.null(tmp.table)){
          tmp.strata=stratas[i]
          tmp.strata=gsub(" = ", " == ",tmp.strata)
          dim.sub.data <- dim(subset(tmp.data,eval(parse(text=tmp.strata))))
          if (dim.sub.data[1] == dim(tmp.data)[1]){
            tmp.table <- PI.bin.table[[i]]
          }
        }
      }
    }
    return(tmp.table)
  }

"setup.PPI" <-
  function(PIlimits,
           PI.mirror,
           tmp.table,
           ...
           ){

    #browser()
    ci.indx <- grep("*.CI.*",names(tmp.table))
    #names(tmp.table)[[ci.indx[[1]]]]
    ci.val <- sub("(\\d*)\\.CI\\..*","\\1",names(tmp.table)[[ci.indx[[1]]]],perl=TRUE)

    sim.bin.table.cols <- paste(PIlimits*100,"sim",sep=".")
    real.bin.table.cols <- paste(PIlimits*100,"real",sep=".")
    sim.ci.upper.bin.table.cols <- paste(ci.val,"CI.for",PIlimits*100,"to",sep=".")
    sim.ci.lower.bin.table.cols <- paste(ci.val,"CI.for",PIlimits*100,"from",sep=".")

    sim.bin.table.cols.50 <- paste(50,"sim",sep=".")
    real.bin.table.cols.50 <- paste(50,"real",sep=".")
    sim.ci.upper.bin.table.cols.50 <- paste(ci.val,"CI.for",50,"to",sep=".")
    sim.ci.lower.bin.table.cols.50 <- paste(ci.val,"CI.for",50,"from",sep=".")
    
    sim.bin.table.cols.mean <- paste("mean","sim",sep=".")
    real.bin.table.cols.mean <- paste("mean","real",sep=".")
    sim.ci.upper.bin.table.cols.mean <- paste(ci.val,"CI.for","mean","to",sep=".")
    sim.ci.lower.bin.table.cols.mean <- paste(ci.val,"CI.for","mean","from",sep=".")

    sim.bin.table.cols.delta.mean <- paste("delta.mean","sim",sep=".")
    real.bin.table.cols.delta.mean <- paste("delta.mean","real",sep=".")
    sim.ci.upper.bin.table.cols.delta.mean <- paste(ci.val,"CI.for","delta.mean","to",sep=".")
    sim.ci.lower.bin.table.cols.delta.mean <- paste(ci.val,"CI.for","delta.mean","from",sep=".")

    mir.bin.table.cols <- NULL
    mir.bin.table.cols.50 <- NULL
    mir.bin.table.cols.mean <- NULL
    mir.bin.table.cols.delta.mean <- NULL
    mir.names.lower <- NULL
    mir.names.upper <- NULL
    mir.names.median <- NULL
    mir.names.mean <- NULL
    mir.names.delta.mean <- NULL
    
    if (!is.null(PI.mirror)) {
      ## what sort of mirror do we have?
      if(is.logical(PI.mirror)) {
        PI.mirror <- 1
      }
      mir.list <- c()
      for (j in 1:PI.mirror){
        mir.list <- c(mir.list,paste("mirror",j,sep="."))
      }
      final.mir.list <- c()
      for (j in 1:length(PIlimits)){
        final.mir.list <- c(final.mir.list,
                            paste(PIlimits[j]*100,mir.list,sep=".")
                            )
      }
      mir.bin.table.cols <- final.mir.list
      mir.bin.table.cols.50 <- paste("50",mir.list,sep=".")
      mir.bin.table.cols.mean <- paste("mean",mir.list,sep=".")
      mir.bin.table.cols.delta.mean <- paste("delta.mean",mir.list,sep=".")
      mir.names.lower <- paste(mir.list,"lower",sep=".")
      mir.names.upper <- paste(mir.list,"upper",sep=".")
      mir.names.median <- paste(mir.list,"median",sep=".")
      mir.names.mean <- paste(mir.list,"mean",sep=".")
      mir.names.delta.mean <- paste(mir.list,"delta.mean",sep=".")
    }


    PPI <- tmp.table[c(sim.bin.table.cols,
                       sim.bin.table.cols.50,
                       sim.ci.lower.bin.table.cols,
                       sim.ci.upper.bin.table.cols,
                       sim.ci.lower.bin.table.cols.50,
                       sim.ci.upper.bin.table.cols.50,
                       real.bin.table.cols,
                       real.bin.table.cols.50,
                       mir.bin.table.cols,
                       mir.bin.table.cols.50,
                       "lower","upper")]
    
    names(PPI) <- c("lower","upper","median",
                    "lower.ci.lower","upper.ci.lower",
                    "lower.ci.upper","upper.ci.upper",
                    "median.ci.lower","median.ci.upper",
                    "real.lower","real.upper","real.median",
                    mir.names.lower,mir.names.upper,mir.names.median,
                    "Xlower","Xupper")
    
    if(length(grep("mean",names(tmp.table)))!=0){
        PPI["mean"] <- tmp.table[c(sim.bin.table.cols.mean)]
        PPI["mean.ci.lower"] <- tmp.table[c(sim.ci.lower.bin.table.cols.mean)]
        PPI["mean.ci.upper"] <- tmp.table[c(sim.ci.upper.bin.table.cols.mean)]
        PPI["real.mean"] <- tmp.table[c(real.bin.table.cols.mean)]
        PPI["mir.names.mean"] <- tmp.table[c(mir.bin.table.cols.mean)]
    } 
    if(length(grep("delta.mean",names(tmp.table)))!=0){
        PPI["delta.mean"] <- tmp.table[c(sim.bin.table.cols.delta.mean)]
        PPI["delta.mean.ci.lower"] <- tmp.table[c(sim.ci.lower.bin.table.cols.delta.mean)]
        PPI["delta.mean.ci.upper"] <- tmp.table[c(sim.ci.upper.bin.table.cols.delta.mean)]
        PPI["real.delta.mean"] <- tmp.table[c(real.bin.table.cols.delta.mean)]
        PPI["mir.names.delta.mean"] <- tmp.table[c(mir.bin.table.cols.delta.mean)]
    } 

    return(PPI)
  }


"get.polygon.regions" <-
  function(PPI,
           PI.mirror,
           ...
           ){
    
    XU <- PPI$Xupper
    XL <- PPI$Xlower
    YU <- PPI$upper
    YL <- PPI$lower
                                        #YMean <- PPI$mean
    Ymed <- PPI$median
    if(length(grep("mean",names(PPI)))!=0) Ymean <- PPI$mean
    if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.mean <- PPI$delta.mean
    
    YUU <- PPI$upper.ci.upper
    YUL <- PPI$upper.ci.lower
    YLU <- PPI$lower.ci.upper
    YLL <- PPI$lower.ci.lower
    YMU <- PPI$median.ci.upper
    YML <- PPI$median.ci.lower
    if(length(grep("mean",names(PPI)))!=0){
        YmeanU <- PPI$mean.ci.upper
        YmeanL <- PPI$mean.ci.lower
    }
    if(length(grep("delta.mean",names(PPI)))!=0){
        Ydelta.meanU <- PPI$delta.mean.ci.upper
        Ydelta.meanL <- PPI$delta.mean.ci.lower
    }
    
    YUR <- PPI$real.upper
    YLR <- PPI$real.lower
    YmedR <- PPI$real.median
    if(length(grep("mean",names(PPI)))!=0) YmeanR <- PPI$real.mean
    if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.meanR <- PPI$real.delta.mean
    
    if (!is.null(PI.mirror)) {
      YUM <- PPI[grep("mirror.*upper",names(PPI))]
      YLM <- PPI[grep("mirror.*lower",names(PPI))]
      YmedM <- PPI[grep("mirror.*median",names(PPI))]
      if(length(grep("mean",names(PPI)))!=0) YmeanM <- PPI[grep("mirror.*mean",names(PPI))]
      if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.meanM <- PPI[grep("mirror.*delta.mean",names(PPI))]
                                        #YUM <- PPI[mir.names.upper]
                                        #YLM <- PPI[mir.names.lower]
                                        #YmedM <- PPI[mir.names.median]
    }
    
    ## Niclas method
    ##tmpx <- c(XU,rev(XU))
    ##tmpy <- c(YU,rev(YL))
    
    if(all(is.na(XL))){ # there are points and not bins
      x.recs <- c(XU,rev(XU))
      y.recs <- c(YU,rev(YL))
      y.up.recs <- c(YUU,rev(YUL))
      y.down.recs <- c(YLU,rev(YLL))
      y.med.recs <- c(YMU,rev(YML))
      if(length(grep("mean",names(PPI)))!=0) y.mean.recs <- c(YmeanU,rev(YmeanL))
      if(length(grep("delta.mean",names(PPI)))!=0) y.delta.mean.recs <- c(Ydelta.meanU,rev(Ydelta.meanL))
    } else { # there are bins
      YU.rec <- YU
      YL.rec <- YL
      YUU.rec <- YUU
      YUL.rec <- YUL
      YLU.rec <- YLU
      YLL.rec <- YLL
      YMU.rec <- YMU
      YML.rec <- YML
      if(length(grep("mean",names(PPI)))!=0){
          YmeanU.rec <- YmeanU
          YmeanL.rec <- YmeanL
      }
      if(length(grep("delta.mean",names(PPI)))!=0){
          Ydelta.meanU.rec <- Ydelta.meanU
          Ydelta.meanL.rec <- Ydelta.meanL
      }
      XU.rec <- XU
      XL.rec <- XL
      
      ## adjust bins of zero length
      if(any(XL==XU)){
        for(i in 1:length(XL)){
          if(XL[i]==XU[i]){
            if(i!=1 & i!=length(XL)){
              XU.rec[i] <- XU[i]+0.05*(XU[i+1]-XU[i])
              XL.rec[i+1] <- XU.rec[i]
              XL.rec[i] <- XL[i]-0.05*(XL[i]-XL[i-1])
              XU.rec[i-1] <- XL.rec[i]
            }
            if(i==1){
              XU.rec[i] <- XU[i]+0.1*(XU[i+1]-XU[i])
              XL.rec[i+1] <- XU.rec[i]
            }
            if(i==length(XL)){
              XL.rec[i] <- XL[i]-0.1*(XL[i]-XL[i-1])
              XU.rec[i-1] <- XL.rec[i]
            }
          }
        }
      }
      
      change.pt <- sort(c(XL,XU))
      X.tot <- NULL
      YU.tot <- NULL
      YL.tot <- NULL
      YU.cur <- NULL
      YL.cur <- NULL
      for(i in 2:length(change.pt)){
        int.val.x <- (change.pt[i]-change.pt[i-1])/2 
        ##if(change.pt[i]==change.pt.cur){
                                        #  X.tot <- c(X.tot,change.pt[i])
                                        #}
        ## get x values
        X.tot <- c(X.tot,change.pt[i])
        
        ## if(i==1 | i==length(change.pt)){
        ##             X.tot <- c(X.tot,change.pt[i])
        ##           } else {
        ##             X.tot <- c(X.tot,change.pt[i])
        ##             YU.tot <- c(YU.tot,YU.cur)
        ##             YL.tot <- c(YL.tot,YL.cur)
        ##           }
        
        ## get y values
        ## check if bin is on or not
        bins.on <- (change.pt[i]<=XU & change.pt[i]>XL)
        if(i==1) bins.on[1] <- (change.pt[i]<=XU[1] & change.pt[i]>=XL[1])
        
        YU.cur <- mean(YU[bins.on])
        YL.cur <- mean(YL[bins.on])
        YU.tot <- c(YU.tot,YU.cur)
        YL.tot <- c(YL.tot,YL.cur)
      }
      
      x.recs <- c(X.tot,rev(X.tot))
      y.recs <- c(YU.tot,rev(YL.tot))
      
      x.recs <- as.vector(t(cbind(XL.rec,XU.rec,XU.rec,XL.rec,NA)))
      y.recs <- as.vector(t(cbind(YU.rec,YU.rec,YL.rec,YL.rec,NA)))
      
      y.up.recs <- as.vector(t(cbind(YUU.rec,YUU.rec,YUL.rec,YUL.rec,NA)))
      y.down.recs <- as.vector(t(cbind(YLU.rec,YLU.rec,YLL.rec,YLL.rec,NA)))
      y.med.recs <- as.vector(t(cbind(YMU.rec,YMU.rec,YML.rec,YML.rec,NA)))
      if(length(grep("mean",names(PPI)))!=0){
          y.mean.recs <- as.vector(t(cbind(YmeanU.rec,YmeanU.rec,YmeanL.rec,YmeanL.rec,NA)))
      }
      if(length(grep("delta.mean",names(PPI)))!=0){
          y.delta.mean.recs <- as.vector(t(cbind(Ydelta.meanU.rec,Ydelta.meanU.rec,Ydelta.meanL.rec,Ydelta.meanL.rec,NA)))
      }
          
    }

    ret <- list(x.recs=x.recs,
                y.recs=y.recs,
                y.up.recs=y.up.recs,
                y.down.recs=y.down.recs,
                y.med.recs=y.med.recs)

    if(length(grep("mean",names(PPI)))!=0) ret <- c(ret, "y.mean.recs"=list(y.mean.recs))
    if(length(grep("delta.mean",names(PPI)))!=0) ret <- c(ret, "y.delta.mean.recs"=list(y.delta.mean.recs))
    
    return(ret)
  }
    
