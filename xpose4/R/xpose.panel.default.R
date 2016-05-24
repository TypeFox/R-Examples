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

"xpose.panel.default" <-
  function(x, y,object,
           subscripts,
           groups = object@Prefs@Xvardef$id,
           grp.col = NULL,
           iplot = NULL,
           inclZeroWRES = FALSE,
           onlyfirst = FALSE,
           samp = NULL,
           
           xvarnam = NULL,
           yvarnam = NULL,
           
           ##xp.xlim = NULL,
           ##xp.ylim = NULL,
           
           ###############################
           ## Prediction interval settings
           ###############################
           PI      = NULL,
           PI.subset=NULL,
           PI.bin.table=NULL,
           PI.real=NULL,  # can be NULL/TRUE
           PI.mirror=NULL,
           PI.ci = NULL,
           PPI = NULL,
           PI.mean = FALSE, # Should the mean be plotted in the VPCs
           PI.delta.mean = FALSE, # Should the delta mean be plotted in the VPCs
           
           PI.limits= c(0.025, 0.975),#object@Prefs@Graph.prefs$PI.limits,
           
           PI.arcol = "lightgreen",#object@Prefs@Graph.prefs$PI.arcol,
           
           PI.up.lty = 2,#object@Prefs@Graph.prefs$PI.up.lty,
           PI.up.type = "l",#object@Prefs@Graph.prefs$PI.up.type,
           PI.up.col = "black",#object@Prefs@Graph.prefs$PI.up.col,
           PI.up.lwd = 2,#object@Prefs@Graph.prefs$PI.up.lwd,
           
           PI.down.lty = 2,#object@Prefs@Graph.prefs$PI.down.lty,
           PI.down.type = "l",#object@Prefs@Graph.prefs$PI.down.type,
           PI.down.col = "black",#object@Prefs@Graph.prefs$PI.down.col,
           PI.down.lwd = 2,#object@Prefs@Graph.prefs$PI.down.lwd,
           
           PI.med.lty = 1,#object@Prefs@Graph.prefs$PI.med.lty,
           PI.med.type = "l",#object@Prefs@Graph.prefs$PI.med.type,
           PI.med.col = "black",#object@Prefs@Graph.prefs$PI.med.col,
           PI.med.lwd = 2,#object@Prefs@Graph.prefs$PI.med.lwd,
           
           PI.mean.lty = 3,#object@Prefs@Graph.prefs$PI.med.lty,
           PI.mean.type = "l",#object@Prefs@Graph.prefs$PI.med.type,
           PI.mean.col = "black",#object@Prefs@Graph.prefs$PI.med.col,
           PI.mean.lwd = 2,#object@Prefs@Graph.prefs$PI.med.lwd,
           
           PI.delta.mean.lty = 3,#object@Prefs@Graph.prefs$PI.med.lty,
           PI.delta.mean.type = "l",#object@Prefs@Graph.prefs$PI.med.type,
           PI.delta.mean.col = "black",#object@Prefs@Graph.prefs$PI.med.col,
           PI.delta.mean.lwd = 2,#object@Prefs@Graph.prefs$PI.med.lwd,
           
           PI.real.up.lty = 2,#object@Prefs@Graph.prefs$PI.real.up.lty,
           PI.real.up.type = "l",#object@Prefs@Graph.prefs$PI.real.up.type,
           PI.real.up.col = "red",#object@Prefs@Graph.prefs$PI.real.up.col,
           PI.real.up.lwd = 2,#object@Prefs@Graph.prefs$PI.real.up.lwd,
           
           PI.real.down.lty = 2,#object@Prefs@Graph.prefs$PI.real.down.lty,
           PI.real.down.type = "l",#object@Prefs@Graph.prefs$PI.real.down.type,
           PI.real.down.col = "red",#object@Prefs@Graph.prefs$PI.real.down.col,
           PI.real.down.lwd = 2,#object@Prefs@Graph.prefs$PI.real.down.lwd,
           
           PI.real.med.lty = 1,#object@Prefs@Graph.prefs$PI.real.med.lty,
           PI.real.med.type = "l",#object@Prefs@Graph.prefs$PI.real.med.type,
           PI.real.med.col = "red",#object@Prefs@Graph.prefs$PI.real.med.col,
           PI.real.med.lwd = 2,#object@Prefs@Graph.prefs$PI.real.med.lwd,
           
           PI.real.mean.lty = 3,#object@Prefs@Graph.prefs$PI.real.med.lty,
           PI.real.mean.type = "l",#object@Prefs@Graph.prefs$PI.real.med.type,
           PI.real.mean.col = "red",#object@Prefs@Graph.prefs$PI.real.med.col,
           PI.real.mean.lwd = 2,#object@Prefs@Graph.prefs$PI.real.med.lwd,
           
           PI.real.delta.mean.lty = 3,#object@Prefs@Graph.prefs$PI.real.med.lty,
           PI.real.delta.mean.type = "l",#object@Prefs@Graph.prefs$PI.real.med.type,
           PI.real.delta.mean.col = "red",#object@Prefs@Graph.prefs$PI.real.med.col,
           PI.real.delta.mean.lwd = 2,#object@Prefs@Graph.prefs$PI.real.med.lwd,
           
           
           
           PI.mirror.up.lty = 2,#object@Prefs@Graph.prefs$PI.mirror.up.lty,
           PI.mirror.up.type = "l",#object@Prefs@Graph.prefs$PI.mirror.up.type,
           PI.mirror.up.col = "darkgreen",#object@Prefs@Graph.prefs$PI.mirror.up.col,
           PI.mirror.up.lwd = 1,#object@Prefs@Graph.prefs$PI.mirror.up.lwd,
           
           PI.mirror.down.lty = 2,#object@Prefs@Graph.prefs$PI.mirror.down.lty,
           PI.mirror.down.type = "l",#object@Prefs@Graph.prefs$PI.mirror.down.type,
           PI.mirror.down.col = "darkgreen",#object@Prefs@Graph.prefs$PI.mirror.down.col,
           PI.mirror.down.lwd = 1,#object@Prefs@Graph.prefs$PI.mirror.down.lwd,
           
           PI.mirror.med.lty = 1,#object@Prefs@Graph.prefs$PI.mirror.med.lty,
           PI.mirror.med.type = "l",#object@Prefs@Graph.prefs$PI.mirror.med.type,
           PI.mirror.med.col = "darkgreen",#object@Prefs@Graph.prefs$PI.mirror.med.col,
           PI.mirror.med.lwd = 1,#object@Prefs@Graph.prefs$PI.mirror.med.lwd,
           
           PI.mirror.mean.lty = 3,#object@Prefs@Graph.prefs$PI.mirror.med.lty,
           PI.mirror.mean.type = "l",#object@Prefs@Graph.prefs$PI.mirror.med.type,
           PI.mirror.mean.col = "darkgreen",#object@Prefs@Graph.prefs$PI.mirror.med.col,
           PI.mirror.mean.lwd = 1,#object@Prefs@Graph.prefs$PI.mirror.med.lwd,
           
           PI.mirror.delta.mean.lty = 3,#object@Prefs@Graph.prefs$PI.mirror.med.lty,
           PI.mirror.delta.mean.type = "l",#object@Prefs@Graph.prefs$PI.mirror.med.type,
           PI.mirror.delta.mean.col = "darkgreen",#object@Prefs@Graph.prefs$PI.mirror.med.col,
           PI.mirror.delta.mean.lwd = 1,#object@Prefs@Graph.prefs$PI.mirror.med.lwd,
           
           PI.ci.up.arcol = "blue",
           PI.ci.up.lty = 3,#object@Prefs@Graph.prefs$PIuplty,
           PI.ci.up.type = "l",#object@Prefs@Graph.prefs$PIuptyp,
           PI.ci.up.col = "darkorange",#object@Prefs@Graph.prefs$PI.up.col,
           PI.ci.up.lwd = 2,#object@Prefs@Graph.prefs$PI.up.lwd,
           
           PI.ci.down.arcol = "blue",
           PI.ci.down.lty = 3,#object@Prefs@Graph.prefs$PIdolty,
           PI.ci.down.type = "l",#object@Prefs@Graph.prefs$PIdotyp,
           PI.ci.down.col = "darkorange",#object@Prefs@Graph.prefs$PI.down.col,
           PI.ci.down.lwd = 2,#object@Prefs@Graph.prefs$PI.down.lwd,
           
           PI.ci.med.arcol = "blue",
           PI.ci.med.lty = 4,#object@Prefs@Graph.prefs$PImelty,
           PI.ci.med.type = "l",#object@Prefs@Graph.prefs$PImetyp,
           PI.ci.med.col = "darkorange",#object@Prefs@Graph.prefs$PI.med.col,
           PI.ci.med.lwd = 2,#object@Prefs@Graph.prefs$PI.med.lwd,
           
           
           PI.ci.mean.arcol = "purple",
           PI.ci.mean.lty = 4,#object@Prefs@Graph.prefs$PImelty,
           PI.ci.mean.type = "l",#object@Prefs@Graph.prefs$PImetyp,
           PI.ci.mean.col = "darkorange",#object@Prefs@Graph.prefs$PI.med.col,
           PI.ci.mean.lwd = 2,#object@Prefs@Graph.prefs$PI.med.lwd,
           
           PI.ci.delta.mean.arcol = "purple",
           PI.ci.delta.mean.lty = 4,#object@Prefs@Graph.prefs$PImelty,
           PI.ci.delta.mean.type = "l",#object@Prefs@Graph.prefs$PImetyp,
           PI.ci.delta.mean.col = "darkorange",#object@Prefs@Graph.prefs$PI.med.col,
           PI.ci.delta.mean.lwd = 2,#object@Prefs@Graph.prefs$PI.med.lwd,
           
           PI.ci.area.smooth=FALSE,
           ###############################
           ## end of PI settings
           ###############################
           
           
           ## Basic plot characteristics
           type = object@Prefs@Graph.prefs$type,
           col  = object@Prefs@Graph.prefs$col,
           pch  = object@Prefs@Graph.prefs$pch,
           cex  = object@Prefs@Graph.prefs$cex,
           lty  = object@Prefs@Graph.prefs$lty,
           lwd  = object@Prefs@Graph.prefs$lwd,
           fill = object@Prefs@Graph.prefs$fill,
           
           ## Text label setting
           ids  = NULL,
           idsmode=object@Prefs@Graph.prefs$idsmode,
           idsext =object@Prefs@Graph.prefs$idsext,
           idscex= object@Prefs@Graph.prefs$idscex,
           idsdir= object@Prefs@Graph.prefs$idsdir,
           
           ## abline settings
           abline= object@Prefs@Graph.prefs$abline,
           abllwd= object@Prefs@Graph.prefs$abllwd,
           abllty= object@Prefs@Graph.prefs$abllty,
           ablcol= object@Prefs@Graph.prefs$ablcol,
           
           smooth= object@Prefs@Graph.prefs$smooth, 
           smlwd = object@Prefs@Graph.prefs$smlwd, 
           smlty = object@Prefs@Graph.prefs$smlty, 
           smcol = object@Prefs@Graph.prefs$smcol, 
           smspan= object@Prefs@Graph.prefs$smspan,
           smdegr= object@Prefs@Graph.prefs$smdegr,
           smooth.for.groups=NULL,
           
           lmline= object@Prefs@Graph.prefs$lmline,
           lmlwd = object@Prefs@Graph.prefs$lmlwd ,
           lmlty = object@Prefs@Graph.prefs$lmlty ,
           lmcol = object@Prefs@Graph.prefs$lmcol ,
           
           suline = object@Prefs@Graph.prefs$suline,
           sulwd  = object@Prefs@Graph.prefs$sulwd ,
           sulty  = object@Prefs@Graph.prefs$sulty ,
           sucol  = object@Prefs@Graph.prefs$sucol ,
           suspan = object@Prefs@Graph.prefs$suspan,
           sudegr = object@Prefs@Graph.prefs$sudegr,
           
           ## Layout parameters
           grid = object@Prefs@Graph.prefs$grid,
           logy = FALSE,
           logx = FALSE,
           
           ## Force x variables to be continuous
           force.x.continuous = FALSE,
           
           ## Categorcal x-variable
           bwhoriz  = object@Prefs@Graph.prefs$bwhoriz,
           bwratio  = object@Prefs@Graph.prefs$bwratio,
           bwvarwid = object@Prefs@Graph.prefs$bwvarwid,
           bwdotpch = object@Prefs@Graph.prefs$bwdotpch,
           bwdotcol = object@Prefs@Graph.prefs$bwdotcol,
           bwdotcex = object@Prefs@Graph.prefs$bwdotcex,
           bwreccol = object@Prefs@Graph.prefs$bwreccol,
           bwrecfill= object@Prefs@Graph.prefs$bwrecfill,
           bwreclty = object@Prefs@Graph.prefs$bwreclty,
           bwreclwd = object@Prefs@Graph.prefs$bwreclwd,
           bwumbcol = object@Prefs@Graph.prefs$bwumbcol,
           bwumblty = object@Prefs@Graph.prefs$bwumblty,
           bwumblwd = object@Prefs@Graph.prefs$bwumblwd,
           bwoutcol = object@Prefs@Graph.prefs$bwoutcol,
           bwoutcex = object@Prefs@Graph.prefs$bwoutcex,
           bwoutpch = object@Prefs@Graph.prefs$bwoutpch,
           autocorr=FALSE,
           
           ## vline settings
           vline= NULL,#object@Prefs@Graph.prefs$abline,
           vllwd= 3,#object@Prefs@Graph.prefs$abllwd,
           vllty= 2,#object@Prefs@Graph.prefs$abllty,
           vlcol= "grey",#object@Prefs@Graph.prefs$ablcol,
           
           ## hline settings
           hline= NULL,#object@Prefs@Graph.prefs$abline,
           hllwd= 3,#object@Prefs@Graph.prefs$abllwd,
           hllty= 1,#object@Prefs@Graph.prefs$abllty,
           hlcol= "grey",#object@Prefs@Graph.prefs$ablcol,
           
           #data,
           pch.ip.sp=pch, # ind.plots single point per individual
           cex.ip.sp=cex, # ind.plots single point per individual
           ...
           
  ) {
    
    ## data should already be passed to the function at this point
    ## this should be changed so that we just use the data passed form the
    ## plotting function
    ##     if(!is.null(samp)) {
    ##       data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,samp=samp)
    ##     } else {
    ##       data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst)
    ##     }
    
    #if(force.x.continuous == FALSE) {
    #  if(length(unique(data[subscripts,xvarnam])) <= object@Prefs@Cat.levels) x <- as.factor(x)
    #}
    
    ## Compute and plot prediction areas if requested.
    ## This needs to be performed here for the area  to appear at
    ## the bottom of the rest.
    if(!is.null(PI) |
         !is.null(PI.real) |
         !is.null(PI.mirror) |
         !is.null(PI.ci)
    ){
      if(is.null(PI.bin.table)){
        if(is.null(PPI)){
          PPI <- computePI(xvarnam,yvarnam,object,logy=logy,logx=logx,limits=PI.limits,
                           onlyfirst=onlyfirst,inclZeroWRES=inclZeroWRES,PI.subset)
        }
      } else {
        if(!is.null(dim(PI.bin.table))){ # there is only one table and no conditioning
          tmp.table <- PI.bin.table
        } else {  # There is a stratification variable
          tmp.table <- find.right.table(object,inclZeroWRES,onlyfirst,samp,PI.subset,
                                        subscripts=subscripts,PI.bin.table,
                                        panel.number=panel.number(),...)
          
          if (is.null(tmp.table)){
            cat(paste("No strata in VPC file found to\n"))
            cat(paste("  match conditioning variables\n"))
            cat(paste("\n"))
            return()
          }
        }
        
        ## now set up PPI table
        PPI <- setup.PPI(PI.limits,PI.mirror,tmp.table,...)
        
      }
      
      XU <- PPI$Xupper
      XL <- PPI$Xlower
      YU <- PPI$upper
      YL <- PPI$lower
      if(length(grep("mean",names(PPI)))!=0) Ymean <- PPI$mean
      if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.mean <- PPI$delta.mean
      Ymed <- PPI$median
      
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
      
    }
    
    if((!is.null(PI) && (PI=="area" | PI=="both")) |
         (!is.null(PI.ci) && (PI.ci=="area" | PI.ci=="both"))) {
      
      poly <- get.polygon.regions(PPI,PI.mirror,...)
      if (!is.null(PI) && (PI=="area" | PI=="both")){
        pi.x.recs <- poly$x.recs
        pi.y.recs <- poly$y.recs
        if(logx) {
          tmp <- is.nan(pi.x.recs)
          pi.x.recs <- log10(pi.x.recs)
          tmp2 <- is.nan(pi.x.recs)
          if(any(tmp!=tmp2)){
            cat(paste("The prediction interval on the x-axis goes below zero.",
                      "This means that taking the log of this prediction",
                      "interval gives non-real numbers.",
                      "The plot will not be created.\n",sep="\n"))
            return(NULL)
          }
        }
        if(logy) {
          tmp <- is.nan(pi.y.recs)
          pi.y.recs <- log10(pi.y.recs)
          tmp2 <- is.nan(pi.y.recs)
          if(any(tmp!=tmp2)){
            cat(paste("The prediction interval on the y-axis goes below zero.",
                      "This means that taking the log of this prediction",
                      "interval gives non-real numbers.",
                      "The plot will not be created.\n",sep="\n"))
            return(NULL)
          }
        }
        grid.polygon(pi.x.recs,pi.y.recs,
                     default.units="native",
                     gp=gpar(fill=PI.arcol,col=NULL,lty=0))
      }
      
      if (!is.null(PI.ci) && (PI.ci=="area" | PI.ci=="both")){
        
        if(PI.ci.area.smooth){
          if(all(is.na(XL))){
            XM <- XU
          } else {
            XM <- (XL+XU)/2
            XM <- c(XL[1],XM,XU[length(XU)])
          }
          xrecs <- c(XM,rev(XM))
          y.up.recs <- c(PPI$upper.ci.upper[1],
                         PPI$upper.ci.upper,
                         PPI$upper.ci.upper[dim(PPI)[1]],
                         PPI$upper.ci.lower[dim(PPI)[1]],
                         rev(PPI$upper.ci.lower),
                         PPI$upper.ci.lower[1]
          )
          y.down.recs <- c(PPI$lower.ci.upper[1],
                           PPI$lower.ci.upper,
                           PPI$lower.ci.upper[dim(PPI)[1]],
                           PPI$lower.ci.lower[dim(PPI)[1]],
                           rev(PPI$lower.ci.lower),
                           PPI$lower.ci.lower[1]
          )
          y.med.recs <- c(PPI$median.ci.upper[1],
                          PPI$median.ci.upper,
                          PPI$median.ci.upper[dim(PPI)[1]],
                          PPI$median.ci.lower[dim(PPI)[1]],
                          rev(PPI$median.ci.lower),
                          PPI$median.ci.lower[1]
          )
          
          if(length(grep("mean",names(PPI)))!=0){
            y.mean.recs <- c(PPI$mean.ci.upper[1],
                             PPI$mean.ci.upper,
                             PPI$mean.ci.upper[dim(PPI)[1]],
                             PPI$mean.ci.lower[dim(PPI)[1]],
                             rev(PPI$mean.ci.lower),
                             PPI$mean.ci.lower[1]
            )
          }
          if(length(grep("delta.mean",names(PPI)))!=0){
            y.delta.mean.recs <- c(PPI$delta.mean.ci.upper[1],
                                   PPI$delta.mean.ci.upper,
                                   PPI$delta.mean.ci.upper[dim(PPI)[1]],
                                   PPI$delta.mean.ci.lower[dim(PPI)[1]],
                                   rev(PPI$delta.mean.ci.lower),
                                   PPI$delta.mean.ci.lower[1]
            )
          }
          
        } else {
          xrecs <- poly$x.recs
          y.up.recs <- poly$y.up.recs
          y.down.recs <- poly$y.down.recs
          y.med.recs <- poly$y.med.recs
          if(length(grep("mean",names(PPI)))!=0){
            y.mean.recs <- poly$y.mean.recs
          }
          if(length(grep("delta.mean",names(PPI)))!=0){
            y.delta.mean.recs <- poly$y.delta.mean.recs
          }
        }
        
        if (logx){
          tmp <- is.nan(xrecs)
          xrecs <- log10(xrecs)
          tmp2 <- is.nan(xrecs)
          if(any(tmp!=tmp2)){
            cat(paste("The PI.ci on the x-axis goes below zero.",
                      "This means that taking the log of this prediction",
                      "interval gives non-real numbers.",
                      "The plot will not be created.\n",sep="\n"))
            return(NULL)
          }
        }
        if(logy){
          tmp <- is.nan(c(y.up.recs,y.down.recs,y.med.recs))
          if(length(grep("mean",names(PPI)))!=0) tmp <- is.nan(c(y.up.recs,y.down.recs,y.med.recs,y.mean.recs))
          if(length(grep("delta.mean",names(PPI)))!=0) tmp <- is.nan(c(y.up.recs,y.down.recs,y.med.recs,y.delta.mean.recs))
          y.up.recs <- log10(y.up.recs)
          y.down.recs <- log10(y.down.recs)
          y.med.recs <- log10(y.med.recs)
          if(length(grep("mean",names(PPI)))!=0) y.mean.recs <- log10(y.mean.recs)
          if(length(grep("delta.mean",names(PPI)))!=0) y.delta.mean.recs <- log10(y.delta.mean.recs)
          tmp2 <- is.nan(c(y.up.recs,y.down.recs,y.med.recs))
          if(length(grep("mean",names(PPI)))!=0) tmp2 <- is.nan(c(y.up.recs,y.down.recs,y.med.recs, y.mean.recs))
          if(length(grep("delta.mean",names(PPI)))!=0) tmp2 <- is.nan(c(y.up.recs,y.down.recs,y.med.recs, y.delta.mean.recs))
          if(any(tmp!=tmp2)){
            cat(paste("The PI.ci on the y-axis goes below zero.",
                      "This means that taking the log of this prediction",
                      "interval gives non-real numbers.",
                      "The plot will not be created.\n",sep="\n"))
            return(NULL)
          }
        }
        
        grid.polygon(xrecs,y.up.recs,
                     default.units="native",
                     gp=gpar(fill=PI.ci.up.arcol,alpha=0.3,col=NULL,lty=0)
        )
        grid.polygon(xrecs,y.down.recs,
                     default.units="native",
                     gp=gpar(fill=PI.ci.down.arcol,alpha=0.3,col=NULL,lty=0)
        )
        grid.polygon(xrecs,y.med.recs,
                     default.units="native",
                     gp=gpar(fill=PI.ci.med.arcol,alpha=0.3,col=NULL,lty=0)
        )
        if(PI.mean){
          if(length(grep("mean",names(PPI)))!=0){
            grid.polygon(xrecs,y.mean.recs,
                         default.units="native",
                         gp=gpar(fill=PI.ci.mean.arcol,alpha=0.3,col=NULL,lty=0)
            )
          }
        }
        if(PI.delta.mean){
          if(length(grep("delta.mean",names(PPI)))!=0){
            grid.polygon(xrecs,y.delta.mean.recs,
                         default.units="native",
                         gp=gpar(fill=PI.ci.delta.mean.arcol,alpha=0.3,col=NULL,lty=0)
            )
          }
        }
        
        ##         grid.polygon(poly$x.recs,poly$y.up.recs,
        ##                      default.units="native",
        ##                      gp=gpar(fill=PI.ci.up.arcol,alpha=0.3,col=NULL,lty=0)
        ##                      )
        ##         grid.polygon(poly$x.recs,poly$y.down.recs,
        ##                      default.units="native",
        ##                      gp=gpar(fill=PI.ci.down.arcol,alpha=0.3,col=NULL,lty=0)
        ##                      )
        ##         grid.polygon(poly$x.recs,poly$y.med.recs,
        ##                      default.units="native",
        ##                      gp=gpar(fill=PI.ci.med.arcol,alpha=0.3,col=NULL,lty=0)
        ##                      )
        
      }
    } # end of make polygon
    
    ## Stuff common to both xy and bw
    if(grid != FALSE) {
      panel.grid(h = -1, v = -1)
    }
    
    ## Line of "identity"
    if(!is.null(abline)) {
      panel.abline(abline,col=ablcol,lwd=abllwd,lty=abllty)
    }
    
    ## vertical Line 
    if(!is.null(vline)) {
      panel.abline(v=vline,col=vlcol,lwd=vllwd,lty=vllty)
    }
    
    ## Horizontal Line 
    if(!is.null(hline)) {
      panel.abline(h=hline,col=hlcol,lwd=hllwd,lty=hllty)
    }
    
    ## for autocorrelation
    if(autocorr){
      auto.ids <- unique(groups)
      auto.n <- 0
      xplt1 <- 0
      xplt2 <- 0
      xgrps <- 0
      for(i in 1:length(auto.ids)) {
        seli <- groups == auto.ids[i]
        nobs <- length(x[seli])
        xplt <- matrix(x[seli], 1, nobs)
        if(nobs > 1) {
          for(j in 1:(nobs - 1)) {
            auto.n <- auto.n + 1
            xplt1[auto.n] <- xplt[1, j]
            xplt2[auto.n] <- xplt[1, j + 1]
            xgrps[auto.n] <- auto.ids[i]
          }
        }
      }
      x <- xplt1
      y <- xplt2
      groups <- xgrps
    }
    
    
    ## Plot the data
    if(!is.factor(x) && !bwhoriz) {
      
      if(any(is.null(groups))) {
        
        panel.xyplot(x,y,
                     col   =col,
                     pch   =pch,
                     lty   =lty,
                     type  =type,
                     cex   = cex,
                     lwd   = lwd,
                     fill = fill
        )
      } else {
        ord <- order(x)
        if((any(!is.null(iplot))) || (is.null(grp.col))) {
          if(length(x)==3){
            #             pch[2]=pch.ip.sp[2]
            #             pch[3]=pch.ip.sp[3]
            #             pch[1]=pch.ip.sp[1]
            #             cex[3]=cex.ip.sp[3]
            #             cex[2]=cex.ip.sp[2]
            #             cex[1]=cex.ip.sp[1]
            pch=pch.ip.sp
            cex=cex.ip.sp
          } 
          panel.superpose(x[ord],
                          y[ord],
                          subscripts[ord],
                          col   =col,
                          pch   =pch,
                          cex   = cex,
                          lty   =lty,
                          type  =type,
                          lwd   = lwd,
                          groups=groups,
                          fill = fill
          )
        } else {
          panel.superpose(x[ord],
                          y[ord],
                          subscripts[ord],
                          #col   =col,
                          pch   =pch,
                          cex   = cex,
                          lty   =lty,
                          type  =type,
                          lwd   = lwd,
                          groups=groups,
                          fill=fill
          )       
        }
      }
      
      ## Add a loess smooth?
      if(!any(is.null(smooth))) {
        if(!is.factor(y)){
          if(!any(is.null(smooth.for.groups)) && !any(is.null(groups))) {
            panel.superpose(x,y,subscripts,groups=groups,
                            span  = smspan,
                            degree= smdegr,
                            col   = smcol,
                            lwd   = smlwd,
                            lty   = smlty,
                            panel.groups="panel.loess")
          } else {
            panel.loess(x,y,
                        span  = smspan, # can change this to 0.75 to match R 
                        degree= smdegr,
                        col   = smcol,
                        lwd   = smlwd,
                        lty   = smlty
            )
          }
        } else { # y is a factor
          ##           panel.linejoin(x, y, fun = median, horizontal = TRUE,
          ##                          lwd=smlwd, lty=smlty, col=smcol,
          ##                          col.line=smcol, type=smlty,
          ##                          ...)
        }
      }
      
      ## Add a lm line?
      if(!any(is.null(lmline))) {
        panel.abline(lm(y~x),
                     col   = lmcol,
                     lwd   = lmlwd,
                     lty   = lmlty
        )
      }
      
      ## Add a superpose smooth?
      if(!any(is.null(suline))) {
        ys <- suline[subscripts]
        xs <- x
        if(logy) ys <- log10(ys)
        if(logx) xs <- log10(xs)
        
        panel.loess(xs,ys,
                    span  = suspan, 
                    degree= sudegr,
                    col   = sucol,
                    lwd   = sulwd,
                    lty   = sulty
        )
      }
      
      ## Add id-numbers as plot symbols
      if(!any(is.null(ids))) {
        if (!is.factor(y)){
          ids <- ids[subscripts]
          addid(x,y,ids=ids,
                idsmode=idsmode,
                idsext =idsext,
                idscex = idscex,
                idsdir = idsdir)
        }
      }
      
      ## Compute and plot prediction intervals if requested.
      ## This needs to be performed here for the lines to appear on
      ## top of the rest.
      if(!is.null(PI) && (PI=="lines" | PI=="both")) {
        if(all(is.na(XL))){
          XM <- XU
        } else {
          XM <- (XL+XU)/2
        }
        if(logx) XM <- log10(XM)
        if(logy){
          YU <- log10(YU)
          YL <- log10(YL)
          Ymed <- log10(Ymed)
          if(length(grep("mean",names(PPI)))!=0) Ymean <- log10(Ymean)
          if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.mean <- log10(Ydelta.mean)
        }
        panel.lines(XM,YU,type=PI.up.type,lty=PI.up.lty,col=PI.up.col,lwd=PI.up.lwd)
        panel.lines(XM,YL,type=PI.down.type,lty=PI.down.lty,col=PI.down.col,lwd=PI.down.lwd)
        panel.lines(XM,Ymed,type=PI.med.type,lty=PI.med.lty,col=PI.med.col,lwd=PI.med.lwd)
        if(PI.mean){
          if(length(grep("mean",names(PPI)))!=0){
            panel.lines(XM,Ymean,type=PI.mean.type,lty=PI.mean.lty,col=PI.mean.col,lwd=PI.mean.lwd)
          }
        }
        if(PI.delta.mean){
          if(length(grep("delta.mean",names(PPI)))!=0){
            panel.lines(XM,Ydelta.mean,type=PI.delta.mean.type,lty=PI.delta.mean.lty,col=PI.delta.mean.col,lwd=PI.delta.mean.lwd)
          }
        }
      }
      if(!is.null(PI.real)) {
        if(all(is.na(XL))){
          XM <- XU
        } else {
          XM <- (XL+XU)/2
        }
        if(logx) XM <- log10(XM)
        if(logy){
          YUR <- log10(YUR)
          YLR <- log10(YLR)
          YmedR <- log10(YmedR)
          if(length(grep("mean",names(PPI)))!=0) YmeanR <- log10(YmeanR)
          if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.meanR <- log10(Ydelta.meanR)
        }
        panel.lines(XM,YUR,type=PI.real.up.type,lty=PI.real.up.lty,col=PI.real.up.col,lwd=PI.real.up.lwd)
        panel.lines(XM,YLR,type=PI.real.down.type,lty=PI.real.down.lty,col=PI.real.down.col,lwd=PI.real.down.lwd)
        panel.lines(XM,YmedR,type=PI.real.med.type,lty=PI.real.med.lty,col=PI.real.med.col,lwd=PI.real.med.lwd)
        if(PI.mean){
          if(length(grep("mean",names(PPI)))!=0){
            panel.lines(XM,YmeanR,type=PI.real.mean.type,lty=PI.real.mean.lty,col=PI.real.mean.col,lwd=PI.real.mean.lwd)
          }
        }
        if(PI.delta.mean){
          if(length(grep("delta.mean",names(PPI)))!=0){
            panel.lines(XM,Ydelta.meanR,type=PI.real.delta.mean.type,lty=PI.real.delta.mean.lty,col=PI.real.delta.mean.col,lwd=PI.real.delta.mean.lwd)
          }
        }
      }
      if(!is.null(PI.ci) && (PI.ci=="lines" | PI.ci=="both")) {
        if(all(is.na(XL))){
          XM <- XU
        } else {
          XM <- (XL+XU)/2
        }
        upper.ci.upper <- PPI$upper.ci.upper
        upper.ci.lower <- PPI$upper.ci.lower
        lower.ci.upper <- PPI$lower.ci.upper
        lower.ci.lower <- PPI$lower.ci.lower
        median.ci.upper <- PPI$median.ci.upper
        median.ci.lower <- PPI$median.ci.lower
        if(length(grep("mean",names(PPI)))!=0){
          mean.ci.upper <- PPI$mean.ci.upper
          mean.ci.lower <- PPI$mean.ci.lower
        }
        if(length(grep("delta.mean",names(PPI)))!=0){
          delta.mean.ci.upper <- PPI$delta.mean.ci.upper
          delta.mean.ci.lower <- PPI$delta.mean.ci.lower
        }
        if(logx) XM <- log10(XM)
        if(logy){
          upper.ci.upper <- log10(upper.ci.upper)
          upper.ci.lower <- log10(upper.ci.lower)
          lower.ci.upper <- log10(lower.ci.upper)
          lower.ci.lower <- log10(lower.ci.lower)
          median.ci.upper <- log10(median.ci.upper)
          median.ci.lower <- log10(median.ci.lower)
          if(length(grep("mean",names(PPI)))!=0){
            mean.ci.upper <- log10(mean.ci.upper)
            mean.ci.lower <- log10(mean.ci.lower)
          }
          if(length(grep("delta.mean",names(PPI)))!=0){
            delta.mean.ci.upper <- log10(delta.mean.ci.upper)
            delta.mean.ci.lower <- log10(delta.mean.ci.lower)
          }
        }
        panel.lines(XM,upper.ci.upper,type=PI.ci.up.type,lty=PI.ci.up.lty,col=PI.ci.up.col,lwd=PI.ci.up.lwd)
        panel.lines(XM,upper.ci.lower,type=PI.ci.up.type,lty=PI.ci.up.lty,col=PI.ci.up.col,lwd=PI.ci.up.lwd)
        panel.lines(XM,lower.ci.upper,type=PI.ci.down.type,lty=PI.ci.down.lty,col=PI.ci.down.col,lwd=PI.ci.down.lwd)
        panel.lines(XM,lower.ci.lower,type=PI.ci.down.type,lty=PI.ci.down.lty,col=PI.ci.down.col,lwd=PI.ci.down.lwd)
        panel.lines(XM,median.ci.upper,type=PI.ci.med.type,lty=PI.ci.med.lty,col=PI.ci.med.col,lwd=PI.ci.med.lwd)
        panel.lines(XM,median.ci.lower,type=PI.ci.med.type,lty=PI.ci.med.lty,col=PI.ci.med.col,lwd=PI.ci.med.lwd)
        if(PI.mean){
          if(length(grep("mean",names(PPI)))!=0){
            panel.lines(XM,mean.ci.upper,type=PI.ci.mean.type,lty=PI.ci.mean.lty,col=PI.ci.mean.col,lwd=PI.ci.mean.lwd)
            panel.lines(XM,mean.ci.lower,type=PI.ci.mean.type,lty=PI.ci.mean.lty,col=PI.ci.mean.col,lwd=PI.ci.mean.lwd)
          }
        }
        if(PI.delta.mean){
          if(length(grep("delta.mean",names(PPI)))!=0){
            panel.lines(XM,delta.mean.ci.upper,type=PI.ci.delta.mean.type,lty=PI.ci.delta.mean.lty,col=PI.ci.delta.mean.col,lwd=PI.ci.delta.mean.lwd)
            panel.lines(XM,delta.mean.ci.lower,type=PI.ci.delta.mean.type,lty=PI.ci.delta.mean.lty,col=PI.ci.delta.mean.col,lwd=PI.ci.delta.mean.lwd)
          }
        }
      }
      if(!is.null(PI.mirror)) {
        if(all(is.na(XL))){
          XM <- XU
        } else {
          XM <- (XL+XU)/2
        }
        if(logx) XM <- log10(XM)
        if(logy){
          YUM <- log10(YUM)
          YLM <- log10(YLM)
          YmedM <- log10(YmedM)
          if(length(grep("mean",names(PPI)))!=0) YmeanM <- log10(YmeanM)
          if(length(grep("delta.mean",names(PPI)))!=0) Ydelta.meanM <- log10(Ydelta.meanM) 
        }
        
        for(jj in 1:PI.mirror){
          panel.lines(XM,YUM[[jj]],type=PI.mirror.up.type,lty=PI.mirror.up.lty,col=PI.mirror.up.col,lwd=PI.mirror.up.lwd)
          panel.lines(XM,YLM[[jj]],type=PI.mirror.down.type,lty=PI.mirror.down.lty,col=PI.mirror.down.col,lwd=PI.mirror.down.lwd)
          panel.lines(XM,YmedM[[jj]],type=PI.mirror.med.type,lty=PI.mirror.med.lty,col=PI.mirror.med.col,lwd=PI.mirror.med.lwd)
          if(PI.mean){
            if(length(grep("mean",names(PPI)))!=0){
              panel.lines(XM,YmeanM[[jj]],type=PI.mirror.mean.type,lty=PI.mirror.mean.lty,col=PI.mirror.mean.col,lwd=PI.mirror.mean.lwd)
            }
          }
          if(PI.delta.mean){
            if(length(grep("delta.mean",names(PPI)))!=0){
              panel.lines(XM,Ydelta.meanM[[jj]],type=PI.mirror.delta.mean.type,lty=PI.mirror.delta.mean.lty,col=PI.mirror.delta.mean.col,lwd=PI.mirror.delta.mean.lwd)
            }
          }
        }
        
      }
      
    } else {
      
      oumbset <- trellis.par.get("box.umbrella")
      on.exit(trellis.par.set("box.umbrella",oumbset),add=T)
      umbset     <- oumbset
      umbset$col <- bwumbcol
      umbset$lty <- bwumblty
      umbset$lwd <- bwumblwd
      trellis.par.set("box.umbrella",umbset)
      
      orecset  <- trellis.par.get("box.rectangle")
      on.exit(trellis.par.set("box.rectangle",orecset),add=T)
      recset     <- orecset
      recset$col <- bwreccol
      recset$lty <- bwreclty
      recset$lwd <- bwreclwd
      recset$fill<- bwrecfill
      trellis.par.set("box.rectangle",recset)
      
      ooutset  <- trellis.par.get("plot.symbol")
      on.exit(trellis.par.set("plot.symbol",ooutset),add=T)
      outset     <- ooutset
      outset$col <- bwoutcol
      outset$pch <- bwoutpch
      outset$cex <- bwoutcex
      trellis.par.set("plot.symbol",outset)
      
      panel.bwplot(x,y,
                   horizontal=bwhoriz,
                   col=bwdotcol,
                   pch=bwdotpch,
                   cex=bwdotcex,
                   ratio=bwratio,
                   varwidth=bwvarwid)
    }
    
  }

