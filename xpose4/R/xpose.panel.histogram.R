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

xpose.panel.histogram <- function(x,
                                  object,
                                  ##data,
                                  ##subscripts,
                                  ##inclZeroWRES = FALSE,
                                  ##onlyfirst = FALSE,
                                  ##samp = NULL,
                                  ##xvarnam=NULL,
                                  breaks=NULL,
                                  dens=TRUE, # density plot on top of histogram?
                                  hidlty = object@Prefs@Graph.prefs$hidlty, 
                                  hidcol = object@Prefs@Graph.prefs$hidcol, 
                                  hidlwd = object@Prefs@Graph.prefs$hidlwd, 
                                  hiborder = object@Prefs@Graph.prefs$hiborder, 
                                  hilty = object@Prefs@Graph.prefs$hilty, 
                                  hicol = object@Prefs@Graph.prefs$hicol, 
                                  hilwd = object@Prefs@Graph.prefs$hilwd,
                                  math.dens=NULL,
                                  
                                  ## vline settings
                                  vline= NULL,#object@Prefs@Graph.prefs$abline,
                                  vllwd= 3,#object@Prefs@Graph.prefs$abllwd,
                                  vllty= 1,#object@Prefs@Graph.prefs$abllty,
                                  vlcol= "grey",#object@Prefs@Graph.prefs$ablcol,

                                  ## hline settings
                                  hline= NULL,#object@Prefs@Graph.prefs$abline,
                                  hllwd= 3,#object@Prefs@Graph.prefs$abllwd,
                                  hllty= 1,#object@Prefs@Graph.prefs$abllty,
                                  hlcol= "grey",#object@Prefs@Graph.prefs$ablcol,

                                  ## allow for different bins in different panels for continuous data
                                  bins.per.panel.equal = TRUE,

                                  showMean = FALSE,
                                  meanllwd= 3,#object@Prefs@Graph.prefs$abllwd,
                                  meanllty= 1,#object@Prefs@Graph.prefs$abllty,
                                  meanlcol= "orange",#object@Prefs@Graph.prefs$ablcol,
                                  
                                  showMedian = FALSE,
                                  medianllwd= 3,#object@Prefs@Graph.prefs$abllwd,
                                  medianllty= 1,#object@Prefs@Graph.prefs$abllty,
                                  medianlcol= "black",#object@Prefs@Graph.prefs$ablcol,

                                  showPCTS = FALSE,
                                  PCTS = c(0.025,0.975), # vector of percentiles to plot, can be any length
                                  PCTSllwd= 2,#object@Prefs@Graph.prefs$abllwd,
                                  PCTSllty= hidlty,#object@Prefs@Graph.prefs$abllty,
                                  PCTSlcol= "black",#object@Prefs@Graph.prefs$ablcol,

                                  ## vline settings different for each histogram
                                  vdline= NULL,#object@Prefs@Graph.prefs$abline,
                                  vdllwd= 3,#object@Prefs@Graph.prefs$abllwd,
                                  vdllty= 1,#object@Prefs@Graph.prefs$abllty,
                                  vdlcol= "red",#object@Prefs@Graph.prefs$ablcol,
                                  ...,
                                  groups) {
  
  ## if(!is.null(samp)) {
##     data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,samp=samp)
##   } else {
##     data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst)
##   }
  if(length(unique(x)) <= object@Prefs@Cat.levels){
    x <- as.factor(x)
  }

  if(is.factor(x)) {
      nint   <- length(levels(x))
      breaks <- seq(0.5, length = length(levels(x))+1)
  } else {
      if(!bins.per.panel.equal){
          nint      <- round(log2(length(x))+1)
          endpoints <- range(x[!is.na(x)])
          #if(is.null(breaks)) breaks <- do.breaks(endpoints, nint)
          breaks <- do.breaks(endpoints, nint)
      }
  }

  panel.histogram(x,
                  breaks=breaks,
                  lty = hilty,
                  lwd = hilwd,
                  col = hicol,
                  border = hiborder,
                  #type = "density",
                  ...
                  )
  
  if(dens){
      if (is.numeric(x)) {## this should be a choice... to plot this not required
          panel.densityplot(x,
                                        #breaks=breaks,
                            lty=hidlty,
                            col=hidcol,
                            lwd=hidlwd,
                            ...)
          
          ## dens <- density(x)
          ##     panel.xyplot(dens$x, dens$y,
          ##                  type="l",
          ##                  lty=hidlty,
          ##                  col=hidcol,
          ##                  lwd=hidlwd,
          ##                  ...
          ##                  )
      }
  }
  if (!is.null(math.dens)){
    panel.mathdensity(dmath = dnorm,
                      args = list(mean=math.dens$mean,sd=math.dens$sd),
                      col.line="black",lwd=3,
                      ...)
  }

  ## vertical Line in histogram
  if(!is.null(vline)) {
    panel.abline(v=vline,col=vlcol,lwd=vllwd,lty=vllty)
  }

  ## Horizontal Line in histogram
  if(!is.null(hline)) {
    panel.abline(h=hline,col=hlcol,lwd=hllwd,lty=hllty)
  }

   ## Horizontal Line in histogram, different for each histogram
  if(!is.null(hline)) {
    panel.abline(h=hline,col=hlcol,lwd=hllwd,lty=hllty)
  }

  if(showMean | showMedian) sp <- summary(x)
  if(showMean) panel.abline(v=sp[4],col=meanlcol,lwd=meanllwd,lty=meanllty)
  if(showMedian) panel.abline(v=sp[3],col=medianlcol,lwd=medianllwd,lty=medianllty)
  if(showPCTS){
      qu <- quantile(x, PCTS, na.rm=T)
      panel.abline(v=qu,col=PCTSlcol,lwd=PCTSllwd,lty=PCTSllty)
  }

  if(!is.null(vdline)) {
      panel.abline(v=vdline[groups = panel.number()], 
                    col=vdlcol,lwd=vdllwd,lty=vdllty)
  }
  
}

  
