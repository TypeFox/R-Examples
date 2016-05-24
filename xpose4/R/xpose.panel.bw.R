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

"xpose.panel.bw" <-
function(x, y, object,
           subscripts,
           groups = NULL,
           inclZeroWRES = FALSE,
           onlyfirst = FALSE,
           samp = NULL,
           
           xvarnam = NULL,
           yvarnam = NULL,
           
           ## Basic plot characteristics
           type = object@Prefs@Graph.prefs$type,
           col  = object@Prefs@Graph.prefs$col,
           pch  = object@Prefs@Graph.prefs$pch,
           cex  = object@Prefs@Graph.prefs$cex,
           lty  = object@Prefs@Graph.prefs$lty,
           fill = object@Prefs@Graph.prefs$col,
           
           ## Text label setting
           ids  = NULL,
           idsmode=object@Prefs@Graph.prefs$idsmode,
           idsext =object@Prefs@Graph.prefs$idsext,
           idscex= object@Prefs@Graph.prefs$idscex,
           idsdir= object@Prefs@Graph.prefs$idsdir,
           
           ## BW settings
           bwhoriz=object@Prefs@Graph.prefs$bwhoriz,
           bwratio=object@Prefs@Graph.prefs$bwratio,
           bwvarwid=object@Prefs@Graph.prefs$bwvarwid,
           bwdotpch= object@Prefs@Graph.prefs$bwdotpch,
           bwdotcol= object@Prefs@Graph.prefs$bwdotcol,
           bwdotcex=object@Prefs@Graph.prefs$bwdotcex,
           bwreccol =object@Prefs@Graph.prefs$bwreccol,
           bwrecfill= object@Prefs@Graph.prefs$bwrecfill,
           bwreclty= object@Prefs@Graph.prefs$bwreclty,
           bwreclwd=object@Prefs@Graph.prefs$bwreclwd,
           bwumbcol =object@Prefs@Graph.prefs$bwumbcol,
           bwumblty= object@Prefs@Graph.prefs$bwumblty,
           bwumblwd= object@Prefs@Graph.prefs$bwumblwd,
           bwoutcol =object@Prefs@Graph.prefs$bwoutcol,
           bwoutcex= object@Prefs@Graph.prefs$bwoutcex,
           bwoutpch= object@Prefs@Graph.prefs$bwoutpch,

           ## Layout parameters
           grid = object@Prefs@Graph.prefs$grid,
           logy = FALSE,
           logx = FALSE,

           ## Force x variables to be continuous
           force.x.continuous = TRUE,
           
           ## bins
           binvar = NULL,
           bins   = 10,
           #xvar   = NULL,
           ...

           ) {
    #cat(x,"\n")
    #cat(str(x))
    
    
    if(!is.null(samp)) {
      data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,samp=samp)
    } else {
      data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst)
    }
    
    ## if lengths disagree, re-read x
    if (length(x) != nrow(data)) {
    #cat(length(x))
    #cat(nrow(data))
      for (i in 1:length(names(data))) {
         if (names(data)[i] == xvarnam) {
         #if (names(data)[i] == binvar) {
          x <- as.vector(as.matrix(data[i]))
        } 
      }
    #cat(length(x))
    #cat(nrow(data))
    } 
     
    #cat(str(x)) 
    
    if(force.x.continuous != FALSE) {
      if(length(unique(data[subscripts,xvarnam])) <= object@Prefs@Cat.levels) x <- as.factor(x)
    }

    ## Stuff common to both xy and bw
    if(grid != "F") {
      panel.grid(h = -1, v = -1)
    }

    y.bw <- xpose.bin(data, binvar, bins)
    bwhoriz <- bwhoriz
    ## Plot the data
    if(!is.factor(x) && !bwhoriz) {

    trellis.par.set(list(box.rectangle = list(col = bwreccol, fill = bwrecfill, lty = bwreclty, lwd = bwreclwd)))
    trellis.par.set(list(box.umbrella = list(col = bwumbcol, lty = bwumblty, lwd = bwumblwd)))
    trellis.par.set(list(box.dot = list(col = bwdotcol, cex = bwdotcex, pch = bwdotpch)))
    trellis.par.set(list(plot.symbol = list(col = bwoutcol, cex = bwoutcex, pch = bwoutpch)))
    
    try(
      if(any(is.null(groups))) {

     #cat(length(x))
     #cat(length(xpdb5@Data$TIME[xpdb5@Data$WRES!=0]))
     #cat(length(y.bw))
        panel.bwplot(x, y.bw,
                     col   =bwdotcol,
                     pch   =bwdotpch,
                     lty   =bwreclty,
                     type  =type,
                     cex   = bwdotcex,
                     varwidth = bwvarwid,
                     box.ratio = bwratio,
                     fill = bwrecfill
                     )
      } else {
        ord <- order(x)
        panel.superpose(x[ord],
                        y.bw[ord],
                        subscripts[ord],
                        col   =bwdotcol,
                        pch   =bwdotpch,
                        cex   =bwdotcex,
                        lty   =bwreclty,
                        type  =type,
                        groups=groups,
                        varwidth = bwvarwid,
                        box.ratio = bwratio,
                        fill = bwrecfill
                        )
      }     ) 
    
      if (ids) {
      ## Add id-numbers as plot symbols
      if(!any(is.null(ids))) {
        ids <- ids[subscripts]
        addid(x,y,ids=ids,
              idsmode=idsmode,
              idsext =idsext,
              idscex = idscex,
              idsdir = idsdir)
      }
      }
      
    } 

  }

