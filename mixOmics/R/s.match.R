# Copyright (C) 2015 
# ADE4 package
# with modification from Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# ---------------
# this function was borrowed from the ade4 package and modified for mixOmics
# --------------


s.match =function (df1xy, df2xy, xax = 1, yax = 2, pch = 20, cpoint = 1, 
                   
                   label = row.names(df1xy), clabel = 1, edge = TRUE, xlim = NULL, 
                   
                   ylim = NULL, grid = FALSE, addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
                   
                   origin = c(0, 0), sub = "", csub = 1.25, possub = "bottomleft", 
                   
                   pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE, col, lty=1) 
  
{
  
  arrow1 <- function(x0, y0, x1, y1, length = 0.1, angle = 15, lty = 1.5, col,
                     
                     edge) {
    
    d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
    
    if (d0 < 1e-07) 
      
      return(invisible())
    
    segments(x0, y0, x1, y1, lty = lty, col=col)  #trace les traits des flches  lwd ici ?
    
    h <- strheight("A", cex = par("cex"))
    
    #  if (d0 > 2 * h) {
    
    x0 <- x1 - h * (x1 - x0)/d0
    
    y0 <- y1 - h * (y1 - y0)/d0
    
    if (edge) 
      
      arrows(x0, y0, x1, y1, angle = angle, length = length, lty = lty, col=col)
    
    #  }
    
  }  #fin function arrow
  
#   df1xy <- data.frame(df1xy)
#   df2xy <- data.frame(df2xy)
  
  # this has been added in case there is a difference of sign between the two df1xy and df2xy.
  # first check the sign
  sigcor = sign(diag(cor(df1xy, df2xy)))
  # then multiply the previous df1xy by the sign from the correlation to 'swap signs'
df1xy = as.data.frame(as.matrix(df1xy) * matrix(rep(sigcor,nrow(df1xy)),ncol=ncol(df1xy),byrow=T)) 


  n <- nrow(df1xy)
  
  
  
  opar <- par(mar = par("mar"))
  
  on.exit(par(opar))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  
  coo <- scatterutil.base(dfxy = rbind.data.frame(df1xy, df2xy), 
                          
                          xax = xax, yax = yax, xlim = xlim, ylim = ylim, grid = grid, 
                          
                          addaxes = addaxes, cgrid = cgrid, include.origin = include.origin, 
                          
                          origin = origin, sub = sub, csub = csub, possub = possub, 
                          
                          pixmap = pixmap, contour = contour, area = area, add.plot = add.plot)
  
  #fleches
  
  for (i in 1:n) {
    
    arrow1(coo$x[i], coo$y[i], coo$x[i + n], coo$y[i + n], lty = lty, edge = edge, col=col[i])
    
  }
  
  if (cpoint > 0)  #met les points 
    
    points(coo$x[1:n], coo$y[1:n], pch = pch, cex = par("cex") *  cpoint, col=col)
  
  if (clabel > 0) { #label
    
    a <- coo$x[1:n] +0.01              #(coo$x[1:n] + coo$x[(n + 1):(2 * n)])/2
    
    b <- coo$y[1:n] + 0.01           #(coo$y[1:n] + coo$y[(n + 1):(2 * n)])#/2
    
    scatterutil.eti(a, b, label, clabel, boxes=FALSE)
    
  }
  
  box()
  
}

