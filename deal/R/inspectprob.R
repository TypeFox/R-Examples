## inspectprob.R
## Author          : Claus Dethlefsen
## Created On      : Sun Feb 03 15:02:14 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Wed Jul 28 09:39:34 2004
## Update Count    : 35
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

inspectprob <-  function(nw,unitscale=20,cexscale=8,
                         arrowlength=.25,xr=c(0,350),yr=xr,...) {

    ## arguments are the same as for plot.network.
  
    par(mfrow=c(1,1))  
    plot(x=nw,unitscale=unitscale,cexscale=cexscale,arrowlength=arrowlength,xr=xr,yr=yr,...)
    title("Inspect/Change initial probability distribution")
    
    xc <- mean(xr)
    yc <- mean(yr)
    
    points(xc,yc,cex=cexscale+4,pch=5)
    text(xc,yc,"Stop")
    
    mode <- "Inspect"
    
    newnet <- nw
    quit   <- FALSE
    unit   <- 2*pi/nw$n
    where <- t(matrix(unlist(lapply(newnet$nodes,
                                    function(x) x$position)),nrow=2)) 
    where <- rbind(where,c(xc,yc))
    
    buttonx <- 20
    buttony <- 30
    where <- rbind(where,c(2*xc-buttonx,2*yc))
    where <- rbind(where,c(2*xc-buttonx,2*yc-buttony))
    
    nlist  <- names(nw$nodes)
    while(!quit) {
        
        if (mode=="Inspect") {
            bgadd <- "black"; fgadd <- "white";
            bgrem <- "white"; fgrem <- "black";
        }
        if (mode=="Change") {
            bgadd <- "white"; fgadd <- "black";
            bgrem <- "black"; fgrem <- "white";
        }
        
        symbols(2*xc-buttonx,2*yc,rectangles=matrix(c(2,1),1),add=TRUE,bg=bgadd)
        text(2*xc-buttonx,2*yc,"Inspect",col=fgadd)
        symbols(2*xc-buttonx,2*yc-buttony,rectangles=matrix(c(2,1),1),add=TRUE,bg=bgrem)
        text(2*xc-buttonx,2*yc-buttony,"Change",col=fgrem)
        
        from <- identify(where[,1],where[,2],rep("",nw$n+3),n=1)
        
        if (from==nw$n+1) break
        if (from==nw$n+2) { mode <- "Inspect"; next }
        if (from==nw$n+3) { mode <- "Change"; next }
        
        
        if (mode=="Change")
        {
          printline()
            cat(mode, "node",nlist[from],"\n")
            print(nw$nodes[[from]]$prob)
            cat("Want to change node",nlist[from],"\n")
            cat("Not yet implemented, sorry...\n")
        }
        else if(mode=="Inspect")
        {
          printline()
            cat(mode, "node",nlist[from],"\n")
            print(nw$nodes[[from]]$prob)
        }
        
        
        plot(newnet,unitscale=unitscale,cexscale=cexscale,arrowlength=arrowlength,xr=xr,yr=yr,...)
        title("Inspect/Change initial probability distribution")
        points(xc,yc,cex=cexscale+4,pch=5)
        text(xc,yc,"Stop")
        
    }
    plot(newnet,unitscale=unitscale,cexscale=cexscale,arrowlength=arrowlength,xr=xr,yr=yr,...)
    
    newnet
}

