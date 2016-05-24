## drawnetwork.R
## Author          : Claus Dethlefsen
## Created On      : Fri Nov 30 22:05:59 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Mon Jan 12 14:31:50 2004
## Update Count    : 292
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


drawnetwork <- function(nw,
                        df,
                        prior,
                        trylist=vector("list",size(nw)),
                        unitscale=20,
                        cexscale=8,
                        arrowlength=.25,
                        nocalc=FALSE,
                        yr=c(0,350),
                        xr=yr,
                        ...)
{  
    
    ## arguments are the similar as for plot.network.
    ## nocalc=T: don't calculate scores (for use with 'specifynetwork')
    
    par(mfrow=c(1,1))  
    plot(nw,unitscale=unitscale,
         cexscale=cexscale,arrowlength=arrowlength,
         showban=TRUE,xr=xr,yr=yr,...)
    
    xc <- mean(xr)
    yc <- mean(yr)
    
    points(xc,yc,cex=cexscale+4,pch=5)
    text(xc,yc,"Stop")
    
    
    mode <- "Add"
    banmode <- FALSE
    movemode <- FALSE
    if (length(nw$banlist)>0)
        banlist <- nw$banlist
    else
        banlist <- matrix(0,0,2)
    
    newnet <- nw
    quit   <- FALSE
    unit   <- 2*pi/nw$n

    nlist  <- names(nw$nodes)
    while(!quit) {

        where <- t(matrix(
                          unlist(
                                 lapply(newnet$nodes,
                                        function(x)x$position)
                                 ), nrow=2))
        buttonx <- 20
        buttony <- 30
        where <- rbind(where,c(xc,yc))
        where <- rbind(where,c(2*xc-buttonx,2*yc))
        where <- rbind(where,c(2*xc-buttonx,2*yc-buttony))
        where <- rbind(where,c(2*xc-buttonx,2*yc-2*buttony))
        where <- rbind(where,c(2*xc-buttonx,2*yc-3*buttony))
        
        
        if (mode=="Add") {
            bgadd <- "black"; fgadd <- "white";
            bgrem <- "white"; fgrem <- "black";
        }
        if (mode=="Remove") {
            bgadd <- "white"; fgadd <- "black";
            bgrem <- "black"; fgrem <- "white";
        }
        if (movemode) {
            bgmove <- "black"; fgmove <- "white";
        }
        else {
            bgmove <- "white"; fgmove <- "black"; }
        
        if (banmode) {
            bgban <- "black"; fgban <- "white";}
        else {
            bgban <- "white"; fgban <- "black"; }
        
        
        symbols(2*xc-buttonx,2*yc,
                rectangles=matrix(c(2,1),1),add=TRUE,bg=bgadd)
        text(2*xc-buttonx,2*yc,"Add",col=fgadd)
        symbols(2*xc-buttonx,2*yc-buttony,
                rectangles=matrix(c(2,1),1),add=TRUE,bg=bgrem)
        text(2*xc-buttonx,2*yc-buttony,"Remove",col=fgrem)
        
        symbols(2*xc-buttonx,2*yc-2*buttony,
                rectangles=matrix(c(2,1),1),add=TRUE,bg=bgban)
        text(2*xc-buttonx,2*yc-2*buttony,"Ban",col=fgban)
        
        symbols(2*xc-buttonx,2*yc-3*buttony,
                rectangles=matrix(c(2,1),1),add=TRUE,bg=bgmove)
        text(2*xc-buttonx,2*yc-3*buttony,"Move",col=fgmove)
        
        from <- identify(where[,1],where[,2],rep("",nw$n+5),n=1)
        
        if (from==nw$n+1) break
        if (from==nw$n+2) { mode <- "Add"; next }
        if (from==nw$n+3) { mode <- "Remove"; next }
        if (from==nw$n+4) { banmode <- !banmode;next }
        if (from==nw$n+5) { movemode <- !movemode;next }
        
        if (movemode) 
            to <- unlist(locator(1))
        else
            to <- identify(where[,1],where[,2],rep("",nw$n+5),n=1)
        
        if (to==nw$n+1) break
        if (to==nw$n+2) { mode <- "Add"; next }
        if (to==nw$n+3) { mode <- "Remove"; next }
        if (to==nw$n+4) { banmode <- !banmode;next }
        if (to==nw$n+5) { movemode <- !movemode;next }
        
        if (!movemode) {
            if (!banmode) {
                if (mode=="Add") {
                    tempnet <-
                        insert(newnet,from,to,df,prior,nocalc,
                               trylist=trylist)
                }
                else if(mode=="Remove")
                    tempnet <- remover(newnet,from,to,df,prior,nocalc,
                                       trylist=trylist)
                
                
                if (length(tempnet$nw)>0) {
                    if (!cycletest(tempnet$nw)) {
                        newnet <- tempnet
                        trylist <- newnet$trylist
                        newnet <- newnet$nw
                    }        
                    else
                        cat("Oh, no - you created a cycle. Try again\n")
                }
                else cat("something happened\n")
            }
            else {
                ##        cat("banmode is on...\n")
                if (mode=="Add") {
                    ##  cat("Trying to add",from,"->",to,"to banlist\n")
                    if (from==to) {
                        cat("Can't add the arrow:",from,"->",to,"\n")
                        next
                    }
                    else if (nw$nodes[[to]]$type=="discrete" &
                             nw$nodes[[from]]$type=="continuous")
                    {
                        cat("Arrow (",from,"->",to,") illegal\n")
                        next
                    }
                    else if (!is.na(match(from,newnet$nodes[[to]]$parents))) {
                        cat("Can't add arrow(",from,"->",to,")\n",
                            "it's already in the graph\n")
                        next
                    }
                    banlist <- rbind(banlist,c(from,to))
                }
                else if(mode=="Remove") {
                    ## cat("Trying to remove",from,"->",to,"from banlist\n")
                    if (!nrow(banlist)>0) {
                        ## cat("nothing in banlist\n")
                        next
                    }
                    idx <- (1:nrow(banlist))[banlist[,1]==from]
                    if (!length(idx)>0) {
                        ## cat("Not in banlist\n")
                        next
                    }
                    if (!is.na(match(to,banlist[idx,2]))) {
                        ## cat("removing from banlist\n")
                        banlist <- banlist[-idx[match(to,banlist[idx,2])],]
                        banlist <- matrix(banlist,ncol=2)
                        next
                    }
                    
                    ##  cat("Its not in the banlist\n")
                }
            }
        }
        else {
            ## cat("changing (",nw$nodes[[from]]$position,") to (",to,")\n")
            newnet$nodes[[from]]$position <- to
        }
        
        
        newnet$banlist <- banlist
        plot(newnet,unitscale=unitscale,cexscale=cexscale,
             arrowlength=arrowlength,showban=TRUE,xr=xr,yr=yr,...)
        points(xc,yc,cex=cexscale+4,pch=5)
        text(xc,yc,"Stop")
    }
    plot(newnet,unitscale=unitscale,
         cexscale=cexscale,arrowlength=arrowlength,
         showban=TRUE,xr=xr,yr=yr,...)

    if (!nocalc) newnet <- learn(newnet,df,prior)$nw
    
  list(nw=newnet,trylist=trylist)
}

