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

add.grid.table <- function(txt, # list of character vectors
                                        # corresponding to
                                        # columns in table
                                        # list(parnam, parval)
                           col.nams=NULL, # collumn names c("col1","col2")
                           ystart, # y-coordinate in viewport to start
                           xstart=unit(0,"npc"), # x coordinate in viewport
                           start.pt=1, # point in lists to start
                           vp, # list of viewport names list(vp1,vp2)
                           vp.num=1, # number of viewport to begin with
                           minrow=5, # minimum rows in each column
                           cell.padding=0.5, # padding in each cell of table
                                           # in character units
                           mult.col.padding=1, # padding between
                                               # new new columns of the table
                           col.optimize=TRUE, # get maximum columns or rows in viewport
                           equal.widths=FALSE, # should all columns have equal widths
                           space.before.table=1,
                           center.table=FALSE,
                           use.rect=FALSE,
                           fill.type=NULL, # all,top,side,both,NULL
                           fill.col="grey",
                           cell.lines.lty=0,
                           ...) {


  ## to do:
  ##
  ## add grid rectangles if asked for
  ## table title

  for(k in vp.num:length(vp)){
    vp.num = k
    
    pushViewport(vp[[k]])
    wdth <- rep(0,length(txt))
    ncols <- length(txt)
    npar <- length(txt[[1]])
    for (i in 1:ncols){
      wdth[i] <- convertWidth(max(unit(rep(1, length(txt[[i]])),"strwidth",
                                       data = as.list(txt[[i]])))+unit(cell.padding,"char"),
                              "npc",valueOnly=TRUE)
      if(!is.null(col.nams)) {
        tmp <- convertWidth(max(unit(rep(1, length(col.nams[i])),"strwidth",
                                     data = as.list(col.nams[i])))
                            +unit(cell.padding,"char"),
                            "npc",valueOnly=TRUE)
        if(tmp>wdth[i]) wdth[i]=tmp
      }
    }
    if (equal.widths){
      wdth[1:length(wdth)]=max(wdth)
    }
    wdth.tot <- sum(wdth) + convertWidth(unit(mult.col.padding,"char"),"npc",valueOnly=TRUE)
    
    vp.lines.left <- convertHeight(ystart,"lines",valueOnly=TRUE)
    vp.width.left <- 1 - convertWidth(xstart,"npc",valueOnly=TRUE)
    
    ncol.grps.avail <- floor(vp.width.left/wdth.tot) # a multiple of ncols
    nrows.avail <- floor(vp.lines.left)
    xtra.space <- vp.width.left - wdth.tot*ncol.grps.avail
    
    upViewport()
    
    ## make sure we have space for table
    if (ncol.grps.avail==0){
      warning("fontsize too large for table to print")
      ret.list <- list(ystart = ystart,  # new starting point for new text
                       stop.pt = NULL, 
                       vp.num = vp.num # the viewport needed for this to work
                       )
      
      return(ret.list)
    }
    if (nrows.avail<2){
      space.avail <- TRUE
      if(!is.null(col.nams)||nrows.avail==0){
        space.avail <- FALSE
        stop.pt <- start.pt
      }
    } else {
      space.avail <- TRUE
    }
    
    ## Decide how many rows in each column of parameters
    if (col.optimize){
      if(npar <= minrow*ncol.grps.avail) {
        ncol.grps.needed <- ceiling(npar/minrow)
        if(ncol.grps.needed==1)  par.in.col <- npar
        if(ncol.grps.needed > 1) par.in.col <- minrow
      } else {
        par.in.col <- ceiling(npar/ncol.grps.avail)
        ncol.grps.needed <- ncol.grps.avail
      }
    } else { # row optimize
      ncol.grps.needed <- ceiling(npar/nrows.avail)
      if(ncol.grps.needed==1)  par.in.col <- npar
      if(ncol.grps.needed > 1) par.in.col <- nrows.avail
    }
    
    ## To make it easier to print the parameters we add empty entries
    ## to fill all columns
    max.par      <- ncol.grps.needed*par.in.col
    extra.pars   <- max.par-npar


    tmp.txt <- txt
    if(extra.pars > 0) {
      for (i in 1:ncols){
        tmp.txt[[i]][(npar+1):max.par] <- ""
      }
    }

    xpose.table <- gTree(name="xpose.table")
    for( i in 1:ncol.grps.needed) {
      range.to.print <-
        (1+(i-1)*par.in.col+(start.pt-1)):(i*par.in.col)
      
      ## Print out table
      wdth.prev.cols <- 0
      for (j in 1:ncols){
        ystart=ystart
        xstart=(i-1)*unit(wdth.tot,"npc")+unit(wdth.prev.cols,"npc")
        if (center.table){
          xstart = xstart + unit(xtra.space/2,"npc") +
            unit((ncol.grps.avail - ncol.grps.needed)*wdth.tot/2,"npc")
        }
        col.txt <- tmp.txt[[j]]
        ## add column names
        if(!is.null(col.nams)) {
          if(!all(col.txt[range.to.print]=="")){
            col.txt <- c(col.nams[j],col.txt[range.to.print])
          } else {
            col.txt <- c(" ",col.txt[range.to.print])
          }
        }
        
        ## add column text
        if(space.avail) {
          fill.type.tmp <- fill.type
          if(!is.null(fill.type)){
            if(fill.type=="top") fill.type.tmp <- "firstonly"
            if(fill.type=="both"){
              fill.type.tmp <- "firstonly"
              if(j==1) fill.type.tmp <- "all"
            }
          }
          txt.marker <- add.grid.text(txt=col.txt,
                                      ystart=ystart,
                                      xstart=xstart,
                                      start.pt=1,
                                      vp=list(vp[[k]]),
                                      vp.num=1,
                                      spaces.before=space.before.table,
                                      wdth=wdth[j],
                                      fill.type=fill.type.tmp,
                                      fill.col=fill.col,
                                      cell.lines.lty=cell.lines.lty,
                                      xpose.table=xpose.table,
                                      use.rect=use.rect,
                                      ...)
          xpose.table <- txt.marker$xpose.table
        
        }
        wdth.prev.cols <- wdth.prev.cols+wdth[j]
      }
      
    }
    if(space.avail) {
      ystart <- txt.marker$ystart
      stop.pt <- txt.marker$stop.pt
    }
    if (!is.null(stop.pt)){
      start.pt=stop.pt
      ystart=unit(1,"npc")
      xstart <- unit(0,"npc")
    } else{
      break
    }
  }
  
  
  ret.list <- list(ystart = ystart,  # new starting point for new text
                   stop.pt = stop.pt, # null if everything gets printed
                   vp.num = vp.num, # the viewport needed for this to work
                   xpose.table=xpose.table # a grob object
                   )
  
  
  return(ret.list)
  
}
  
