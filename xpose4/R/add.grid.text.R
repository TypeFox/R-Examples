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

add.grid.text <- function(txt,
                          ystart, # y-coordinate in viewport to start (top)
                          xstart=unit(0,"npc"), # x coordinate in viewport (left)
                          start.pt=1, # point in text to start
                          vp, # list of viewport names list(vp1,vp2)
                          vp.num=1, # number of viewport to begin with
                          spaces.before=NULL,
                          spaces.after=NULL,
                          v.space.before=0,
                          v.space.after=0,
                          use.rect=FALSE,
                          wdth=NULL, #the width of the column of data
                          fill.type=NULL, # all, firstonly, NULL
                          fill.col="grey",
                          cell.lines.lty=0,
                          xpose.table=gTree(), # a grob object
                          ...) {

  ystart.tmp <- ystart
  if (!is.null(spaces.before)){
    tmp <- rep(" ",spaces.before)
    txt <- c(tmp,txt)
  }
  if (!is.null(spaces.after)){
    tmp <- rep(" ",spaces.after)
    txt <- c(txt,tmp)
  }
  for(j in vp.num:length(vp)){
    vp.num.final = j
    pushViewport(vp[[j]])

    ## need size of each line and draw a grid around it
          
    ## determine line hight of first row and space left 
    #line.ht <- convertHeight(unit(1,"lines"),"npc")
    #space.left <- convertHeight(ystart - line.ht,"npc",valueOnly=TRUE)
    
    ## loop through text
    stop.pt <- NULL
    #xpose.cells <- gTree(name="xpose.cells")

    for (i in start.pt:length(txt)){
      new.txt <- txt[i]
      new.txt.lines <- gregexpr("\n",new.txt)
      if(all(new.txt.lines[[1]]==-1)) {
        num.lines=1
      } else {
        num.lines=length(new.txt.lines[[1]])+1
      }
      line.ht <- convertHeight(unit(num.lines,"lines"),"npc")
      over.space <- convertHeight(unit(v.space.before,"lines"),"npc")
      under.space <- convertHeight(unit(v.space.after,"lines"),"npc")
      tot.ht <- line.ht+over.space+under.space
      space.left <- convertHeight(ystart - tot.ht,"npc",valueOnly=TRUE)
      if(space.left>=0){
        txt.width <- convertWidth(stringWidth(new.txt),"npc",valueOnly=TRUE)
        while (txt.width >= 1) {
          new.txt <- substr(new.txt,1,nchar(new.txt)-1)
          txt.width <- convertWidth(stringWidth(new.txt),"npc",valueOnly=TRUE)
        }
        
        fill.col.tmp=fill.col
        if(is.null(fill.type)){
          fill.col.tmp <- NULL
        } else {
          if(fill.type=="firstonly"){
            fill.col.tmp <- NULL
#            fill.list <- c(1:(spaces.before+1))
            fill.list <- spaces.before+1
            if(any(i==fill.list)) fill.col.tmp <- fill.col 
          }
          if(fill.type=="all"){
            fill.col.tmp <- NULL
#            fill.list <- c(1:(spaces.before+1))
            fill.list <- c(1:spaces.before)
            if(!any(i==fill.list)) fill.col.tmp <- fill.col 
          }
        }
        if(use.rect){
          xpose.rect <- grid.rect(x=xstart,y=ystart,width=wdth,
                                  height=tot.ht,
                                  just=c("left","top"),
                                  gp=gpar(fill=fill.col.tmp,
                                    lty=cell.lines.lty,
                                    ...),
                                  ...)
          xpose.table <- addGrob(gTree=xpose.table, child=xpose.rect)
          #xpose.line <- grid.lines()
          #x=c(xstart,wdth),y=ystart,default.units="npc",
                                  #just=c("left","top"),
                                  #gp=gpar(fill=fill.col.tmp,
                                  #  lty=cell.lines.lty,
                                  #  ...),
           #                       ...)
          #xpose.table <- addGrob(gTree=xpose.table, child=xpose.line)

        }
        xpose.text <- grid.text(new.txt,x=xstart,y=(ystart-over.space),just=c("left","top"),check.overlap=TRUE,...)
        xpose.table <- addGrob(gTree=xpose.table, child=xpose.text)
        ystart <- ystart-tot.ht
      } else {
        stop.pt = i
        break
      }
    }
    
    upViewport()

    
    if (!is.null(stop.pt)){
      start.pt=stop.pt
      #ystart=ystart.tmp
      ystart=unit(1,"npc")
      vp.num.final=j+1
    } else{
      break
    }
        

  }

  ret.list <- list(ystart = ystart,  # new starting point for new text
                   stop.pt = stop.pt, # null if everything gets printed
                   vp.num = vp.num.final, # the viewport needed for this to work
                   xpose.table=xpose.table # the grob object created
                   )
  
  return(ret.list)
  
}
