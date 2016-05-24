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

runsum <-
  function(object,
           dir="",
           modfile=paste(dir,"run",object@Runno,".mod",sep=""),
           listfile=paste(dir,"run",object@Runno,".lst",sep=""),
           main=NULL,
           subset=xsubset(object),
           show.plots=TRUE,
           txt.cex=0.7,
           txt.font=1,
           show.ids=FALSE,
           param.table=TRUE,
           txt.columns=2,
           force.wres=FALSE,
           ...)
{
  
  
  ## Read model file
  if(is.readable.file(modfile)) {
    modfile <- scan(modfile,sep="\n",what=character(),quiet=TRUE)
    mod.file.lines <- length(modfile)
  } else {
    cat(paste("model file",modfile,"not found, run summary not created!\n"))
    return()
  }

  ## Global settings concerning number of lines, number of columns etc.

                                        #txtnrow  <- 63       # Number of rows in each column

  parameter.list <- create.parameter.list(listfile)
  #attach(parameter.list,warn.conflicts=F)


  ## Set up screen
  grid.newpage()
  gr.width <- par("din")[1]
  gr.height <- par("din")[2]
  title.size <- 0.05
  graph.size <- 0.25
  graph2.size <- 0
  if (gr.width < gr.height){
    graph.size <- graph.size*gr.width/gr.height
    graph2.size <- graph2.size*gr.width/gr.height
  }
  text.size <- 1 - (graph.size + graph2.size + title.size)
  title.vp <- viewport(x=0, y=1, just=c("left","top"),
                       width=1, height=title.size,
                       name="title.vp")

  graph.1.vp <-(viewport(x=0, y=1-title.size, just=c("left","top"),
                         width=.25, height=graph.size,
                                        #layout=grid.layout(1,4),
                         name="graph.1.vp"))
  graph.2.vp <-(viewport(x=.25, y=1-title.size, just=c("left","top"),
                         width=.25, height=graph.size,
                                        #layout=grid.layout(1,4),
                         name="graph.2.vp"))
  graph.3.vp <-(viewport(x=.50, y=1-title.size, just=c("left","top"),
                         width=.25, height=graph.size,
                                        #layout=grid.layout(1,4),
                         name="graph.3.vp"))
  graph.4.vp <-(viewport(x=.75, y=1-title.size, just=c("left","top"),
                         width=.25, height=graph.size,
                                        #layout=grid.layout(1,4),
                         name="graph.4.vp"))
  graph.5.vp <-(viewport(x=0, y=1-title.size-graph.size, just=c("left","top"),
                         width=1, height=graph2.size,
                                        #layout=grid.layout(1,4),
                         name="graph.5.vp"))


  ## create text column viewports
  textColumnList <- vector("list",txt.columns) # empty list for viewports
  for(col.num in 1:txt.columns){
    txt.margin <- 0.015
    x.val <- 1/txt.columns*(col.num-1) + txt.margin
    w.val <- 1/txt.columns - (txt.margin)
    ##cat(paste(x.val,w.val,"\n"))
    textColumnList[[col.num]] <- viewport(x=x.val,
                                          y=text.size,
                                          just=c("left","top"),
                                          width=w.val,
                                          height=text.size,
                                          gp=gpar(lineheight=1.0,
                                            cex=txt.cex,font=txt.font
                                            ),
                                          name=paste("text",col.num,"vp",
                                            sep="."))
  }

  ## text.1.vp <- viewport(x=0.015, y=text.size, just=c("left","top"),
  ##                         w=0.485, h=text.size,gp=gpar(lineheight=1.0,
  ##                                                cex=txt.cex,font=2),
  ##                         name="text.1.vp")
  ##   text.2.vp <- viewport(x=0.515, y=text.size, just=c("left","top"),
  ##                         w=0.485, h=text.size,gp=gpar(lineheight=1.0,
  ##                                                cex=txt.cex,font=2),
  ##                         name="text.2.vp")

  ##   ## to look at how page is set up:
  ##   ##
  ##    grid.show.viewport(
  ##                        viewport(x=0.515, y=text.size, just=c("left","top"),
  ##                         w=0.485, h=text.size,gp=gpar(lineheight=1.0,
  ##                                                cex=txt.cex,font=2),
  ##                         name="text.2.vp")
  ##                       )


  ## add the title
  pushViewport(title.vp)

  if(!is.null(subset)) {
    maintit <- paste("Summary of run ",object@Runno,
                     ", ",subset,sep="")
  } else {
    maintit <- paste("Summary of run ",object@Runno,sep="")
  }

  title.gp=gpar(cex=1.2,fontface="bold") # title fontsize
  grid.text(maintit,gp=title.gp)

  upViewport()


  ## Add the plots
  if(show.plots){

    ## create plots
    lw <- list(left.padding = list(x = -0.05, units="snpc"))
    lw$right.padding <- list(x = -0.05,units="snpc")
    lh <- list(bottom.padding = list(x = -0.05,units="snpc"))
    lh$top.padding <- list(x = -0.05,units="snpc")


    if(show.ids) plt.id=TRUE else plt.id=FALSE

    ##grid.rect(gp=gpar(col="grey"))
    pushViewport(graph.1.vp)

    ##xplot1 <- dv.vs.pred(object,runsum=TRUE,
    xplot1 <- dv.vs.pred(object,main=NULL,xlb=NULL,ylb=NULL,
                         ##main=list("DV vs PRED",cex=0.00005),
                         aspect="fill",
                         subset=subset,
                         type="b",
                         ids=plt.id,
                         lty=8,
                         abllwd=2,
                         ##xlb=list("",cex=0.00001),ylb=list("",cex=0.00001),
                         cex=0.5,lwd=0.1,
                         scales=list(cex=0.7,tck=c(0.3),y=list(rot=90)),
                         lattice.options = list(layout.widths = lw,
                           layout.heights = lh),
                         ...)

    print(xplot1,newpage=FALSE)
    grid.text("DV vs PRED",x=0.5,y=1,just=c("center","top"),gp=gpar(cex=0.5))
    upViewport()

    pushViewport(graph.2.vp)
    xplot2 <- dv.vs.ipred(object,
                          main=NULL,xlb=NULL,ylb=NULL,
                                        #runsum=TRUE,
                          ##main=list("DV vs PRED",cex=0.00005),
                          aspect="fill",
                          subset=subset,
                          type="b",
                          ids=plt.id,
                          lty=8,
                          abllwd=2,
                          ##xlb=list("",cex=0.00001),
                          ##ylb=list("",cex=0.00001),
                          cex=0.5,lwd=0.1,
                          scales=list(cex=0.7,tck=c(0.3),
                            y=list(rot=90)),
                          lattice.options = list(layout.widths = lw,
                            layout.heights = lh),
                          ...)
    print(xplot2,newpage=FALSE)
    grid.text("DV vs IPRED",x=0.5,y=1,
              just=c("center","top"),
              gp=gpar(cex=0.5))
    upViewport()

    pushViewport(graph.3.vp)
    xplot3 <- absval.iwres.vs.ipred(object,
                                    main=NULL,xlb=NULL,ylb=NULL,
                                    ##runsum=TRUE,
                                    ##main=list("DV vs PRED",cex=0.00005),
                                    aspect="fill",
                                    subset=subset,
                                    ##type="b",
                                    ids=F,
                                    ##lty=8,
                                    ##abllwd=2,
                                    ##xlb=list("",cex=0.00001),
                                    ##ylb=list("",cex=0.00001),
                                    cex=0.5,lwd=0.1,
                                    scales=list(cex=0.7,tck=c(0.3),
                                      y=list(rot=90)),
                                    lattice.options = list(layout.widths = lw,
                                      layout.heights = lh),
                                    ...)
    print(xplot3,newpage=FALSE)
    grid.text("|IWRES| vs IPRED",x=0.5,y=1,
              just=c("center","top"),gp=gpar(cex=0.5))
    upViewport()

    pushViewport(graph.4.vp)
    use.cwres=TRUE
    if(force.wres){
      use.cwres=FALSE
    } else {
      if(is.null(check.vars(c("cwres"),object,silent=TRUE))) {
        use.cwres=FALSE
      }
    }
    if(use.cwres){
      xplot4 <- cwres.vs.idv(object,
                             main=NULL,xlb=NULL,ylb=NULL,
                             ##runsum=TRUE,
                             ##main=list("DV vs PRED",cex=0.00005),
                             aspect="fill",
                             subset=subset,
                             type="b",
                             ids=plt.id,
                             lty=8,
                             abllwd=2,
                             ##xlb=list("",cex=0.00001),
                             ##ylb=list("",cex=0.00001),
                             cex=0.5,lwd=0.1,
                             scales=list(cex=0.7,tck=c(0.3),y=list(rot=90)),
                             lattice.options = list(layout.widths = lw,
                               layout.heights = lh),
                             ...)
      res.txt <- "CWRES"

    }else{
      xplot4 <- wres.vs.idv(object,
                            main=NULL,xlb=NULL,ylb=NULL,
                            ##runsum=TRUE,
                            ##main=list("DV vs PRED",cex=0.00005),
                            aspect="fill",
                            subset=subset,
                            type="b",
                            ids=plt.id,
                            lty=8,
                            abllwd=2,
                            ##xlb=list("",cex=0.00001),
                            ##ylb=list("",cex=0.00001),
                            cex=0.5,lwd=0.1,
                            scales=list(cex=0.7,tck=c(0.3),y=list(rot=90)),
                            lattice.options = list(layout.widths = lw,
                              layout.heights = lh),
                            ...)
      res.txt <- "WRES"
    }
      print(xplot4,newpage=FALSE)
      grid.text(paste(res.txt,"vs",xvardef("idv",object)),
                x=0.5,y=1,
                just=c("center","top"),gp=gpar(cex=0.5))
      upViewport()


    ##       pushViewport(graph.5.vp)
    ##       xplot5 <- dv.preds.vs.idv(object,
    ##                             runsum=TRUE,
    ##                             ##main=list("DV vs PRED",cex=0.00005),
    ##                             aspect="fill",
    ##                             subset=subset,
    ##                             type="b",
    ##                             ids=plt.id,
    ##                             lty=8,
    ##                             abllwd=2,
    ##                             ##xlb=list("",cex=0.00001),ylb=list("",cex=0.00001),
    ##                             cex=0.5,lwd=0.1,
    ##                             #scales=list(cex=0.7,tck=c(0.3),y=list(rot=90)),
    ##                             lattice.options = list(layout.widths = lw, layout.heights = lh),
    ##                             ...)
    ##       print(xplot5,newpage=FALSE)
    ##       grid.text(paste("WRES vs",xvardef("idv",xpdb)),x=0.5,y=1,just=c("center","top"),gp=gpar(cex=0.5))
    ##       upViewport()


  } # end show plots

  ## add text
                                        #text.vp.list <- list(text.1.vp,text.2.vp)
  text.vp.list <- textColumnList
  ystart <- unit(1,"npc")
  vp.num <- 1
  space.avail <- TRUE

  ## Add the termination messages
  
  if(parameter.list$seenterm == 1 && space.avail) {

    termtxt <- parameter.list$term

    txt.marker <- add.grid.text(txt=termtxt,
                                ystart=ystart,
                                vp=text.vp.list,
                                vp.num=vp.num,
                                spaces.before=1,
                                ...)
    ystart <- txt.marker$ystart
    space.avail <- is.null(txt.marker$stop.pt)
    vp.num <- txt.marker$vp.num
  }

  ## Add objective
  if(parameter.list$seenobj == 1 && space.avail) {
    obj.txt <- paste("Objective:",parameter.list$ofv)
    txt.marker <- add.grid.text(txt=obj.txt,
                                ystart=ystart,
                                vp=text.vp.list,
                                vp.num=vp.num,
                                spaces.before=1,
                                ...)
    ystart <- txt.marker$ystart
    space.avail <- is.null(txt.marker$stop.pt)
    vp.num <- txt.marker$vp.num
  }


###############################
  ## Table of parameters and RSEs
################################

  table.txt <- list(parameter.list$parnam,format(parameter.list$parval,digits=3))
  table.col.names <- c("Par","Val")

  have.ses   <- 0
  if(parameter.list$seenseth ==1 || parameter.list$seenseom==1 || parameter.list$seensesi==1) {
    have.ses     <- 1
    table.txt <- list(parameter.list$parnam,format.default(parameter.list$parval,digits=3),parameter.list$separval)
    table.col.names <- c("Par","Val","RSE")
  }

  ##ystart <- unit(3,"lines")
  txt.marker <- add.grid.table(table.txt,
                               col.nams=table.col.names,
                               ystart=ystart,
                               vp=text.vp.list,
                               vp.num=vp.num,
                               ##center.table=TRUE,
                               ##col.optimize=FALSE,
                               ##equal.widths=TRUE,
                               ##mult.col.padding=2,
                               ...)
  ystart <- txt.marker$ystart
  space.avail <- is.null(txt.marker$stop.pt)
  vp.num <- txt.marker$vp.num





  ## Add model file
  if(space.avail) {
    txt.marker <- add.grid.text(txt=modfile,
                                ystart=ystart,
                                vp=text.vp.list,
                                vp.num=vp.num,
                                spaces.before=1,
#                                spaces.before=0,
                                ...)
    ystart <- txt.marker$ystart
    space.avail <- is.null(txt.marker$stop.pt)
    vp.num <- txt.marker$vp.num
  }

  #detach(parameter.list)
  invisible()

                                        #return()

}
