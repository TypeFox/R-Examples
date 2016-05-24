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

"xpose.multiple.plot.default" <-
  function(plotList,
           plotTitle=NULL,
           prompt=FALSE,
           new.first.window=FALSE,
           max.plots.per.page=4, # absolute max is 9
           #title.size=0.1, # title size
           title    = list(
             title.x = unit(0.5, "npc"),
             title.y = unit(0.5, "npc"),
             title.gp= gpar(cex=1.2,fontface="bold"),#,font=2),
             title.just = c("center","center")
             ),
#           title.x=unit(0.5, "npc"),  # title placement
#           title.y=unit(0.5, "npc"), # title placement
#           title.just=c("center","center"), # title placement
#           title.gp=gpar(cex=1.2,fontface="bold"), # title fontsize
           mirror=FALSE,
           ##record=TRUE,
           ##main=NULL,
           ##object,
           ##main = NULL,
           ##xlb  = NULL,
           ##ylb  = NULL,
           ##onlyfirst=TRUE,
           ##inclZeroWRES=FALSE,
           ##subset=xsubset(object),
           ## abline=c(0,1),
           ##smooth=TRUE,
           ##abllwd=2,
           bql.layout=FALSE,
           page.numbers=TRUE,
           ...) {


    ## Extract title graphical parameters
    title.x    <- title$title.x
    title.y    <- title$title.y
    title.gp   <- title$title.gp
    title.just <- title$title.just

    ## flatten plotList if we have lists of lists
    if (mirror) {
      if(length(plotList[[1]])==2 | length(plotList[[1]])==4) {
        plotList <- unlist(plotList,recursive=FALSE)
      }
    }

    ## plots per page
    absolute.max.plots.per.page = 9
    if (max.plots.per.page > absolute.max.plots.per.page) {
      max.plots.per.page = absolute.max.plots.per.page
    }

    if(bql.layout) max.plots.per.page=2

    ## split the pages and find the number of pages needed
    ## should use n2mfrow() here!
    tot.pages <- ceiling(length(plotList)/max.plots.per.page)

    if (max.plots.per.page==1) splt = c(1,1)
    if (max.plots.per.page==2) {
      if (length(plotList)==1) splt = c(1,1)
      if (length(plotList) > 1) splt = c(2,1)
    }
    if (max.plots.per.page==3 ||
        max.plots.per.page==4) {
      if (length(plotList)==1) splt = c(1,1)
      if (length(plotList)==2) splt = c(2,1)
      if (length(plotList)>2) splt = c(2,2)
    }
    if (max.plots.per.page==5 ||
        max.plots.per.page==6) {
      if (length(plotList)==1) splt = c(1,1)
      if (length(plotList)==2) splt = c(2,1)
      if (length(plotList)==3) splt = c(2,2)
      if (length(plotList)==4) splt = c(2,2)
      if (length(plotList)>4)  splt = c(3,2)
    }
    if (max.plots.per.page==7 ||
        max.plots.per.page==8 ||
        max.plots.per.page==9) {
      if (length(plotList)==1) splt = c(1,1)
      if (length(plotList)==2) splt = c(2,1)
      if (length(plotList)==3) splt = c(2,2)
      if (length(plotList)==4) splt = c(2,2)
      if (length(plotList)==5) splt = c(3,2)
      if (length(plotList)==6) splt = c(3,2)
      if (length(plotList) >6) splt = c(3,3)
    }


    if(mirror) { # beginning of Mirror stuff
      ## Decide the layout of the graphs
      if(!is.logical(mirror)) {
        if(mirror != 1 && mirror !=3) {
          cat("The mirror should either be logical, 1 or 3!\n")
          invisible()
          return()
        }
      } else {
        mirror <- 1
      }
      tot.pages <- ceiling(length(plotList)/(mirror+1))
      max.plots.per.page = mirror+1
      if(mirror==1) {
        splt <- c(1,2)
      } else {
        splt <- c(2,3)
      }
    } # end of Mirror stuff



    ## Start recording (may not work in X11)
    ##if(dev.cur()==1) {
    ##  get(getOption("device"))(record=TRUE)
    ##} else {
    ##  dev.off()
    ##  get(getOption("device"))(record=TRUE)
    ##}

    ##if ((theme=="windows") || (theme=="x11") || (theme=="pdf") || (theme=="postscript")) {
    ##  theme = theme
    ##} else {
    ## theme = "windows"
    ##}


    ## set up the title
    if (!is.null(plotTitle)){
      plot.title <- textGrob(plotTitle,
                             x=title.x,
                             y=title.y,
                             just=title.just,
                             gp=title.gp)
      plot.height <- grobHeight(plot.title)
    }

    ## Loop over the terms
    j <- 1
    page.num <- 1
    for(i in 1:length(plotList)) {
      if (j==(max.plots.per.page + 1)) {
        j <- 1
        page.num <- page.num+1
        if (prompt == TRUE) {
          cat("Next plot:  page", page.num, "of", tot.pages,
              "- Press RETURN to continue...\n",  sep=" ")
          readline()
        }
      }
      if (j==1){
        devcur <- names(dev.cur())
        if(dev.cur() == 1 | new.first.window==TRUE) { # if a device is not open
          if(tot.pages==1){
            xpose.dev.new(...)
            grid.newpage()
            #trellis.device(new=FALSE,...)#, theme = canonical.theme(theme))
            ##trellis.par.set(theme = col.whitebg())
          } else { # turn on recording if there are more than one page to print
            xpose.dev.new(record=TRUE,...) # record only pased to windows
            grid.newpage()
            #trellis.device(new=FALSE,...)#, theme = canonical.theme(theme))
            ##trellis.par.set(theme = col.whitebg())
          }
        } else { # if another graphics device is open
          seen <- 0
          if (devcur == "windows") {
            seen <- 1
            if (tot.pages==1 | i!=1){
              grid.newpage()
              #plot.new()
              #trellis.device(new=FALSE)#, theme = canonical.theme(theme))
              ##trellis.par.set(theme = col.whitebg())
            } else {
              ##options(graphics.record=FALSE)
              ##grid.newpage(recording=TRUE)
              ##dev.control("enable")
              dev.off()
              xpose.dev.new(record=TRUE,...)
              grid.newpage()
              #trellis.device(new=FALSE,...)#, theme = canonical.theme(theme))
              ##trellis.par.set(theme = col.whitebg())
            }
          }
##          if ((devcur == "x11") | (devcur == "X11") | (devcur=="quartz")) {
##             seen <- 1
##             if (tot.pages==1 | i!=1){
##               grid.newpage()
##               ##trellis.device(new=FALSE,...)#, theme = canonical.theme(theme))
##               ##trellis.par.set(theme = col.whitebg())
##             } else {
##               ##grid.newpage(recording=TRUE)
##               ##dev.control("enable")
##               ##dev.off()
##               get(getOption("device"))()
##               grid.newpage()
##               trellis.device(new=FALSE,...)#, theme = canonical.theme(theme))
##               ##trellis.par.set(theme = col.whitebg())
##             }
##           }
          if (seen!=1) {
            grid.newpage()
          }
        }

        if (is.null(plotTitle)){
          if (tot.pages>1){
            lvp <- viewport(y=0,height=unit(1, "npc") - unit(.025, "npc"),
                            just="bottom", name="lvp")
            tvp <- viewport(y=1, height=unit(.025, "npc"),
                            just="top", name="tvp")
          } else {
            lvp <- viewport(y=0,height=unit(1, "npc"),
                            just="bottom", name="lvp")
          }
        } else {
          if(length(plotList)>1 | any(class(plotList[[i]])=="grob")){
            lvp <- viewport(y=0,
                            height=unit(1, "npc") - plot.height*1.1,
                            just="bottom", name="lvp")

            ##           lvp <- viewport(y=0,height=unit(1, "npc") - unit(title.size, "npc"),
            ##                           just="bottom", name="lvp")
                                        #tvp <- viewport(y=1, height=unit(title.size, "npc"),
                                        #                just="top", name="tvp",
                                        #                gp=gpar(cex=1.2,fontface="bold")
                                        #                )

            ## find how many \n there are in the string
            ##           tvp <- viewport(y=1, height=grobHeight(grid.text(plotTitle)),
            ##                           #stringHeight(plotTitle),
            ##                           just="top", name="tvp",
            ##                           gp=gpar(cex=1.2,fontface="bold")
            ##                           )

            tvp <- viewport(y=1, height=plot.height*1.1,
                                        #stringHeight(plotTitle),
                            just="top", name="tvp"#,
                                        #gp=gpar(cex=1.2,fontface="bold")
                            )
                                        #grid.show.viewport(lvp)
                                        #browser()
                                        #pushViewport(lvp)
                                        #grid.rect()
                                        #upViewport()
                                        #pushViewport(tvp)
                                        #grid.rect()
                                        #upViewport()
            #for(jj in 1:length(plotList)){
            #    plotList[[jj]] <- update(plotList[[jj]],main$cex <- 0.5)
            #}
          } else {
            lvp <- viewport(y=0,height=unit(1, "npc"),
                            just="bottom", name="lvp")
            plotList[[i]] <- update(plotList[[i]],main=plotTitle)
            plotTitle <- NULL
          }

        }
      }
      mre=TRUE
      if (i==length(plotList)) mre=FALSE
      pushViewport(lvp)
      if(any(class(plotList[[i]])=="grob")){
        grid.draw(plotList[[i]])
      } else {
        if (mirror){
          if (j==1){
            if(mirror==1) {
              print(plotList[[i]],split=c(1,1,splt),more=mre,newpage=FALSE)
            } else {
              print(plotList[[i]],split=c(1,2,splt),more=mre,newpage=FALSE)
            }
          } else {
            if(mirror==1) {
              print(plotList[[i]],split=c(1,j,splt),more=mre,newpage=FALSE)
            } else {
              print(plotList[[i]],split=c(2,j-1,splt),more=mre,newpage=FALSE)
            }
          }
        } else {
          if(bql.layout){
            if(j==1) print(plotList[[i]],position=c(0,0.25,1,1),more=mre,newpage=FALSE)
            if(j==2) print(plotList[[i]],position=c(0,0,1,0.33),more=mre,newpage=FALSE)
          } else {
            if (j==1) print(plotList[[i]],split=c(1,1,splt),more=mre,newpage=FALSE)
            if (j==2) print(plotList[[i]],split=c(2,1,splt),more=mre,newpage=FALSE)
            if (j==3) print(plotList[[i]],split=c(1,2,splt),more=mre,newpage=FALSE)
            if (j==4) print(plotList[[i]],split=c(2,2,splt),more=mre,newpage=FALSE)
            if (j==5) print(plotList[[i]],split=c(3,1,splt),more=mre,newpage=FALSE)
            if (j==6) print(plotList[[i]],split=c(3,2,splt),more=mre,newpage=FALSE)
            if (j==7) print(plotList[[i]],split=c(1,3,splt),more=mre,newpage=FALSE)
            if (j==8) print(plotList[[i]],split=c(2,3,splt),more=mre,newpage=FALSE)
            if (j==9) print(plotList[[i]],split=c(3,3,splt),more=mre,newpage=FALSE)
          }
        }
      }
      upViewport()

      if (j==max.plots.per.page || i==length(plotList)){

        if (!is.null(plotTitle) || tot.pages >1 ){
          pushViewport(tvp)
        }


        if (!is.null(plotTitle)){

          grid.draw(plot.title)

        }

        if (tot.pages > 1){
          if(page.numbers){
            plot.page.num <- paste("page", page.num, "of", tot.pages, sep=" ")
            grid.text(plot.page.num, x=unit(.98, "npc"),
                      y=unit(.98, "npc"),
                      just=c("right","top"),
                      gp=gpar(cex=0.8))
          }
        }

        if (!is.null(plotTitle) || tot.pages >1 ){
          upViewport()
        }
      }

      j <- j+1

    }

    invisible()
  }


