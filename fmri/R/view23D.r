# plots slices of a given data set ttt (spm, pvalue or data)
# it plots a coronal, a sagittal and an axial slice at the same time
#
# is called from plot.fmridata, if type=3d
# needs tkrplot
fmri.view3d <- function(ttt, sigma=NULL,type = "data", col = grey(0:255/255), ext = 1, weights =
                        c(1,1,1), scale=c(0,1), scalecol = col,
                        hrf=rep(0,100), quant =3, maxpvalue = 0.05,pos=c(-1,-1,-1)) {
  if (!require(tkrplot))
    stop("required package mytkrplot not found. Please install from cran.r-project.org")
# some basic data properties
  dt <- dim(ttt)
  zlim <- range(ttt, na.rm = TRUE)
  label <- c("x", "y", "z", "t", "signal cut-off")
# center position with Tcl objects
  if (pos[1] == -1) {
    pos <- c(round(dt[1:3])/2, 1, scale[1])
  } else {
    pos <- c(pos,1,scale[1])
  }
# needed for posv, creates a tclVar
# vector of tclVars, which contains the slider position of every (of the three) slice
  posv <- lapply(pos, helpFunc) 


# creates a slice and plots it by calling mytkrplot
  fmri.image <- function(which, factor) {
    switch(which, x = {# depending on the viweing direction
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") 
           thresh <- (as.numeric(tclvalue(posv[[5]])) - scale[1])/diff(scale)  
# plot image
        if (length(dt) == 4) {# data 
          slice <- ttt[pos[1],,,pos[4]]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[2],1:dt[3],slice, col=col, zlim=zlim)
        } else {
          slice <- ttt[pos[1],,]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[2],1:dt[3],slice, col=col, zlim=zlim)
        }
# mark position
        lines(c(pos[2],pos[2]), c(0,dt[3])+0.5, col=2)
        lines(c(0,dt[2])+0.5, c(pos[3],pos[3]), col=2)
      }
    }, y = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") 
           thresh <- (as.numeric(tclvalue(posv[[5]])) - scale[1])/diff(scale)  
# plot image
        if (length(dt) == 4) {
          slice <- ttt[,pos[2],,pos[4]]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[3],slice, col=col, zlim=zlim)
        } else {
          slice <- ttt[,pos[2],]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[3],slice, col=col, zlim=zlim)
        }
# mark position
        lines(c(pos[1],pos[1]), c(0,dt[3])+0.5, col=2)
        lines(c(0,dt[1])+0.5, c(pos[3],pos[3]), col=2)
      }
    }, z = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") 
           thresh <- (as.numeric(tclvalue(posv[[5]])) - scale[1])/diff(scale)  
# plot image
        if (length(dt) == 4) {
##          slice <- ttt[,dt[2]:1,pos[3],pos[4]]
          slice <- ttt[,,pos[3],pos[4]]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[2],slice, col=col, zlim=zlim)
        } else {
##          slice <- ttt[,dt[2]:1,pos[3]]
          slice <- ttt[,,pos[3]]
          if (type == "spm") slice[slice<thresh] <- 0          
          image(1:dt[1],1:dt[2],slice, col=col, zlim=zlim)
        }
# mark position
        lines(c(pos[1],pos[1]), c(0,dt[2])+0.5, col=2)
        lines(c(0,dt[1])+0.5, c(dt[2]-pos[2]+1,dt[2]-pos[2]+1), col=2)
      }
    })      
# create the Tk-widget
    mytkrplot(tt3, f, hscale=ext, vscale=factor*ext,typeV=3)
  }
# creates a slider (position of the slice in its viewing direction) for a given slice (indicated by the number i)
# defines the reaction on a slider move:
# each slice gets replotted, position of slice i gets updated     
  fmri.slider <- function(i) {
    f <- function(...) {# f() gets called, if a slider is moved
      current <- as.numeric(tclvalue(posv[[i]]))
      if (current != pos[i]) {
        pos[i] <<- current
        tkrplot::tkrreplot(img[[1]])
        tkrplot::tkrreplot(img[[2]])
        tkrplot::tkrreplot(img[[3]])
        if (i == 4) tkrplot::tkrreplot(img[[4]])
        tkconfigure(label2, text=pos[i])
      }
    }
# slider created and added to toplevel tt3
    fr <- tkframe(tt3,bg="#AACCDB")
    s <- tkscale(fr, command=f, from=1, to=dt[i], resolution=1, 
                 variable=posv[[i]], showvalue=FALSE, orient="horiz",bg="#BBDDEC")
    label1 <- tklabel(fr, text=label[i],bg="#BBDDEC")
    label2 <- tklabel(fr, text=pos[i],bg="#BBDDEC")
    tkgrid(label1, s, label2)
    fr
  }

# creates a slider for the positioning of the threshold
# if the value of the threshold is changed the slices will be replotted
  fmri.threshold <- function(i) {
    f <- function(...) {
      current <- as.numeric(tclvalue(posv[[i]]))
      if (current != pos[i]) {
        pos[i] <<- current
        tkrplot::tkrreplot(img[[1]])
        tkrplot::tkrreplot(img[[2]])
        tkrplot::tkrreplot(img[[3]])
        tkconfigure(label2, text=pos[i])
      }
    }
    fr <- tkframe(tt3,bg="#BBDDEC")
    s <- tkscale(fr, command=f, from=scale[1], to=scale[2], resolution=diff(scale)/100, 
                 variable=posv[[i]], showvalue=FALSE, orient="horiz",bg="#BBDDEC")
    label1 <- tklabel(fr, text=label[i],bg="#BBDDEC")
    label2 <- tklabel(fr, text=pos[i],bg="#BBDDEC")
    tkgrid(label1, s, label2)
    fr
  }

# adds a scale for each type of data (spm, pvalue, data)
  fmri.scale <- function(which,scale=scale, scalecol=scalecol) {
     switch(which, "data" = {
      f <- function() {
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        layout(matrix(1:2,2,1,byrow=TRUE),widths=c(200),heights=c(160,40))
        on.exit(par(oldpar))
# plot timeseries
        plot(ttt[pos[1],pos[2],pos[3],], xlab="Scan", ylab="BOLD signal")
# mark scan number position
        lines(c(pos[4], pos[4]),range(ttt[pos[1],pos[2],pos[3],]),col=2)
# draw scale 
        image(seq(scale[1],scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
              matrix(rep(seq(scale[1],scale[2],length=100),10),100,10),
              yaxt="n",xlab="", ylab="",zlim=scale, col=scalecol)
      }
# create the Tk-widget
      mytkrplot(tt3, f, hscale=ext, vscale=ext,typeV=3)
    }, "spm" = {
      f <- function() {
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        layout(matrix(1:2,2,1,byrow=TRUE),widths=c(200),heights=c(160,40))
        on.exit(par(oldpar))
# draw something
        if (!is.null(sigma)) {
          value <- scale[1]+ttt[pos[1],pos[2],pos[3]]*diff(scale)          
          plot(c(1,length(hrf)),range(c(value*hrf,(value-3*sigma[pos[1],pos[2],
               pos[3]])*hrf,(value+3*sigma[pos[1],pos[2],pos[3]])*hrf)),type="n",
               xlab="Scan",ylab="Parameter estimate")
          xx <- c(1:length(hrf),length(hrf):1)
          yy <- c((value-quant*sigma[pos[1],pos[2],pos[3]])*hrf,
               rev((value+quant*sigma[pos[1],pos[2],pos[3]])*hrf))
          polygon(xx,yy,col="gray",lty=1)
          lines(value*hrf)
          lines(c(1,length(hrf)),c(0,0))
#          text(0.1,0.5,paste("Parameter:",signif(ttt[pos[1],pos[2],pos[3]],3)),pos=4,cex=1.5)
        } else {
          value <- scale[1]+ttt[pos[1],pos[2],pos[3]]*diff(scale)          
          plot(c(0,1),c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
          text(0.1,0.5,paste("t-value:",signif(value,3)),pos=4,cex=1.5)
        }
# draw scale
        image(seq(scale[1],scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
              matrix(rep(seq(scale[1],scale[2],length=100),10),100,10),
              yaxt="n",xlab="", ylab="", zlim=scale, col=scalecol)
        lines(c(value,value),scale,col="white")
      }
# create the Tk-widget
      mytkrplot(tt3, f, hscale=ext, vscale=ext,typeV=3)
    }, "pvalue" = {
      f <- function() {
        if (ttt[pos[1],pos[2],pos[3]] <= 0.5) {
          value <- -log(maxpvalue)
        } else {
          value <- scale[1]+2*(ttt[pos[1],pos[2],pos[3]]-0.5)*diff(scale)          
        }
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        layout(matrix(1:2,2,1,byrow=TRUE),widths=c(200),heights=c(160,40))
        on.exit(par(oldpar))
# draw something
        plot(c(0,1),c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",type="n",
             bty="n",bg="#AACCDB")
        if (value == -log(maxpvalue)) {
          text(0.2,0.5,paste("p-value: <",signif(exp(-value),3)),
                pos=4,cex=1.5,bg="#AACCDB")
        } else if (value == scale[2]) {
          text(0.2,0.5,paste("p-value: <",signif(exp(-value),3)),
                pos=4,cex=1.5,bg="#AACCDB")
        } else {
          text(0.2,0.5,paste("p-value:",signif(exp(-value),3)),
                pos=4,cex=1.5,bg="#AACCDB")
        }
# draw scale
#
#   test if there are any significant voxel
#
        if(scale[2] > -log(maxpvalue)){
          image(seq(-log(maxpvalue),scale[2],length=100),
                seq(scale[1],scale[2],length=10)/10,
                matrix(rep(seq(-log(maxpvalue),scale[2],length=100),10),100,10),
                yaxt="n",xaxt="n",xlab="", ylab="",zlim=c(-log(maxpvalue),scale[2]),
                col=scalecol,bg="#AACCDB")
          lines(c(value,value),scale,col=1)
          lines(c(-log(0.01),-log(0.01)),scale,col=2)
          text(-log(0.01),scale[1]+0.01*diff(scale),pos=4,"1e-2")
          lines(c(-log(0.001),-log(0.001)),scale,col=2)
          text(-log(0.001),scale[1]+0.01*diff(scale),pos=4,"1e-3")
          lines(c(-log(0.0001),-log(0.0001)),scale,col=2)
          text(-log(0.0001),scale[1]+0.01*diff(scale),pos=4,"1e-4")
          lines(c(-log(0.00001),-log(0.00001)),scale,col=2)
          text(-log(0.00001),scale[1]+0.01*diff(scale),pos=4,"1e-5")
          lines(c(-log(0.000001),-log(0.000001)),scale,col=2)
          text(-log(0.000001),scale[1]+0.01*diff(scale),pos=4,"1e-6")
          lines(c(-log(0.0000001),-log(0.0000001)),scale,col=2)
          text(-log(0.0000001),scale[1]+0.01*diff(scale),pos=4,"1e-7")
          lines(c(-log(0.00000001),-log(0.00000001)),scale,col=2)
          text(-log(0.00000001),scale[1]+0.01*diff(scale),pos=4,"1e-8")
          lines(c(-log(0.000000001),-log(0.000000001)),scale,col=2)
          text(-log(0.000000001),scale[1]+0.01*diff(scale),pos=4,"1e-9")
        } else {
           cat("No significant voxel\n")
        }
      }
# create the Tk-widget
      mytkrplot(tt3, f, hscale=ext, vscale=ext,typeV=3)
     })
  }
  
# create window
  tt3 <- tktoplevel(bg="#AACCDB")
  tkwm.title(tt3, "FMRI - View Results 3d")

# create slider and images
  if (type == "data") {
    s <- lapply(1:4, fmri.slider)
  } else if (type == "spm") {
    s <- c(lapply(1:3, fmri.slider),fmri.threshold(5))
  } else {
    s <- lapply(1:3, fmri.slider)
  }
  
  img <- list(fmri.image("x",dt[3]/dt[1]*weights[3]),
              fmri.image("y",dt[3]/dt[2]*weights[3]),
              fmri.image("z",1),
              fmri.scale(type,scale,scalecol))

# place the images and scales
  tkgrid(img[[2]], img[[1]])
  tkgrid(s[[2]], s[[1]])
  tkgrid(img[[3]], img[[4]])
  if (type == "data") {
    tkgrid(s[[3]], s[[4]])
  } else if (type == "spm") {
    tkgrid(s[[3]], s[[4]])
  } else {
    tkgrid(s[[3]])
  }
  
# return the window object to the master
  tt3
}


# creates a window, in which you can view your results (spm, pvalue or data)
# it creates an image of every possible slice
# the choice which slices are shown can be made by the user
# more on the ths will be described below 
# 
# includes a method to save your results via adimpro
#
# is called from plot.fmridata, if type=2d
# needs tkrplot
fmri.view2d <- function(ttt, sigma=NULL,type = "data", col = grey(0:255/255), ext = 1, 
    weights = c(1,1,1), scale=c(0,1), scalecol = col, hrf=rep(0,100), quant =3, 
    maxpvalue = 0.05,posNew=c(-1,-1,-1), localx, inputStuff) {
# check whether Tk/Tcl environment is present
  if (!require(tkrplot))
    stop("required package mytkrplot not found. Please install from cran.r-project.org")

# some basic data properties
  dt <- dim(ttt)
  zlim <- range(ttt, na.rm = TRUE)
  label <- c("x", "y", "z", "t", "signal cut-off")
  sliceslist <- c()# the vector which later contains all slices
  hscaleList <- c()# the vector which later contains the horziontal scale of every slice
  vscaleList <- c()# the vector which later contains the vertical scale of every slice 
  label1 <- list()# the list which later contains the dimension of every slice
  label2 <- list()
  # the list which later contains the position of every slice in its dimension (changes if the slider is moved)
  label1Th <- list()
  # lists, which all later contain a single element (the position of the threshold/time or the string label[5]/label[4]) 
  label2Th <- list() 
  label1Ti <- list()
  label2Ti <- list()
# tclVars which later contain the postion of every slice resp. the time or threshold  
  posThreshold <- 1
  posTime <- 1
  posvNew <- lapply(posNew, helpFunc) 
  posvThreshold <- lapply(posThreshold,helpFunc)# only needed if type=="spm"
  posvTime <- lapply(posTime,helpFunc)# only needed if type=="data"  

# creates a slice and plots it by calling mytkrplot  
  fmri.image <- function(which,hscale,vscale, slicenr=-1,toplevel,overview=FALSE) {
    switch(which, x = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") 
           thresh <- (as.numeric(tclvalue(posvThreshold[[1]])) - scale[1])/diff(scale)  
# plot image
        if (length(dt) == 4) {
          if (!overview) slice <- ttt[posNew[slicenr],,,posTime]
          else slice <- ttt[slicenr,,,posTime]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[2],1:dt[3],slice, col=col, zlim=zlim)
        } else {# ERSTMAL NUR FUER dim(SPM)=3 bearbeitet, rest ToDo
          if (!overview) slice <- ttt[posNew[slicenr],,]
          else slice <- ttt[slicenr,,]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[2],1:dt[3],slice, col=col, zlim=zlim)
        }
      }
    }, y = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") 
           thresh <- (as.numeric(tclvalue(posvThreshold[[1]])) - scale[1])/diff(scale)  
# plot image
        if (length(dt) == 4) {
          if (!overview) slice <- ttt[,posNew[slicenr+dt[1]],,posTime]
          else slice <- ttt[,slicenr,,posTime]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[3],slice, col=col, zlim=zlim)
        } else {
          if (!overview) slice <- ttt[,posNew[slicenr+dt[1]],]
          else slice <- ttt[,slicenr,]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[3],slice, col=col, zlim=zlim)
        }
      }
    }, z = {
      f <- function() {
      oldpar <- par(mar=c(0,0,0,0))
      on.exit(par(oldpar))
      if (type == "spm") 
         thresh <- (as.numeric(tclvalue(posvThreshold[[1]])) - scale[1])/diff(scale)
# plot image
      if (length(dt) == 4) {
##        if (!overview) slice <- ttt[,dt[2]:1,posNew[slicenr+dt[1]+dt[2]],posTime]
        if (!overview) slice <- ttt[,,posNew[slicenr+dt[1]+dt[2]],posTime]
##        else slice <- ttt[,dt[2]:1,slicenr,posTime]
        else slice <- ttt[,,slicenr,posTime]
        if (type == "spm") slice[slice<thresh] <- 0
        image(1:dt[1],1:dt[2],slice, col=col, zlim=zlim)
      } else {
##        if (!overview) slice <- ttt[,dt[2]:1,posNew[slicenr+dt[1]+dt[2]]]
        if (!overview) slice <- ttt[,,posNew[slicenr+dt[1]+dt[2]]]
##        else slice <- ttt[,dt[2]:1,slicenr]
        else slice <- ttt[,,slicenr]
        if (type == "spm") slice[slice<thresh] <- 0          
        image(1:dt[1],1:dt[2],slice, col=col, zlim=zlim)
      }
    }
  })      
# create the Tk-widget
    mytkrplot(toplevel, f, hscale=hscale, vscale=vscale,slicenr,which,
              frameSlicesList[[1]],2,overview)        
  }
# creates a slider (position of the slice in its viewing direction) for a given slice (indicated by the number i)
# defines the reaction on a slider move:
# each slice gets replotted, position of slice i gets updated   
  fmri.slider <- function(i) {
    if (i==-1){
      f <- function(...) {
        current <- as.numeric(tclvalue(posvTime[[1]]))
        posTimeHelp[currPage] <<- current  
        if (current != posTime) {
          posTime <<- current
          for (j in ((currPage-1)*nrslicesVec[1]+1):
                    ((currPage-1)*nrslicesVec[1]+nrslicesVec[currPage])){
            tkrreplot(img[[sliceslist[j]]],hscale=hscaleList[sliceslist[j]],
                                           vscale=vscaleList[sliceslist[j]])    
            tkconfigure(label2[[sliceslist[j]]], text=posNew[sliceslist[j]])
          }
          tkconfigure(label2Ti[[1]], text=posTime)
        }
      }    
      fr <- tkframe(tt,bg="#BBDDEC")
      s <- tkscale(fr, command=f, from=1, to=dt[4], resolution=1, 
             variable=posvTime[[1]], showvalue=FALSE, orient="horiz",bg="#BBDDEC")
      label1Ti[[1]] <<- tklabel(fr, text=label[4],bg="#BBDDEC")
      label2Ti[[1]] <<- tklabel(fr, text=posTime,bg="#BBDDEC")
      tkgrid(label1Ti[[1]], s, label2Ti[[1]],pady=5)
      fr
    } else {  
      islice = i
      if (i<=dt[1]) i=1 else if (i<=dt[1]+dt[2]) i=2 else i=3
      f <- function(...) {
        current <- as.numeric(tclvalue(posvNew[[islice]]))
        if (current != posNew[islice]) {
          posNew[islice] <<- current
          tkrreplot(img[[islice]],hscale=hscaleList[islice],vscale=vscaleList[islice])
          tkconfigure(label2[[islice]], text=posNew[islice])
        }
      }    
        fr <- tkframe(frameSlicesList[[1]],bg="#BBDDEC")     
        s <- tkscale(fr, command=f, from=1, to=dt[i], resolution=1, 
            variable=posvNew[[islice]], showvalue=FALSE, orient="horiz",bg="#BBDDEC")
        label1[[islice]] <<- tklabel(fr, text=label[i],bg="#BBDDEC")
        label2[[islice]] <<- tklabel(fr, text=posNew[islice],bg="#BBDDEC")
        tkgrid(label1[[islice]], s, label2[[islice]])
        fr
    }
  }
# creates a slider for the positioning of the threshold
# if the value of the threshold is changed the slices will be replotted
  fmri.threshold <- function() {
    f <- function(...) {
      current <- as.numeric(tclvalue(posvThreshold[[1]]))
      posThresHelp[currPage] <<- current  
      if (current != posThreshold) {
        posThreshold <<- current
        for (i in ((currPage-1)*nrslicesVec[1]+1):
                  ((currPage-1)*nrslicesVec[1]+nrslicesVec[currPage])){
          tkrreplot(img[[sliceslist[i]]],hscale=hscaleList[sliceslist[i]],
                                         vscale=vscaleList[sliceslist[i]])    
          tkconfigure(label2[[sliceslist[i]]], text=posNew[sliceslist[i]])
        }
          tkconfigure(label2Th[[1]], text=posThreshold)
      }
    }
      fr <- tkframe(tt,bg="#BBDDEC")
      s <- tkscale(fr, command=f, from=scale[1], to=scale[2], resolution=diff(scale)/100, 
            variable=posvThreshold[[1]], showvalue=FALSE, orient="horiz",bg="#BBDDEC")
      label1Th[[1]] <<- tklabel(fr, text=label[5],bg="#BBDDEC")
      label2Th[[1]] <<- tklabel(fr, text=posThreshold,bg="#BBDDEC")
      tkgrid(label1Th[[1]], s, label2Th[[1]],pady=5)
      fr
  }
# adds a scale for each type of data (spm, pvalue, data)
  fmri.scale <- function(which,scale=scale, scalecol=scalecol) {
    switch(which, "data" = {
      f <- function() {
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        on.exit(par(oldpar))
# draw scale
        par(mai=c(0,0,0,0))
        image(seq(scale[1],scale[2],length=200),seq(scale[1],scale[2],length=20)/20,
              matrix(rep(seq(scale[1],scale[2],length=200),20),200,20),
              yaxt="n",xlab="", ylab="",zlim=scale, col=scalecol)
        lines(c(scale[1]+0.25*scale[2],scale[1]+0.25*scale[2]),scale,col=2)
        text(scale[1]+0.25*scale[2],scale[1]+0.01*diff(scale),
             pos=4,0.01*round(100*(scale[1]+0.25*scale[2])),col=2)
        lines(c(scale[1]+0.5*scale[2],scale[1]+0.5*scale[2]),scale,col=2)
        text(scale[1]+0.5*scale[2],scale[1]+0.01*diff(scale),
             pos=4,0.01*round(100*(scale[1]+0.5*scale[2])),col=2)
        lines(c(scale[1]+0.75*scale[2],scale[1]+0.75*scale[2]),scale,col=2)
        text(scale[1]+0.75*scale[2],scale[1]+0.01*diff(scale),
             pos=4,0.01*round(100*(scale[1]+0.75*scale[2])),col=2)
      }
# create the Tk-widget
      mytkrplot(tt, f, hscale=0.8, vscale=0.05,typeV=3)
    }, "spm" = {
        f <- function() {
          oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
          on.exit(par(oldpar))
# draw scale
          par(mai=c(0,0,0,0))
          image(seq(scale[1],scale[2],length=200),seq(scale[1],scale[2],length=20)/20,
                matrix(rep(seq(scale[1],scale[2],length=200),20),200,20),
                yaxt="n",xlab="", ylab="", zlim=scale, col=scalecol)
          lines(c(0,0),scale,col=2)
          text(0,0,pos=4,col=2,0)#scale[1]+0.25*scale[2],col=2)
          lines(c(0.5*scale[2],0.5*scale[2]),scale,col=2)
          text(0.5*scale[2],0,pos=4,0.01*round(100*0.5*scale[2]),col=2)
          lines(c(0.5*scale[1],0.5*scale[1]),scale,col=2)
          text(0.5*scale[1],0,pos=4,0.01*round(100*0.5*scale[1]),col=2)
        }
# create the Tk-widget
        mytkrplot(tt, f, hscale=0.8, vscale=0.05,typeV=3)
    }, "pvalue" = {
        f <- function() {
          value <- -log(maxpvalue)
          oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
          on.exit(par(oldpar))
# draw scale
          par(mai=c(0,0,0,0))
          if(-log(maxpvalue) < scale[2]){
            image(seq(-log(maxpvalue),scale[2],length=200),seq(scale[1],scale[2],
                  length=20)/20,matrix(rep(seq(-log(maxpvalue),scale[2],length=200),20),
                  200,20),yaxt="n",xaxt="n",xlab="", ylab="",
                  zlim=c(-log(maxpvalue),scale[2]), col=scalecol)
            lines(c(value,value),scale,col=1)
            i <- trunc(-log(maxpvalue)/log(10))+1
            while (-log(1/10^i)<0.95*scale[2]){
              lines(c(-log(1/10^i),-log(1/10^i)),scale,col=2)
              text(-log(1/10^i),0.01*diff(scale),pos=4,paste("1e-",i,sep=""))  
              i = i + 1
            }
          } else {
            cat("No significant voxel\n")
          }
        }
# create the Tk-widget
          mytkrplot(tt, f, hscale=0.8, vscale=0.05,typeV=3)
      })
    }
# create toplevel window
    tt <- tktoplevel(bg="#AACCDB")
    tkwm.title(tt, "FMRI - View Results 2d")   
# frame in which all slices will be created
    frameSlicesList <- list()       
    frameSlicesList[[1]] <- tkframe(tt,bg="#AACCDB")
# slices which will be shown at the start  
    if (dt[3]>4){
      nrslices = 4  
      for (i in 1:4){  
        sliceslist[i] = dt[1]+dt[2]+i*round(0.2*dt[3])
      }
    } else {
      nrslices = 1
      sliceslist[1] = dt[1] + dt[2] + 1 
    }
# determine scaling factors tor the start slices
    hsc = (0.9*as.numeric(tkwinfo("screenwidth", tt)))/(480*2)
    if (type=="pvalue"){
      vsc = (0.8*as.numeric(tkwinfo("screenheight", tt))-38-2*21)/(480*2)
    } else {
      vsc = (0.8*as.numeric(tkwinfo("screenheight", tt))-38-36-2*21)/(480*2)
    }
    if (hsc>=1.2*vsc) hsc <- 1.2*vsc# keep nearly quadratic
    if (vsc>=1.2*hsc) vsc <- 1.2*hsc  
# create images and sliders
    img <- list()
    print("Creating sagittal slices ...")  
    for (i in 1:dt[1]){
      img[[i]] = fmri.image("x",1,2*dt[3]/dt[1]*weights[3],i,toplevel=tt)
      hscaleList[i] = 1
      vscaleList[i] = 2*dt[3]/dt[1]*weights[3]  
    }
    print("Creating axial slices ...")    
    for (i in 1:dt[2]){
      img[[i+dt[1]]] = fmri.image("y",1,2*dt[3]/dt[2]*weights[3],i,toplevel=tt)
      hscaleList[i+dt[1]] = 1
      vscaleList[i+dt[1]] = 2*dt[3]/dt[2]*weights[3]  
    }
    print("Creating coronal slices ...")  
    for (i in 1:dt[3]){
      img[[i+dt[1]+dt[2]]] = fmri.image("z",hsc,vsc,i,toplevel=tt)
      hscaleList[i+dt[1]+dt[2]] = hsc
      vscaleList[i+dt[1]+dt[2]] = vsc  
  }
    s <- lapply(1:(dt[1]+dt[2]+dt[3]), fmri.slider)
# applying fmri.slider to every element of s 
# place the images and scales (of start situation (2x2 slices))
    for (i in 1:round(nrslices/2)){
      if (i == round(nrslices/2) && round(nrslices/2) != (nrslices/2)) {
        tkgrid(img[[sliceslist[2*i-1]]])
        tkgrid(s[[sliceslist[2*i-1]]])
      } else {
        tkgrid(img[[sliceslist[2*i-1]]],img[[sliceslist[2*i]]])
        tkgrid(s[[sliceslist[2*i-1]]],s[[sliceslist[2*i]]])
      }
    }
    tkgrid(frameSlicesList[[1]],padx=5,pady=5)
    scaleL <- list(fmri.scale(type,scale,scalecol))
    tkgrid(scaleL[[1]],pady=5)
    scaleHeight <- (2*5 + as.numeric(tkwinfo("height",scaleL[[1]])))# height of the scale
    rbValue <- tclVar(3)# current viewing direction
    t2 <- list()
# frames for the buttons and the scale
    frame0 <- tkframe(tt,bg="#AACCDB")
    frame1 <- tkframe(tt,bg="#AACCDB")
    frame2 <- tkframe(tt,bg="#AACCDB")
# viewAll is called by the button "View all slices"
# it first sets the number of slices to all and the calls the changeFunction
    viewAll <- function(){
      tclvalue(nrSlices) <<- dt[as.numeric(tclvalue(rbValue))]
      changeFunction()
    }  
# most important function of view2d
# it calculates the scaling factors for the slices of each page
    changeFunction <- function(){
      currPage <<- 1
      oldPage <<- 1
      view = as.numeric(tclvalue(rbValue))
      if (min(as.integer(tclvalue(nrSlicespp)),dt[view])>50) tclvalue(cbVar) <<- 1
      if (as.integer(tclvalue(nrSlices)) < 1) tclvalue(nrSlices) <<- 1
      if (as.integer(tclvalue(nrSlicespp)) < 1) tclvalue(nrSlicespp) <<- 1      
      cbVal <- ((as.integer(tclvalue(cbVar))+1)%%2)
      screenwidth <<- max(as.numeric(tkwinfo("screenwidth", tt)),screenwidth)
      screenheight <<- max(as.numeric(tkwinfo("screenheight", tt)),screenheight)
      sliderheight <<- max(as.numeric(tkwinfo("height",s[[sliceslist[1]]])),sliderheight)
      sliderwidth <<- max(as.numeric(tkwinfo("width",s[[sliceslist[1]]])),sliderwidth) 
      if (type=="spm") sliderBottomHeight <<-
         max((as.numeric(tkwinfo("height",thres[[1]])) + 5),sliderBottomHeight)
      if (type=="data") sliderBottomHeight <<-
         max((as.numeric(tkwinfo("height",time[[1]])) + 5),sliderBottomHeight)
      scaleHeight <<- max((2*5 + as.numeric(tkwinfo("height",scaleL[[1]]))),scaleHeight)
# calculates the optimal number of rows and colums for a given number of slices (nrSlicesTmp)
# this methods is especially senseful, if the aspect ratio is kept !!
# returns the optimal number of rows (the optimal number of colums is determined by this)
      findOptimum <- function(nrSlicesTmp){
# determine vector of the non chosen viewing directions      
        helpVec = c(-1,-1)    
        index = 1      
        for (i in 1:3){
          if (i!=view)  { 
            helpVec[index] = i
            index = index + 1  
          }
        }
# calculate for all possible (senseful) row-col combinations the scaling factor and by this find combination where most
# area is covered (not always the most, but nearly)
        matSize <- ceiling(3*sqrt(nrSlicesTmp))
        valueMax = 0
        indexMaxRow = 1  
        for (colsTmp in 1:matSize){ 
          for (rowsTmp in 1:matSize){
            if ((colsTmp-1)*rowsTmp < nrSlicesTmp && colsTmp*rowsTmp >= nrSlicesTmp){
              hscaleTmp = (0.9*screenwidth)/(480*colsTmp)
              vscaleTmp = (0.8*screenheight-sliderBottomHeight-scaleHeight-
                               cbVal*rowsTmp*sliderheight+rowsTmp*3)/(480*rowsTmp)
              if (as.integer(tclvalue(ksVar))==1) {# keep aspect ratio          
                if (dt[helpVec[1]]>dt[helpVec[2]]){
                  if (hscaleTmp*(dt[helpVec[2]]/dt[helpVec[1]]) > vscaleTmp){
                    hscaleTmp <- (dt[helpVec[1]]/dt[helpVec[2]])*vscaleTmp
                  } else vscaleTmp = hscaleTmp*(dt[helpVec[2]]/dt[helpVec[1]]) 
                } else {
                  if (hscaleTmp*(dt[helpVec[1]]/dt[helpVec[2]]) > vscaleTmp){
                    hscaleTmp <- (dt[helpVec[2]]/dt[helpVec[1]])*vscaleTmp
                  } else vscaleTmp = hscaleTmp*(dt[helpVec[1]]/dt[helpVec[2]]) 
                }            
              } else {# nearly quadratic
                if (hscaleTmp>=1.2*vscaleTmp) hscaleTmp <- 1.2*vscaleTmp
                if (vscaleTmp>=1.2*hscaleTmp) vscaleTmp <- 1.2*hscaleTmp
              }  
              val1 <- 480*hscaleTmp*colsTmp
              val2 <- 480*vscaleTmp*rowsTmp+sliderBottomHeight+scaleHeight+cbVal*
                            rowsTmp*sliderheight+rowsTmp*3
              quot1 = val1/(screenwidth*0.9)  
              quot2 = val2/(screenheight*0.8)
              sum <- quot2+quot1
              if (sum > valueMax){
                if (as.integer(tclvalue(cbVar))==1 ||
                    (sliderwidth+20)*colsTmp<screenwidth){
                  valueMax = sum
                  indexMaxRow = rowsTmp
                }
              }                  
            }        
          }      
        }
        return (indexMaxRow)# returns the optimal number of rows
      }
      for (i in 1:nrslices){
        tkgrid.forget(img[[sliceslist[i]]])
        tkgrid.forget(s[[sliceslist[i]]])
      }
      nrslices <<- as.integer(tclvalue(nrSlices))
      nrslicespp <- as.integer(tclvalue(nrSlicespp))
      if (nrslicespp > nrslices) nrslicespp <- nrslices
# chosen nr of slices per page unnecessary big
      if (nrslices>dt[view]){# chosen nr of slices unnecessary bi
        nrslices <<- dt[view]
        print(append("Maximum number of slices:",dt[view]),quote=FALSE)    
      }
      nrpages <<- ceiling(nrslices/nrslicespp)# number of pages
      for (i in 1:nrpages){
        if (i!=nrpages) nrslicesVec[i] <<- nrslicespp
        else nrslicesVec[i] <<- nrslices-(nrpages-1)*nrslicespp     
      }
# read all buttons + arrow-buttons
      tkgrid.forget(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
                    nrSlicesppEntry,changeButton,viewAllButton,nextPageButton)
      if (nrpages>1) tkgrid(nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,nrSlicesppEntry,
                    changeButton,viewAllButton,nextPageButton,padx=10,pady=5)
      else tkgrid(nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,nrSlicesppEntry,
                    changeButton,viewAllButton,padx=10,pady=5)
# determine all chosen slices (number of them is nrslices)
      doubles = 0  
      for (i in 1:nrslices){  
        if (view==1) sliceslist[i] <<- round(i*dt[view]/(nrslices+1))
        if (view==2) sliceslist[i] <<- dt[1] + round(i*dt[view]/(nrslices+1))
        if (view==3) sliceslist[i] <<- dt[1] + dt[2] + round(i*dt[view]/(nrslices+1))
        if (i>1 && sliceslist[i-1]==sliceslist[i]) doubles = 1
      }
      if (doubles==1){
        for (i in 1:nrslices){  
          if (view==1) sliceslist[i] <<- round(i*dt[view]/(nrslices))
          if (view==2) sliceslist[i] <<- dt[1] + round(i*dt[view]/(nrslices))
          if (view==3) sliceslist[i] <<- dt[1] + dt[2] + round(i*dt[view]/(nrslices))
        }    
      }
      if (nrslices < length(sliceslist)){
        for (i in (nrslices+1):length(sliceslist)){
        sliceslist[i] <<- -1
      }
    }  
# determine the optimal number of rows and cols for the first (and every further) and the last page
    rows <- findOptimum(nrslicesVec[1])
    cols <- ceiling(nrslicesVec[1]/rows)  
    rowshelp <- findOptimum(nrslicesVec[nrpages])
    colshelp <- ceiling(nrslicesVec[nrpages]/rowshelp)  
# save choice of rows and cols for every page
    for (i in 1:(nrpages-1)){
      colsVec[i] <<- cols
      rowsVec[i] <<- rows
    } 
    colsVec[nrpages] <<- colshelp
    rowsVec[nrpages] <<- rowshelp
# calculate standard scaling factors for the first and the last page
    hscaleNew = (0.9*screenwidth)/(480*colsVec[currPage])
    vscaleNew = (0.8*screenheight-sliderBottomHeight-scaleHeight-
         ((as.integer(tclvalue(cbVar))+1)%%2)*rowsVec[currPage]*sliderheight-
         rowsVec[currPage]*3)/(480*rowsVec[currPage])
    hscaleNewLast = (0.9*screenwidth)/(480*colshelp)
    vscaleNewLast = (0.8*screenheight-sliderBottomHeight-scaleHeight-
         ((as.integer(tclvalue(cbVar))+1)%%2)*rowshelp*sliderheight-
         rowshelp*3)/(480*rowshelp)
# if aspect ratio has to be kept: adjust scaling factors
    if (as.integer(tclvalue(ksVar))==1) {  
      helpVec = c(-1,-1)    
      index = 1      
      for (i in 1:3){# determine vector of the non chosen viewing directions
        if (i!=view)  { 
          helpVec[index] = i
          index = index + 1  
        }
      }
      if (dt[helpVec[1]]>dt[helpVec[2]]){
        if (hscaleNew*(dt[helpVec[2]]/dt[helpVec[1]]) > vscaleNew){
          hscaleNew <- (dt[helpVec[1]]/dt[helpVec[2]])*vscaleNew
          hscaleNewLast <- (dt[helpVec[1]]/dt[helpVec[2]])*vscaleNewLast
        } else {
          vscaleNew = hscaleNew*(dt[helpVec[2]]/dt[helpVec[1]]) 
          vscaleNewLast = hscaleNewLast*(dt[helpVec[2]]/dt[helpVec[1]]) 
        }
      } else {
        if (hscaleNew*(dt[helpVec[1]]/dt[helpVec[2]]) > vscaleNew){
          hscaleNew <- (dt[helpVec[2]]/dt[helpVec[1]])*vscaleNew
          hscaleNewLast <- (dt[helpVec[2]]/dt[helpVec[1]])*vscaleNewLast
        } else {
          vscaleNew = hscaleNew*(dt[helpVec[1]]/dt[helpVec[2]]) 
          hscaleNewLast = vscaleNewLast*(dt[helpVec[1]]/dt[helpVec[2]])        
        }
      }  
    }
# minimal rescale of the slices, if the scaling factors are unchanged (compared to the last run) (otherwise: little error)  
    curr = 1
    while (curr!=(nrslices+1)){
      if ((curr <= (nrslicesVec[1]*(nrpages-1)) && 
          (hscaleNew == hscaleList[sliceslist[curr]]) && 
          (vscaleNew == vscaleList[sliceslist[curr]])) || 
          (curr > (nrslicesVec[1]*(nrpages-1)) && 
          (hscaleNewLast == hscaleList[sliceslist[curr]]) && 
          (vscaleNewLast == vscaleList[sliceslist[curr]]))) {
        hscaleNew     = 0.99*hscaleNew  
        vscaleNew     = 0.99*vscaleNew
        hscaleNewLast = 0.99*hscaleNewLast
        vscaleNewLast = 0.99*vscaleNewLast
      } else curr = curr + 1
    }
# deviance from quadrat not too big
    if (as.integer(tclvalue(ksVar))==0) {
      if (hscaleNew>=1.2*vscaleNew) hscaleNew <- 1.2*vscaleNew
      if (vscaleNew>=1.2*hscaleNew) vscaleNew <- 1.2*hscaleNew
      if (hscaleNewLast>=1.2*vscaleNewLast) hscaleNewLast <- 1.2*vscaleNewLast
      if (vscaleNewLast>=1.2*hscaleNewLast) vscaleNewLast <- 1.2*hscaleNewLast
    }        
    if (type=="spm") {# reset threshold value
      posThreshold <<- as.numeric(tclvalue(posvThreshold[[1]]))
      for (i in 1:nrpages) posThresHelp[i] <<- posThreshold
    }
    if (type=="data") {# reset time value
      posTime <<- as.numeric(tclvalue(posvTime[[1]]))
      for (i in 1:nrpages) posTimeHelp[i] <<- posTime
    }
# replot new slices with scaling factors; slider positons new set
    for (i in 1:(nrslicesVec[1]*(nrpages-1))){# pages till nrpages-1
      if ((nrslicesVec[1]*(nrpages-1)) != 0){
        if (view==1) {
          posNew[sliceslist[i]] <<- sliceslist[i]  
          tclvalue(posvNew[[sliceslist[i]]]) <<- sliceslist[i]  
        }  
        if (view==2) { 
          posNew[sliceslist[i]]<<- (sliceslist[i]-dt[1])
          tclvalue(posvNew[[sliceslist[i]]]) <<- (sliceslist[i]-dt[1])
        }
        if (view==3) {
          posNew[sliceslist[i]] <<- (sliceslist[i]-sum(dt[1:2]))
          tclvalue(posvNew[[sliceslist[i]]]) <<- (sliceslist[i]-sum(dt[1:2]))
        }
        tkrreplot(img[[sliceslist[i]]],hscale=hscaleNew,vscale=vscaleNew)
        hscaleList[sliceslist[i]] <<- hscaleNew  
        vscaleList[sliceslist[i]] <<- vscaleNew  
      }
    }
    for (i in (nrslicesVec[1]*(nrpages-1)+1):nrslices){# Page number nrpages
      if (view==1) {
        posNew[sliceslist[i]] <<- sliceslist[i]  
        tclvalue(posvNew[[sliceslist[i]]]) <<- sliceslist[i]  
      }  
      if (view==2) { 
        posNew[sliceslist[i]] <<- (sliceslist[i]-dt[1])
        tclvalue(posvNew[[sliceslist[i]]]) <<- (sliceslist[i]-dt[1])
      }
      if (view==3) {
        posNew[sliceslist[i]] <<- (sliceslist[i]-sum(dt[1:2]))
        tclvalue(posvNew[[sliceslist[i]]]) <<- (sliceslist[i]-sum(dt[1:2]))
      }
      tkrreplot(img[[sliceslist[i]]],hscale=hscaleNewLast,vscale=vscaleNewLast)
      hscaleList[sliceslist[i]] <<- hscaleNewLast  
      vscaleList[sliceslist[i]] <<- vscaleNewLast  
    }
# order all slices
    for (j in 1:rowsVec[1]){
      for (l in 1:colsVec[1]){
        if ((j-1)*colsVec[1]+l-1 < nrslicesVec[1]){
          tkgrid.configure(img[[sliceslist[(j-1)*colsVec[1]+l]]],column=l,row=2*j-1)
          tkgrid.configure(s[[sliceslist[(j-1)*colsVec[1]+l]]],column=l,row=2*j)
        }
      }
    }
# technical things
    for (i in 1:nrslices){
      if (view==1) tclvalue(posvNew[[sliceslist[i]]]) <<- sliceslist[i]    
      if (view==2) tclvalue(posvNew[[sliceslist[i]]]) <<- (sliceslist[i]-dt[1])
      if (view==3) tclvalue(posvNew[[sliceslist[i]]]) <<- (sliceslist[i]-sum(dt[1:2]))
      tkconfigure(label2[[sliceslist[i]]],
         text=as.integer(tclvalue(posvNew[[sliceslist[i]]])))# label infroamtion updated
      tkconfigure(label1[[sliceslist[i]]],text=label[view])        
    }  
    if (as.integer(tclvalue(cbVar))==1) hide()# hide sliders if wished
  }  
# is called by the right-arrow-button
# goes to next page and shows the next slices
  nextPage <- function(){
    oldPage <<- currPage
    currPage <<- currPage + 1
    if (type=="spm") { 
      tclvalue(posvThreshold[[1]]) <<- posThresHelp[currPage]
      tkconfigure(label2Th[[1]],text=posThresHelp[currPage])
    }
    if (type=="data") {
      tclvalue(posvTime[[1]]) <<- posTimeHelp[currPage]
      tkconfigure(label2Ti[[1]],text=posTimeHelp[currPage])
    }
# remove old slices and slider
    for (j in ((oldPage-1)*nrslicesVec[1]+1):
                  ((oldPage-1)*nrslicesVec[1]+nrslicesVec[oldPage])){
      tkgrid.forget(img[[sliceslist[j]]])
      tkgrid.forget(s[[sliceslist[j]]])    
    }
# remove old buttons and scales and readd buttons and scales (necessary, since otherwise: displacement) 
    if (nrpages>1){
      if (oldPage==1 || oldPage==nrpages || currPage==1 || currPage==nrpages){
        tkgrid.forget(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
                      nrSlicesppEntry,changeButton,viewAllButton,nextPageButton)
        if (currPage!=1){
          if (currPage!=nrpages){
            tkgrid(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
               nrSlicesppEntry,changeButton,viewAllButton,nextPageButton,padx=10,pady=5)    
          } else {
            tkgrid(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
               nrSlicesppEntry,changeButton,viewAllButton,padx=10,pady=5)
          }      
        } else {
          tkgrid(nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,nrSlicesppEntry,
                 changeButton,viewAllButton,nextPageButton,padx=10,pady=5)
        }
      }
    }
# add new slices and sliders
    for (j in 1:rowsVec[currPage]){
      for (l in 1:colsVec[currPage]){
        if ((j-1)*colsVec[currPage]+l-1 < nrslicesVec[currPage]){
          tkgrid.configure(img[[sliceslist[(currPage-1)*nrslicesVec[1]+
                           (j-1)*colsVec[currPage]+l]]],column=l,row=2*j-1) 
# note, that despite the last page every page has the same number of slices
          tkgrid.configure(s[[sliceslist[(currPage-1)*nrslicesVec[1]+
                           (j-1)*colsVec[currPage]+l]]],column=l,row=2*j)
        }
      }
    }
    hide()# sliders hidden, if necessary
  }
# is called by the left-arrow-button
# goes to previous page and shows the previous slices
  prevPage <- function(){
    oldPage <<- currPage
    currPage <<- currPage - 1
    if (type=="spm") { 
      tclvalue(posvThreshold[[1]]) <<- posThresHelp[currPage]
      tkconfigure(label2Th[[1]],text=posThresHelp[currPage])
    }
    if (type=="data") {
      tclvalue(posvTime[[1]]) <<- posTimeHelp[currPage]
      tkconfigure(label2Ti[[1]],text=posTimeHelp[currPage])
    }
# remove old slices and slider
    for (j in ((oldPage-1)*nrslicesVec[1]+1):
               ((oldPage-1)*nrslicesVec[1]+nrslicesVec[oldPage])){
      tkgrid.forget(img[[sliceslist[j]]])
      tkgrid.forget(s[[sliceslist[j]]])    
    }
# remove old buttons and scales and readd buttons and scales (necessary, since otherwise: displacement) 
    if (nrpages>1){
      if (oldPage==1 || oldPage==nrpages || currPage==1 || currPage==nrpages){
        tkgrid.forget(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
                      nrSlicesppEntry,changeButton,viewAllButton,nextPageButton)
        if (currPage!=1){
          if (currPage!=nrpages){
            tkgrid(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
               nrSlicesppEntry,changeButton,viewAllButton,nextPageButton,padx=10,pady=5)    
          } else {
            tkgrid(prevPageButton,nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,
                      nrSlicesppEntry,changeButton,viewAllButton,padx=10,pady=5)
          }      
        } else {
          tkgrid(nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,nrSlicesppEntry,
                 changeButton,viewAllButton,nextPageButton,padx=10,pady=5)
        }
      }
    }
# add new slices and sliders
    for (j in 1:rowsVec[currPage]){
      for (l in 1:colsVec[currPage]){
        if ((j-1)*colsVec[currPage]+l-1 < nrslicesVec[currPage]){
          tkgrid.configure(img[[sliceslist[(currPage-1)*nrslicesVec[1]+
                           (j-1)*colsVec[currPage]+l]]],column=l,row=2*j-1) 
#note, that despite the last page every page has the same number of slices
          tkgrid.configure(s[[sliceslist[(currPage-1)*nrslicesVec[1]+
                           (j-1)*colsVec[currPage]+l]]],column=l,row=2*j)
        }
      }
    }
    hide()# sliders hidden, if necessary
  }
  
# calls the 3-dimensional view
# is called by the button "View 3d"  
  call3d <- function(){
    plot.fmridata(localx,ttt,type="3d")
  }

  adjustContrast <- function(){
    okAdjCon <- function(){ 
      tkdestroy(ttAdjCon) 
      tkdestroy(tt)
      plot(localx,inputStuff[[1]],maxpvalue=inputStuff[[2]],cutOff=
               c(as.numeric(tclvalue(minVal)),as.numeric(tclvalue(maxVal))))      
    }
    quitAdjCon <- function(){ tkdestroy(ttAdjCon) }
    ttAdjCon <- tktoplevel(bg=wiasblue)
    tkwm.title(ttAdjCon, "Adjust contrast")
    adjConFrame1 <- tkframe(ttAdjCon,bg=wiasblue)
    adjConFrame2 <- tkframe(ttAdjCon,bg=wiasblue)
    adjConFrame3 <- tkframe(ttAdjCon,bg=wiasblue)
    adjConFrame4 <- tkframe(ttAdjCon,bg=wiasblue)
    minVal <- tclVar(0)
    maxVal <- tclVar(1)
    minEntry <- tkentry(adjConFrame2,textvariable=minVal,bg="#FFF",width=8)
    maxEntry <- tkentry(adjConFrame3,textvariable=maxVal,bg="#FFF",width=8)
    okAdjConButton <- tkbutton(adjConFrame4,text="Ok",command=okAdjCon,bg=wiaslightblue)
    quitAdjConButton <- tkbutton(adjConFrame4,text="Quit",
                                 command=quitAdjCon,bg=wiaslightblue)  
    tkgrid(tklabel(adjConFrame1,text=
       "Determine the value of the lower and the upper cutoff depending on the maximum 
       value: \n Re estimation can last a moment.",bg=wiasblue),padx=10,pady=10)
    tkgrid(tklabel(adjConFrame2,text="Lower Cutoff ",
                   bg=wiasblue),minEntry,padx=10,pady=10)
    tkgrid(tklabel(adjConFrame3,text="Upper Cutoff ",
                   bg=wiasblue),maxEntry,padx=10,pady=10)
    tkgrid(okAdjConButton,quitAdjConButton,padx=20,pady=10)
    tkgrid(adjConFrame1)
    tkgrid(adjConFrame2)
    tkgrid(adjConFrame3)
    tkgrid(adjConFrame4)
  }
  adjustMask <- function(){
    dataFile  <- ""
                nrrow <- numeric(1)
                nrcol <- numeric(1)           
    selectDataFile <- function(){
      tclvalue(dataFileTcl) <- tkgetOpenFile(filetypes =
      "{{ANALYZE} {.IMG .Img .img .HDR .Hdr .hdr}} {{AFNI} {.BRIK .Brik .brik .HEAD .Head .head}} {{NIFTI} {.NII .Nii .nii .HDR .Hdr .hdr}} {{All files} *}",
      title="Select data")  
    }
    loadDataHelp <- function(){
      if (as.character(tclvalue(dataFileTcl))=="") {
        quitttError <- function(){ tkdestroy(ttError) }
        ttError <- tktoplevel(bg=wiasblue)
        tkwm.title(ttError, "Error")
        tkgrid(tklabel(ttError,text="No file was selected.",bg=wiasblue),padx=10,pady=10)
        tkgrid(tkbutton(ttError,text="Ok",command=quitttError,
                        bg=wiaslightblue),padx=10,pady=10)
      } else loadData()
    }
    loadData <- function(){
      dataLoaded <<- TRUE
      dataFile <- tclvalue(dataFileTcl)
      help <- tolower(unlist(strsplit(dataFile,"")))
      help2 <- unlist(strsplit(dataFile,""))
      nrChars <- length(help)
      if (help[nrChars-2]=="i"&& help[nrChars-1]=="m" && help[nrChars]=="g")
        dataType = "ANALYZE"        
      else if (help[nrChars-2]=="n"&& help[nrChars-1]=="i" && help[nrChars]=="i")
        dataType = "NIFTI"
      else if ((help[nrChars-3]=="h"&& help[nrChars-2]=="e" && help[nrChars-1]=="a" &&
                help[nrChars]=="d")||(help[nrChars-3]=="b"&& help[nrChars-2]=="r" &&
                help[nrChars-1]=="i" && help[nrChars]=="k"))
        dataType = "AFNI"
      else if (help[nrChars-2]=="h"&& help[nrChars-1]=="d" && help[nrChars]=="r"){
        if (file.info(dataFile)$size == 348) dataType = "ANALYZE"
        else dataType = "NIFTI"
      } else dataType = "unknown"  
      if (dataType=="AFNI") data <<- read.AFNI(dataFile)
      if (dataType=="ANALYZE") data <<- read.ANALYZE(dataFile)
      if (dataType=="NIFTI") data <<- read.NIFTI(dataFile)
      if (dataType=="unknown"){
        print("The data type is unknown !!")
        print("Please check your path or press 'help'.")
      }
      tclvalue(thresValTcl) <- round(quantile(extract.data(data),0.75),2)
    }
    viewMaskHelp <- function(){
      if (!dataLoaded) {
        quitttError <- function(){ tkdestroy(ttError) }
        ttError <- tktoplevel(bg=wiasblue)
        tkwm.title(ttError, "Error")
        tkgrid(tklabel(ttError,text="No file was loaded.",bg=wiasblue),padx=10,pady=10)
        tkgrid(tkbutton(ttError,text="Ok",command=quitttError,
                        bg=wiaslightblue),padx=10,pady=10)
      } else viewMask()
    }
    viewMask <- function(){
      ttt <- extract.data(data)
      ddim <- dim(ttt)    
      dev.new(width=12,height=7)
      if (round(sqrt(ddim[3]))==sqrt(ddim[3])){
        nrrow <<- sqrt(ddim[3])
        nrcol <<- sqrt(ddim[3])  
      } else {
        if ((ceiling(sqrt(ddim[3]))-1)*ceiling(sqrt(ddim[3])) >= ddim[3]) {
          nrrow <<- ceiling(sqrt(ddim[3]))-1
          nrcol <<- ceiling(sqrt(ddim[3]))
        } else {
          nrrow <<- ceiling(sqrt(ddim[3]))
          nrcol <<- ceiling(sqrt(ddim[3])) 
        }
      }
      mat = matrix(0,nrrow,nrcol+1)
      for (i in 1:nrrow) for (j in 1:(nrcol+1)) if ((i-1)*(nrcol+1)+j-(i-1) <= ddim[3])
          mat[i,j]=(i-1)*(nrcol+1)+j-(i-1)   
      for (i in 1:nrrow){ mat[i,nrcol+1] = ddim[3]+1 }
      widthsvec = c(1:nrcol+1)
      for (i in 1:nrcol){ widthsvec[i]=0.5/nrcol }
      widthsvec[nrcol+1] = 0.5
      layout(mat,widthsvec)
      par(mar=c(0.5,0.5,0.5,0.5))
      for (i in 1:ddim[3])
         image(ttt[,,i,1]>as.numeric(tclvalue(thresValTcl)),yaxt="n",xaxt="n")
      par(mar=c(5,5,3,1))
      bwV = diff(range(ttt))/(length(ttt[,,,1])/1200)
      d0 <- density(ttt[,,,1],bw=bwV)  
      d1 <- density(ttt[round((1/8)*ddim[1]):round((7/8)*ddim[1]),round((1/8)*ddim[2]):
               round((7/8)*ddim[2]),round((1/8)*ddim[3]):round((7/8)*ddim[3]),1],bw=bwV)
      d2 <- density(ttt[round((2/8)*ddim[1]):round((6/8)*ddim[1]),round((2/8)*ddim[2]):
               round((6/8)*ddim[2]),round((2/8)*ddim[3]):round((6/8)*ddim[3]),1],bw=bwV)  
      d3 <- density(ttt[round((3/8)*ddim[1]):round((5/8)*ddim[1]),round((3/8)*ddim[2]):
               round((5/8)*ddim[2]),round((3/8)*ddim[3]):round((5/8)*ddim[3]),1],bw=bwV)  
      plot(d0,main="")
      title(main="Density plots",cex.main=1.5)
      lines(d1,col=2)
      lines(d2,col=3)
      lines(d3,col=4)    
      lines(c(as.numeric(tclvalue(thresValTcl)),as.numeric(tclvalue(thresValTcl))),
            range(d0$y),col=6)
      legend(0.55*max(d0$x),0.98*max(d0$y),
         c("Data","Centered 75% of data","Centered 50% of data","Centered 25% of data",
         "Threshold line"),text.col=c(1,2,3,4,6),pch=c(1,1,1,1,1),col=c(1,2,3,4,6),
         title="Density of",cex=1.5)
    }
    okAdjMask <- function(){
      if (!dataLoaded){
        quitttError <- function(){ tkdestroy(ttError) }
        ttError <- tktoplevel(bg=wiasblue)
        tkwm.title(ttError, "Error")
        tkgrid(tklabel(ttError,text="No file was loaded.",bg=wiasblue),padx=10,pady=10)
        tkgrid(tkbutton(ttError,text="Ok",command=quitttError,
                        bg=wiaslightblue),padx=10,pady=10)
      }  else {
        tkdestroy(ttAdjMask) 
        quantile <- as.numeric(tclvalue(thresValTcl))
        anatomic=extract.data(data)[,,,1]
        anatomic[anatomic<quantile] <- 0   
        plot(localx,anatomic,maxpvalue=inputStuff[[2]],cutOff=inputStuff[[3]])
        tkdestroy(tt)
      }  
    }
    quitAdjMask <- function(){ tkdestroy(ttAdjMask) }
    ttAdjMask <- tktoplevel(bg=wiasblue)
    tkwm.title(ttAdjMask, "Adjust mask")
    dataFileTcl <- tclVar("")
    data <- list()
    dataLoaded <- FALSE
    adjMaskFrame1 <- tkframe(ttAdjMask,bg=wiasblue)
    adjMaskFrame2 <- tkframe(ttAdjMask,bg=wiasblue)
    adjMaskFrame3 <- tkframe(ttAdjMask,bg=wiasblue)
    adjMaskFrame4 <- tkframe(ttAdjMask,bg=wiasblue)
    adjMaskFrame5 <- tkframe(ttAdjMask,bg=wiasblue)
    adjMaskFrame6 <- tkframe(ttAdjMask,bg=wiasblue)  
    adjMaskFrame7 <- tkframe(ttAdjMask,bg=wiasblue)  
    objFileL    <- tklabel(adjMaskFrame1,text="Load data",
                           bg=wiasblue,font="Arial 13 bold")  
    objFileE1   <- tkentry(adjMaskFrame2, textvariable = dataFileTcl, 
                           width = 40, bg = "#ffffff")
    objFileB1   <- tkbutton(adjMaskFrame2, text = "Select file", width = 15, 
                     command = selectDataFile, bg = wiaslightblue, anchor = "c")  
    objFileLoad <- tkbutton(adjMaskFrame4, text = "Load", width = 15, 
                     command = loadDataHelp, bg = wiaslightblue)
    helpLabel1 <- tklabel(adjMaskFrame3,text="",
                     bg=wiasblue,width=0,font="Arial 1")
    thresLabel0 <- tklabel(adjMaskFrame5,text="Determine threshold",
                     bg=wiasblue,font="Arial 13 bold")  
    thresValTcl <- tclVar()
    thresEntry <- tkentry(adjMaskFrame6,textvariable=thresValTcl,width=8)
    thresLabel  <- tklabel(adjMaskFrame6,text="Threshold", 
                     bg=wiasblue,font="Arial 12 bold")
    viewMaskButton <- tkbutton(adjMaskFrame6,text="View mask",
                     command=viewMaskHelp, bg=wiaslightblue)
    okAdjMaskButton <- tkbutton(adjMaskFrame7,text="Ok",command=okAdjMask,
                     bg=wiaslightblue)
    quitAdjMaskButton <- tkbutton(adjMaskFrame7,text="Quit",
                     command=quitAdjMask,bg=wiaslightblue)  
    tkgrid(objFileL,pady = 10, padx = 10,sticky="ew")
    tkgrid(objFileB1,objFileE1,pady = 10, padx = 10)
    tkgrid(objFileLoad,pady=10,padx=15)
    tkgrid(thresLabel0,pady = 10, padx = 10,sticky="ew")
    tkgrid(thresLabel,thresEntry,viewMaskButton,padx=10,pady=10)
    tkgrid(okAdjMaskButton,quitAdjMaskButton,padx=15,pady=10)
    tkgrid(adjMaskFrame1)
    tkgrid(adjMaskFrame2)
    tkgrid(adjMaskFrame4)
    tkgrid(adjMaskFrame5)
    tkgrid(adjMaskFrame6)
    tkgrid(adjMaskFrame7)
  }
# save images with the help of adimpro
# is called via the button "Extract Images"
# first step: choose slices (from current viewing direction)
# second step: choose file format (only if ImageMagick is installed)
# third step; choose filename and directory
# if the filename is exists, antoher one can be chosen
  extractImages <- function(){
    if (!require(adimpro)) 
       stop("required package adimpro not found. Please install from cran.r-project.org")  
    quitttImgSelect <- function(){
      tkdestroy(ttImgSelect)
    }
    currentFileName = "unspecified"
    filePathFinal = "unspecified"
    chosenSlices <- c()
    slicesCounter <- 0    
# filename chosen by user
# only once called
    chooseFileNameGeneral <- function(){
      imageFile <- tclVar()  
      if (as.character(tclvalue(rbFF)) == "ppm")  tclvalue(imageFile) <- 
            tkgetSaveFile(filetypes ="{{Image file} {.ppm}}",title="Save Images")  
      if (as.character(tclvalue(rbFF)) == "jpeg") tclvalue(imageFile) <- 
            tkgetSaveFile(filetypes ="{{Image file} {.jpeg}}",title="Save Images")  
      if (as.character(tclvalue(rbFF)) == "png")  tclvalue(imageFile) <- 
            tkgetSaveFile(filetypes ="{{Image file} {.png}}",title="Save Images")  
      if (as.character(tclvalue(rbFF)) == "tif")  tclvalue(imageFile) <- 
            tkgetSaveFile(filetypes ="{{Image file} {.tif}}",title="Save Images")
# check file extension        
      help <- tolower(unlist(strsplit(as.character(tclvalue(imageFile)),"")))
      nrChars <- length(help)    
      if (nrChars <=4 || !((help[nrChars-3]=="." && help[nrChars-2]=="p" &&
      help[nrChars-1]=="p" && help[nrChars]=="m") || (help[nrChars-3]=="." &&
      help[nrChars-2]=="p" && help[nrChars-1]=="n" && help[nrChars]=="g") || 
      (help[nrChars-3]=="." && help[nrChars-2]=="t" && help[nrChars-1]=="i" &&
      help[nrChars]=="f") || (help[nrChars-4]=="." && help[nrChars-3]=="j" && 
      help[nrChars-2]=="p" && help[nrChars-1]=="e" && help[nrChars]=="g"))){
        filePath <- paste(as.character(tclvalue(imageFile)),
                          paste(".",as.character(tclvalue(rbFF)),sep=""),sep="")
      } else {
        filePath <- as.character(tclvalue(imageFile))
      }  
      help <- unlist(strsplit(filePath,""))
      nrChars <- length(help)
      index = -1
      if (help[nrChars-3]==".") index = (nrChars-3) else index = (nrChars-4)
      filePathFinal = ""
      for (j in 1:index-1) filePathFinal <- paste(filePathFinal,help[j],sep="")
      if (viewAxis==3) filePathFinal <- paste(filePathFinal,"Axial",sep="")
      if (viewAxis==2) filePathFinal <- paste(filePathFinal,"Sagittal",sep="")
      if (viewAxis==1) filePathFinal <- paste(filePathFinal,"Coronal",sep="")  
      filePathFinal <- paste(filePathFinal,"Slice",sep="")
      if (type=="data"){
        tmp <- paste("Time",tclvalue(posvTime[[1]]),sep="")
        filePathFinal <<- paste(filePathFinal,tmp,sep="")
      }
      if (type=="spm"){
        tmp <- paste("Threshold",tclvalue(posvThreshold[[1]]),sep="")
        filePathFinal <<- paste(filePathFinal,tmp,sep="")
      }
      if (type=="pvalue") filePathFinal <<- filePathFinal
      ending = ""
      for (j in (index+1):nrChars) ending <- paste(ending,help[j],sep="")
      tclvalue(rbFF) = ending
    }
# check if the filename is exists, if it is offer rechoosing of filename
# calls writeImage
# is called for every slices
    chooseFileNameExactly <- function(overwrite=FALSE,nrofChosen=1){
#nrofchosen is the position of the current processed slice in chosenSlices
      if (viewAxis==3) slicenr <- chosenSlices[nrofChosen]-dt[1]-dt[2]  
      if (viewAxis==2) slicenr <- chosenSlices[nrofChosen]-dt[1]
      if (viewAxis==1) slicenr <- chosenSlices[nrofChosen]      
      filePathFinal <- paste(filePathFinal,slicenr,sep="")
      filePathFinal <- paste(filePathFinal,
          paste(".",as.character(tclvalue(rbFF)),sep=""),sep="")#.fileformat appended
      if (overwrite == FALSE){# check if name exists
        pathSplitted <- unlist(strsplit(filePathFinal,""))
        noSeparator = TRUE
        indexSep = length(pathSplitted)              
        while (noSeparator){
          if (pathSplitted[indexSep]==.Platform$file.sep) noSeparator = FALSE
          indexSep = indexSep - 1
        }
        pathCutted = ""
        cuttedPart = ""
        for (k in 1:indexSep) pathCutted <- paste(pathCutted,pathSplitted[k],sep="")
        for (k in (indexSep+2):length(pathSplitted)) 
          cuttedPart <- paste(cuttedPart,pathSplitted[k],sep="")
        listFiles <- list.files(pathCutted)
        conflict = FALSE
        for (k in 1:length(listFiles)){
          if (listFiles[k]==cuttedPart) {
            conflict = TRUE
            quitttWarn <- function(){
              tkdestroy(ttWarn)
              currentFileName <<- filePathFinal
              writeImage(nrofChosen)
              overwrite = TRUE                    
            }
            otherName <- function(){
              tkdestroy(ttWarn)
              chooseFileNameGeneral()
              chooseFileNameExactly(FALSE,nrofChosen)
            }
            print("Name exists !!")
            ttWarn <- tktoplevel(bg=wiasblue)
            tkwm.title(ttWarn, "Warning")
            framettWarn1 <- tkframe(ttWarn,bg=wiasblue)
            framettWarn2 <- tkframe(ttWarn,bg=wiasblue)
            contButton <- tkbutton(framettWarn2,text="Continue",command=quitttWarn,
                                   bg=wiaslightblue)
            alternButton <- tkbutton(framettWarn2,text="Other filename",
                                   command=otherName,bg=wiaslightblue)
            tkgrid(tklabel(framettWarn1,bg=wiasblue,text= 
                       paste("Do you really want to overwrite the file?",cuttedPart))) 
            tkgrid(framettWarn1,pady=5)
            tkgrid(contButton,alternButton,padx=10)
            tkgrid(framettWarn2,padx=10,pady=5)
          }                
        }
        if (conflict == FALSE) {
          currentFileName <<- filePathFinal  
          writeImage(nrofChosen)  
        }
      } else {
        currentFileName <<- filePathFinal    
        writeImage(nrofChosen)        
      }              
    }
# write the current slice into a previously determined file
# is called from chooseFileNameExactly    
    writeImage <- function(nrofChosen){
      if (length(dt)==4){
        if (type=="data"){
          timeVal <- as.numeric(tclvalue(posvTime[[1]]))
          if (viewAxis==3) tmp <- make.image(ttt[,dt[2]:1,
                 chosenSlices[nrofChosen]-dt[1]-dt[2],timeVal],gammatype="ITU")
          if (viewAxis==2) tmp <- make.image(ttt[,
                 chosenSlices[nrofChosen]-dt[1],,timeVal],gammatype="ITU")
          if (viewAxis==1) tmp <- make.image(ttt[chosenSlices[nrofChosen],,,timeVal],
                                             gammatype="ITU")
        } else {# threshold korrektur koennte man noch einbauen
          if (viewAxis==3) tmp <- make.image(ttt[,dt[2]:1,
                      chosenSlices[nrofChosen]-dt[1]-dt[2],1],gammatype="ITU")
          if (viewAxis==2) tmp <- make.image(ttt[,
                           chosenSlices[nrofChosen]-dt[1],,1],gammatype="ITU")
          if (viewAxis==1) tmp <- make.image(ttt[chosenSlices[nrofChosen],,,1],
                                              gammatype="ITU")
        }
      } else {  
        if (type=="spm"){
          thresVal <- (as.numeric(tclvalue(posvThreshold[[1]])) - scale[1])/diff(scale)
          tttHelp <- ttt
          tttHelp[tttHelp < thresVal] <- 0
          if (viewAxis==3) tmp <- make.image(tttHelp[,dt[2]:1,
                           chosenSlices[nrofChosen]-dt[1]-dt[2]],gammatype="ITU")
          if (viewAxis==2) tmp <- make.image(tttHelp[,
                               chosenSlices[nrofChosen]-dt[1],],gammatype="ITU")
          if (viewAxis==1) tmp <- make.image(tttHelp[chosenSlices[nrofChosen],,],
                                             gammatype="ITU")
        } else{
          if (viewAxis==3) tmp <- ttt[,dt[2]:1,chosenSlices[nrofChosen]-dt[1]-dt[2]]
          if (viewAxis==2) tmp <- ttt[,chosenSlices[nrofChosen]-dt[1],]
          if (viewAxis==1) tmp <- ttt[chosenSlices[nrofChosen],,]
          rgbcolors <-col2rgb(col)/255
          ncolors <- length(col)
          ctmp <- array(0,c(dim(tmp),3))
          ctmp[,,1] <- rgbcolors[1,trunc(tmp*(ncolors-1)+1)]
          ctmp[,,2] <- rgbcolors[2,trunc(tmp*(ncolors-1)+1)]
          ctmp[,,3] <- rgbcolors[3,trunc(tmp*(ncolors-1)+1)]
          tmp <- make.image(ctmp, gammatype="ITU")
        }
      }
      write.image(tmp,currentFileName)
      if (nrofChosen < slicesCounter) chooseFileNameExactly(FALSE,nrofChosen+1)
    }
# is called via the next button after choosing the file format
# determines the filename with the help of chooseFileNameGeneral and chooseFileNameExactly
    onNext <- function(){              
      for (i in 1:dt[viewAxis]){
        if (as.integer(tclvalue(cbVars[[i]]))==1) {
          slicesCounter <<- slicesCounter +1
          if (viewAxis==1) chosenSlices[slicesCounter] <<- i
          if (viewAxis==2) chosenSlices[slicesCounter] <<- dt[1] + i
          if (viewAxis==3) chosenSlices[slicesCounter] <<- dt[1] + dt[2] + i  
        }
      }
      if (slicesCounter>0){
        chooseFileNameGeneral()
        chooseFileNameExactly(FALSE,1)
      } else {# error
        quitttError <- function() tkdestroy(ttError) 
        ttError <- tktoplevel(bg=wiasblue)
        tkwm.title(ttError, "Error")
        tkgrid(tklabel(ttError,text="No file was selected.",bg=wiasblue),padx=10,pady=10)
        tkgrid(tkbutton(ttError,text="Ok",command=quitttError,
                          bg=wiaslightblue),padx=10,pady=10)
      }
#writeImage called from chooseFileName
    }
# is called after choosing the slices
# gives the user the possibility to choose the file format
# needs imgmagick, otherwise: file format == ppm (automatically)
    toFileFormat <- function(){
      tkdestroy(ttImgSelect)# Menu1 closed 
      onNextHelp <- function(){
        tkdestroy(ttFF)# Menu2 closed
        onNext()      
      }
      imgmagickInstalled = FALSE
      if (.Platform$OS.type == "windows") {
        a <- system(paste(convert.path,"-version"),FALSE)
        if (a >= 0) {
          imgmagickInstalled = TRUE  
        }    
      } else {# unix
        if (.Platform$OS.type != "unix") 
           warning("never tested this OS. maybe we cannot proceed here.\n")
        a <- system("convert -version",TRUE,TRUE)
        if (length(grep("ImageMagick",a,ignore.case=TRUE)) > 0) imgmagickInstalled = TRUE          
      }    
      if (imgmagickInstalled){# choose file format
        ttFF = tktoplevel(bg=wiasblue)
        tkwm.title(ttFF, "Extract images - file format")  
        frameFF1 <- tkframe(ttFF,relief="groove",borderwidth=0,bg=wiasblue)
        frameFF2 <- tkframe(ttFF,relief="groove",borderwidth=0,bg=wiasblue)
        frameFF3 <- tkframe(ttFF,relief="groove",borderwidth=0,bg=wiasblue)
        rbFF1 <- tkradiobutton(frameFF2,bg=wiasblue)
        rbFF2 <- tkradiobutton(frameFF2,bg=wiasblue)
        rbFF3 <- tkradiobutton(frameFF2,bg=wiasblue)
        rbFF4 <- tkradiobutton(frameFF2,bg=wiasblue)
        nextButton <- tkbutton(frameFF3,text="next ",command=onNextHelp,
                               bg=wiaslightblue)
        tkconfigure(rbFF1,variable=rbFF,value="ppm",bg=wiasblue)
        tkconfigure(rbFF2,variable=rbFF,value="jpeg",bg=wiasblue)
        tkconfigure(rbFF3,variable=rbFF,value="png",bg=wiasblue)
        tkconfigure(rbFF4,variable=rbFF,value="tif",bg=wiasblue)  
        tkgrid(tklabel(frameFF1,text="Choose your preferred file format",
                       bg = wiasblue,font="Arial 12 bold"))
        tkgrid(frameFF1,sticky="ew")
        tkgrid(tklabel(frameFF2,text="ppm ",bg = wiaslightblue,font="Arial 12"),
               rbFF1,tklabel(frameFF2,text="jpeg ",bg = wiaslightblue,font="Arial 12"),
               rbFF2, padx=30,pady=10)
        tkgrid(tklabel(frameFF2,text="png ",bg = wiaslightblue,font="Arial 12"),
               rbFF3,tklabel(frameFF2,text="tif ",bg = wiaslightblue,font="Arial 12"),
               rbFF4, padx=30,pady=10)
        tkgrid(frameFF2,sticky="ew")
        tkgrid(nextButton,padx=30,pady=10)
        tkgrid(frameFF3)  
      } else {
        tclvalue(rbFF) <- "ppm"
        onNext()
      }
    }
# can be called in step1 by the button "Select All"
    selectAll <- function() for (i in 1:dt[viewAxis]) tclvalue(cbVars[[i]]) <- 1
# can be called in step1 by the button "Deselect All"
    deselectAll <- function() for (i in 1:dt[viewAxis]) tclvalue(cbVars[[i]]) <- 0
    ttImgSelect <- tktoplevel(bg=wiasblue)# toplevel window for the slices choice
    tkwm.title(ttImgSelect, "Extract images - slice selection")       
    imgSelectFrom <- list()
    viewAxis = as.integer(tclvalue(rbValue))
# determine number of rows and cols in layout
    cols <- rows <- ceiling(sqrt(dt[viewAxis]+1))
    while ((rows-1)*cols >= (dt[viewAxis]+1)) rows <- rows-1
    if (((cols-1)*rows==dt[viewAxis])) cols <- cols-1
    help <- 0    
    if (cols<rows){
      help = cols
      cols <- rows
      rows <- help
    }         
    hsc = (0.9*as.numeric(tkwinfo("screenwidth", ttImgSelect)))/(480*cols)
    vsc = (0.8*as.numeric(tkwinfo("screenheight", ttImgSelect))-rows*20)/(480*rows)
    imgFrame <- tkframe(ttImgSelect,bg=wiasblue)
# create necessary slices
    if (viewAxis == 1) for (i in 1:dt[1]) 
      imgSelectFrom[[i]] = fmri.image("x",hsc,vsc,i,toplevel=imgFrame,overview=TRUE)
    if (viewAxis == 2) for (i in 1:dt[2]) 
      imgSelectFrom[[i]] = fmri.image("y",hsc,vsc,i,toplevel=imgFrame,overview=TRUE)
    if (viewAxis == 3) for (i in 1:dt[3]) 
      imgSelectFrom[[i]] = fmri.image("z",hsc,vsc,i,toplevel=imgFrame,overview=TRUE)
    helpFunc1 <- function(a) a <- tclVar("0")  
    helpFunc2 <- function(a) a <- tkcheckbutton(imgFrame,
                     text="",variable=cbVars[[as.integer(a)]],bg="#BBDDEC")
    cbVars <- lapply(1:dt[viewAxis], helpFunc1)
    cboxSelects <- lapply(1:dt[viewAxis], helpFunc2)    
# currently viewed slices preselected
    nrslicesCurr <- as.integer(tclvalue(nrSlices))  
    for (i in 1:dt[viewAxis])
      for (j in 1:nrslicesCurr) 
        if (tclvalue(posvNew[[sliceslist[j]]])==i) tclvalue(cbVars[[i]]) = 1              
# add all necessary slices
    for (j in 1:rows){
      for (l in 1:cols){
        if ((j-1)*cols+l-1 < dt[viewAxis]){
          tkgrid.configure(imgSelectFrom[[(j-1)*cols+l]],column=l,row=2*j-1) 
          tkgrid.configure(cboxSelects[[(j-1)*cols+l]],column=l,row=2*j)
        }
      }
    }
    tkgrid.configure(tklabel(imgFrame,text="",bg=wiasblue),column=1,row=2*rows+1)
    tkgrid(imgFrame)
# buttons
    bottomFrame <- tkframe(ttImgSelect,bg=wiasblue)
    next2Button <- tkbutton(bottomFrame,text="Next",command = toFileFormat,
                            bg=wiaslightblue)
    selAllButton <- tkbutton(bottomFrame,text="Select all",command = selectAll,
                            bg=wiaslightblue)
    deselAllButton <- tkbutton(bottomFrame,text="Deselect all",command = deselectAll, 
                            bg=wiaslightblue)
    quitttImgButton <- tkbutton(bottomFrame,text="Quit",command = quitttImgSelect,
                            bg=wiaslightblue)
    tkgrid(next2Button,selAllButton,deselAllButton,quitttImgButton,padx=17,pady=5)
    tkgrid(bottomFrame)
    tkgrid.configure(tklabel(ttImgSelect,text="",bg=wiasblue),column=1,row=2*rows+3)
    rbFF <- tclVar("ppm")  
  }
# is called by the button "View Help" 
  helpFunction <- function(){  
    onQuit <- function() tkdestroy(ttHelp)
    ttHelp = tktoplevel(bg=wiasblue)
    tkwm.title(ttHelp, "Help")
    helptextIntr = "Generic function plot for objects of classes ''fmridata'' (fmri data objects), ''spm'' (statistical parametric maps) and ''fmripvalue'' \n \n \n"  
    helptextViewOpt="Options \n 
      \n You can choose between coronal, sagittal and axial slices. Furthermore you may 
      \n decide how many slices shall be selected and how many of these are shown at once. If you 
      \n choose the number of slices shown at once smaller than the number of selected slices, there 
      \n will appear an arrow on the right side. With the help of this arrow you can go to next 
      \n page, and the next slices will be presented. To confirm your choice on the number of slices 
      \n and the viewing direction press the button 'Change View/Slices'. If you are not on the first
      \n page, you can go back to this by using the arrow on the left side. Below each
      \n slice there is a slider, which you can use to slide between the slices of the current
      \n viewing direction. To remove the sliders press 'Hide Sliders'. By selecting 'Keep aspect 
      \n ratio' the original heigth-to-width ratio is rebuilt. If this is unselected, the slice will
      \n be plotted nearly quadratic. To navigate through coronal, sagittal and axial slices at the same time press
      \n 'View 3d'. In the first row below the slices there is printed a scale, which interpolates
      \n between the colours red(lowest value) and white(highest value) or grey(lowest value) and
      \n white (highest value). If you have plotted the statistical parametric map, there will be a slider in the 
      \n last row which can be used to determine the threshold. For fmri data objects  the time can be 
      \n chosen by a slider.
      \n Using the button 'Adjust Mask' you can determine a threshold and create an underlying mask. 
      \n If you have problems with the brightness you can use 'Adjust contrast'. This may lead too a 
      \n more balanced colour distribution. \n 
      \n"
    helptextSave="Save results \n 
      \n If you want to save your results the button 'Extract Images' gives you a comfortable 
      \n possibility to do so. First you have to select all desired slices, To select/ deselect all
      \n slices use the corresponding button. To continue use 'Next'. If Image Magick is installed on
      \n your system you can choose between the image types ppm, jpeg, png and tif. If not ppm will
      \n automatically be taken. In the last step you have to choose a filename and filedirectory.
      \n To ensure not to overwrite existing data, it will be checked that this filedirectory/
      \n filename does not exist. If it exists a warning is issued with an option to change the filename.
      \n The images will then be created in the chosen directory as [chosen 
      \n filename][viewing direction]Slice[slicenumber].[filetype]."    
    helptext = paste(paste(helptextIntr,helptextViewOpt,sep=""),helptextSave,sep="")
    helpFrame1 = tkframe(ttHelp,bg=wiasblue)
    helpLabel = tklabel(helpFrame1,text=helptext,font="Arial 13",bg=wiasblue)
    helpFrame2 = tkframe(ttHelp,bg=wiasblue)  
    helpB1 = tkbutton(helpFrame2,text="Quit",command=onQuit,bg=wiaslightblue)
    tkgrid(helpLabel)
    tkgrid(helpB1,padx=10,pady=10)
    tkgrid(helpFrame1)  
    tkgrid(helpFrame2)  
  }

# is called by the checkbutton "Hide Sliders"
# pressing this checkbutton changes the value of cbVar
# changeFunction changes necessary things with higher local value of cbVar  
  hidesliders <- function(){
    helpPage <- currPage
    changeFunction()# change all necessary things
    if (helpPage>1){# go to current page
      for (i in 1:(helpPage-1)) nextPage()
    }
#hide() wird in changeFunction aufgreufen
  }

# hide sliders if so selected
# can be called from changeFunction, nextPage and currentPage  
  hide <- function(){
    if (as.integer(tclvalue(cbVar))==1){# remove sliders
      for (i in ((currPage-1)*nrslicesVec[1]+1):
      ((currPage-1)*nrslicesVec[1]+nrslicesVec[currPage])){
        tkgrid.forget(s[[sliceslist[i]]])  
      }
    } else {# add sliders
      for (j in 1:rowsVec[currPage]){
        for (l in 1:colsVec[currPage]){
          if ((j-1)*colsVec[currPage]+l-1 < nrslicesVec[currPage]){
            tkgrid.configure(img[[sliceslist[(currPage-1)*nrslicesVec[1]+
                             (j-1)*colsVec[currPage]+l]]],column=l,row=2*j-1) 
#note, that despite the last page every page has the same number of slices
            tkgrid.configure(s[[sliceslist[(currPage-1)*nrslicesVec[1]+
                             (j-1)*colsVec[currPage]+l]]],column=l,row=2*j)
          }
        }
      }  
    }  
  }
# is called from checkbutton "Keep Aspect Ratio"
# pressing this checkbutton changes the value of ksVar
# changeFunction changes necessary things with higher local value of ksVar  
  keepsides <- function() changeFunction()
  keep <- function(){
}
  convert.path <- paste(Sys.getenv("ImageMagick"),"convert",sep="")
  posTimeHelp <- c(1)
  posThresHelp <- c(0)
  wiasblue <- "#AACCDB"
  wiaslightblue <- "#BBDDEC"
  rowsVec <- c(2)
  colsVec <- c(2)
  nrpages <- 0
  oldPage <- 1  
  currPage <- 1  
  nrslicesVec <- c(4)
  screenwidth <- 0
  screenheight <- 0
  sliderheight <- 0
  sliderwidth <- 0
  sliderBottomHeight <- 0
  scaleHeight <- 0
  prevPageButton <- tkbutton(frame1,text="<-",command=prevPage,bg="#BBDDEC")
  nextPageButton <- tkbutton(frame1,text="->",command=nextPage,bg="#BBDDEC")
  changeButton <- tkbutton(frame1,text=
                   "Change View/Slices",command=changeFunction,bg="#BBDDEC")
  view = ""
  var = as.numeric(tclvalue(rbValue))
  textViewAll = paste("View all ","",view,sep="")
  textViewAll = paste(textViewAll,""," slices",sep="")
  viewAllButton <- tkbutton(frame1,text=textViewAll,command=viewAll,bg="#BBDDEC")
  nrSlicesLabel <- tklabel(frame1,text="Number of slices",bg="#BBDDEC")
  nrSlices <- tclVar("4")
  nrSlicesEntry <- tkentry(frame1,textvariable=nrSlices,bg="#FFF",width=8)
  nrSlicesppLabel <- tklabel(frame1,text="Number of slices per page",bg="#BBDDEC")
  nrSlicespp <- tclVar("4")
  nrSlicesppEntry <- tkentry(frame1,textvariable=nrSlicespp,bg="#FFF",width=8)  
  rb1 <- tkradiobutton(frame0)
  rb2 <- tkradiobutton(frame0)
  rb3 <- tkradiobutton(frame0)
  tkconfigure(rb1,variable=rbValue,value=1,bg="#BBDDEC")
  tkconfigure(rb2,variable=rbValue,value=2,bg="#BBDDEC")
  tkconfigure(rb3,variable=rbValue,value=3,bg="#BBDDEC")
  cbVar <- tclVar("0")  
  ksVar <- tclVar("0")
  cboxSlider <- tkcheckbutton(frame2,text="Hide Sliders",variable=cbVar,
                              bg="#BBDDEC",command=hidesliders)  
  dButton <- tkbutton(frame2,text="View 3d",command = call3d,bg="#BBDDEC")
  extractButton <- tkbutton(frame2,text="Extract Images", 
                              command = extractImages,bg="#BBDDEC")
  adjustCButton <- tkbutton(frame2,text="Adjust Contrast", 
                              command = adjustContrast,bg="#BBDDEC")
  adjustMButton <- tkbutton(frame2,text="Adjust Mask", command = adjustMask,bg="#BBDDEC")
  helpButton <- tkbutton(frame2,text="View Help", command = helpFunction,bg="#BBDDEC")
  keepSides <- tkcheckbutton(frame2,text="Keep Aspect Ratio",variable=ksVar,
                             bg="#BBDDEC",command=keepsides)  
  tkgrid(tklabel(frame0,text="Coronal",bg="#BBDDEC"),rb1,
         tklabel(frame0,text="Sagittal",bg="#BBDDEC"),rb2,
         tklabel(frame0,text="Axial",bg="#BBDDEC"),rb3,padx=10,pady=5)
  tkgrid(frame0)
  tkgrid(nrSlicesLabel,nrSlicesEntry,nrSlicesppLabel,nrSlicesppEntry,
         changeButton,viewAllButton,padx=10,pady=5)
  tkgrid(frame1)  
  tkgrid(cboxSlider,keepSides,dButton,extractButton,adjustMButton,
         adjustCButton,helpButton,padx=10,pady=5) 
  tkgrid(frame2)  
  sliderBottomHeight <- 0
  if (type=="spm") {
    thres <- c(fmri.threshold())
    tkgrid(thres[[1]])
    sliderBottomHeight <- (as.numeric(tkwinfo("height",thres[[1]])) + 5)
  }
  if (type=="data") {
    time <- c(fmri.slider(-1.))
    tkgrid(time[[1]])
    sliderBottomHeight <- (as.numeric(tkwinfo("height",time[[1]])) + 5)
  }    
  for (i in 1:nrslices){ 
    posNew[sliceslist[i]] = i*round(0.2*dt[3])
    tclvalue(posvNew[[sliceslist[i]]]) = i*round(0.2*dt[3])  
    tkrreplot(img[[sliceslist[i]]])
    tkconfigure(label2[[sliceslist[i]]],text=posNew[sliceslist[i]])
  }
# return the window object to the master
  tt
}
# s <- list() OUTCOMMENTING THIS COULD HAVE LEAD TO MISTAKES
# show a slice of pvalues with anatomical overlay!
# should this really use adimpro???
show.slice <- function(x, anatomic, maxpvalue = 0.05, slice = 1, view = "axial", col.u, col.o, zlim.u =
                    NULL, zlim.o = NULL) {
  pvalue <- x$pvalue
  pvalue[pvalue>maxpvalue] <- 1
  pvalue[pvalue == 0] <- min(pvalue[pvalue>0])
  pvalue <- -log(pvalue)
  mask <- pvalue > 0
  ind2pos.ana <- conv.ip(anatomic, what="i2p")
  pos2ind.ana <- conv.ip(anatomic, what="p2i")
  ind2pos.func <- conv.ip(x, what="i2p")
  pos2ind.func <- conv.ip(x, what="p2i")
  pixdim.ana <- pixdim(anatomic$header,anatomic$format)
  pixdim.func <- pixdim(x$header,x$format)
  ttt.ana <- extract.data(anatomic)
  ddim.ana <- dim(ttt.ana) <- dim(ttt.ana)[1:3]
  if (view == "axial") {
    dfunc <- dim(pvalue)[1:2]
    if ((slice >= 1) & (slice <= dim(pvalue)[3])) {
      imgdata.o <- pvalue[,,slice]
      mask <- mask[,,slice]
    } else {
      mask <- imgdata.o <- array(0,dim=dfunc)
    }
    scale <- ceiling(max(abs(pixdim.func[1:2]))/min(abs(pixdim.ana)))
  } else if (view == "coronal") {
    dfunc <- dim(pvalue)[c(1,3)]
    if ((slice >= 1) & (slice <= dim(pvalue)[2])) {
      imgdata.o <- pvalue[,slice,]
      mask <- mask[,slice,]
    } else {
      mask <- imgdata.o <- array(0,dim=dfunc)
    }
    scale <- ceiling(max(abs(pixdim.func[c(1,3)]))/min(abs(pixdim.ana)))
  } else if (view == "sagittal") {
    dfunc <- dim(pvalue)[c(2,3)]
    if ((slice >= 1) & (slice <= dim(pvalue)[1])) {
      imgdata.o <- pvalue[slice,,]
      mask <- mask[slice,,]
    } else {
      mask <- imgdata.o <- array(0,dim=dfunc)
    }
    scale <- ceiling(max(abs(pixdim.func[2:3]))/min(abs(pixdim.ana)))
  } else {
    stop("unknown view",view)
  }
  imgdata.n <- array(0,dim=c(scale*dim(imgdata.o)))
  mask.n <- array(FALSE,dim=c(scale*dim(imgdata.o)))
  for (i in 1:dim(imgdata.o)[1]) {
    for (j in 1:dim(imgdata.o)[2]) {
      imgdata.n[(i-1)*scale+c(1:scale),(j-1)*scale+c(1:scale)] <- imgdata.o[i,j]
      mask.n[(i-1)*scale+c(1:scale),(j-1)*scale+c(1:scale)] <- imgdata.o[i,j]
    }
  }
  imgdata.o <- imgdata.n
  mask <- mask.n
  imgdata.u <- array(0, dim=dfunc*scale)
  for (i in 1:(dfunc[1]*scale)) {
    for (j in 1:(dfunc[2]*scale)) {
      if (view == "axial") {
        pos <- ind2pos.func( c(x$roixa+(2*i-1)/(2*scale)-0.5,
                            x$roiya+(2*j-1)/(2*scale)-0.5, x$roiza + slice - 1) )
      } else if (view == "coronal") {
        pos <- ind2pos.func( c(x$roixa+(2*i-1)/(2*scale)-0.5, 
                            x$roiya + slice - 1, x$roiza+(2*j-1)/(2*scale)-0.5) )
      } else if (view == "sagittal") {
        pos <- ind2pos.func( c(x$roixa + slice -1, 
                   x$roiya+(2*i-1)/(2*scale)-0.5, x$roiza+(2*j-1)/(2*scale)-0.5) )
      }
      ind.ana <- pos2ind.ana(pos)# this is real(!) index for anatomic image
      ii <- ind.ana[1]
      jj <- ind.ana[2]
      kk <- ind.ana[3]
      iint <- ceiling(ind.ana[1])# these are the integer indices
      jint <- ceiling(ind.ana[2])
      kint <- ceiling(ind.ana[3])
      if ((iint >= 1) & (jint >= 1) & (kint >= 1) &
          (iint <= ddim.ana[1]) & (jint <= ddim.ana[2]) & (kint <= ddim.ana[3])) {
        imgdata.u[i,j] <- ttt.ana[iint,jint,kint] * (ii - iint + 1) * 
                                        (jj - jint + 1) * (kk - kint + 1)
        if (kint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint,kint-1] * 
                            (ii - iint + 1) * (jj - jint + 1) * (kint - kk)
        if (jint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint-1,kint] * 
                             (ii - iint + 1) * (jint - jj) * (kk - kint + 1)
        if (iint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint-1,jint,kint] *
                             (iint - ii) * (jj - jint + 1) * (kk - kint + 1)
        if ((iint > 1) & (jint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] +
           ttt.ana[iint-1,jint-1,kint] * (iint - ii) * (jint - jj) * (kk - kint + 1)
        if ((iint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] +
           ttt.ana[iint-1,jint,kint-1] * (iint - ii) * (jj - jint + 1) * (kint - kk)
        if ((jint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] + 
           ttt.ana[iint,jint-1,kint-1] * (ii - iint + 1) * (jint - jj) * (kint - kk)
        if ((iint > 1) & (jint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] +  
           ttt.ana[iint-1,jint-1,kint-1] * (iint - ii) * (jint - jj) * (kint - kk)
      }
    }
  }
  if (is.null(zlim.o)) {
    zlim.o <- range(imgdata.o)
  } else {
    if (length(zlim.o) != 2) stop("zlim.o not length 2")
    if (zlim.o[2] < zlim.o[1]) stop("zlim.o[2] < zlim.o[1]")
    imgdata.o[imgdata.o > zlim.o[2]] <- zlim.o[2]
    imgdata.o[imgdata.o < zlim.o[1]] <- zlim.o[1]
  }
    if (is.null(zlim.u)) {
    zlim.u <- range(imgdata.u)
  } else {
    if (length(zlim.u) != 2) stop("zlim.u not length 2")
    if (zlim.u[2] < zlim.u[1]) stop("zlim.u[2] < zlim.u[1]")
    imgdata.u[imgdata.u > zlim.u[2]] <- zlim.u[2]
    imgdata.u[imgdata.u < zlim.u[1]] <- zlim.u[1]
  }
  img <- array(0, dim=c(dim(imgdata.u),3))
  for (i in 1:dim(imgdata.u)[1]) {
    for (j in 1:dim(imgdata.u)[2]) {
      if (mask[i,j]) {# use overlay
        level <- length(col.o) * (imgdata.o[i,j] - zlim.o[1]) / diff(zlim.o)
        level <- ceiling(level)# now in 0:length(col.o)
        if (is.na(level)) level <- 1
        if (level == 0) level <- 1# now in 1:length(col.o)
        img[i,j,] <- as.integer(col2rgb(col.o[level])) * 256
      } else {# use underlay
        level <- length(col.u) * (imgdata.u[i,j] - zlim.u[1]) / diff(zlim.u)
        level <- ceiling(level)# now in 0:length(col.u)
        if (is.na(level)) level <- 1
        if (level == 0) level <- 1# now in 1:length(col.u)
        img[i,j,] <- as.integer(col2rgb(col.u[level])) * 256
      }
    }
  }
  invisible(img)
}

pixdim <- function(header,format) {
  if (format == "NIFTI") {
    return(header$pixdim[2:4])
  } else if (format == "ANALYZE") {
    return(header$pixdim[2:4])
  } else if (format == "HEAD/BRIK") {
    return(header$DELTA)
  } else {
    stop("Not implemented for this data format:", format)
  }
}

conv.ip <- function(data, what="i2p") {
  if (!("fmridata" %in% class(data))) 
       stop("Cannot evaluate real-space position for this dataset. Not type fmridata!")
  if (data$format == "NIFTI") {
    if (data$header$qform > 0) {
      origin <- c(data$header$qoffsetx, data$header$qoffsety, data$header$qoffsetz)
      b <- data$header$quaternb
      c <- data$header$quaternc
      d <- data$header$quaternd
      a <- sqrt(pmax(0,1-b*b-c*c-d*d))
      R <- t(matrix(c(a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c,
                    2*b*c+2*a*d, a*a+c*c-b*b-d*d, 2*c*d -2*a*b,
                    2*b*d-2*a*c, 2*c*d+2*a*b, a*a+d*d-c*c-b*b),3,3))
      pixdim <- data$header$pixdim[2:4]
      qfac <- data$header$pixdim[1]
      if (what == "i2p") {
        return(function(ind) R %*% (c(1,1,qfac) * pixdim * (ind-1)) + origin)
      } else {
        return(function(pos) (solve(R) %*% (pos - origin))/(c(1,1,qfac) * pixdim) + 1)
      }
    } else if (data$header$sform > 0) {
      origin <- c(data$header$srowx[4],data$header$srowy[4],data$header$srowz[4])
      SR <- matrix(c(data$header$srowx[1],data$header$srowy[1],data$header$srowz[1],
                     data$header$srowx[2],data$header$srowy[2],data$header$srowz[2],
                     data$header$srowx[3],data$header$srowy[3],data$header$srowz[3]),3,3)
      if (what == "i2p") {
        return(function(ind) SR %*% (ind-1) + origin)
      } else {
        return(function(pos) solve(SR) %*% (pos - origin) + 1)
      }
    } else if (data$header$qform == 0) {
      warning("This method is specified only for compatibility reasons to
      ANALYZE 7.5. May not deliver useful results")
      pixdim <- data$header$pixdim[2:4]
      if (what == "i2p") {
        return(function(ind) pixdim * (ind-1))
      } else {
        return(function(pos) pos/pixdim + 1)
      }
    } else {
      stop("Neither Method 1, 2, nor 3 for real-space position applicable. 
      See NIFTI specification!")
    }
  } else if (data$format == "ANALYZE") {
    stop("Not yet implemented real-space position evaluation for this data format:",
             data$format)
  } else if (data$format == "HEAD/BRIK") {
    orientation <- data$header$ORIENT_SPECIFIC
    if (any(sort((orientation)%/%2) != 0:2)) 
          stop("invalid orientation",orientation,"found! \n")
    rxyz <- (orientation)%/%2+1
    xyz <- rxyz[rxyz]
    pixdim <- data$header$DELTA[xyz]
    origin <- data$header$ORIGIN[xyz]
    if (what == "i2p") {
      return(function(ind) pixdim * (ind[xyz]-1) + origin)
    } else {
      return(function(pos) ((pos-origin)/pixdim + 1)[rxyz])
    }
  } else if (data$format == "DICOM") {
    stop("Not yet implemented real-space position evaluation for this data format:",
          data$format) 
  } else {
    stop("Not implemented real-space position evaluation for this data format 
          (not in fmri package):", data$format)
  }
}

show.segmentslice <- function(x, 
		                       anatomic, 
							   slice = 1, 
							   view = "axial", 
							   col.u, 
							   col.o, 
							   zlim.u, 
							   zlim.o) {

    if (length(col.o) != 128) stop( "length of overlay color scale has to be 128!")
	if (length(col.u) != 128) stop( "length of underlay color scale has to be 128!")
	
	## we want to display the estimated parameter values					   
    value <- x$cbeta
	value[ x$segm == 0] <- NA ## no significant voxel

	## create the index2position (and back) conversion functions
	ind2pos.ana <- conv.ip( anatomic, what = "i2p")
	pos2ind.ana <- conv.ip( anatomic, what = "p2i")
	ind2pos.func <- conv.ip( x, what = "i2p")
	pos2ind.func <- conv.ip( x, what = "p2i")
	
	## get the voxel extensions
	pixdim.ana <- pixdim( anatomic$header, anatomic$format)
	pixdim.func <- pixdim( x$header, x$format)

	## extract the anatomical data from compressed object
	ttt.ana <- extract.data( anatomic)
	ddim.ana <- dim(ttt.ana) <- dim(ttt.ana)[1:3]
	
	## select correct overlay slice according to view
	if (view == "axial") {
		dim.o <- dim(value)[ 1:2]
		if ((slice >= 1) & (slice <= dim(value)[3])) {
			imgdata.o <- value[ , , slice]
		} else {
			imgdata.o <- array( NA, dim = dim.o)
		}
		scale <- ceiling( max( abs( pixdim.func[ 1:2]))/min( abs( pixdim.ana)))
		aspect <- pixdim.func[2]/pixdim.func[1]
	} else if (view == "coronal") {
		dim.o <- dim(value)[ c( 1, 3)]
		if ((slice >= 1) & (slice <= dim(value)[2])) {
			imgdata.o <- value[ , slice, ]
		} else {
			imgdata.o <- array( NA,dim = dim.o)
		}
		scale <- ceiling( max( abs( pixdim.func[ c( 1, 3)]))/min( abs( pixdim.ana)))
		aspect <- pixdim.func[3]/pixdim.func[1]
	} else if (view == "sagittal") {
		dim.o <- dim(value)[ c( 2, 3)]
		if ((slice >= 1) & (slice <= dim(value)[1])) {
			imgdata.o <- value[ slice, , ]
		} else {
			imgdata.o <- array(NA, dim = dim.o)
		}
		scale <- ceiling( max( abs( pixdim.func[ 2:3]))/min( abs( pixdim.ana)))
		aspect <- pixdim.func[3]/pixdim.func[2]
	} 
	
	## upscale overlay data slice
	imgdata.n <- array( 0, dim = c( scale*dim(imgdata.o)))
	for (i in 1:dim(imgdata.o)[1]) {
		for (j in 1:dim(imgdata.o)[2]) {
			imgdata.n[ (i-1)*scale + c(1:scale), (j-1)*scale + c(1:scale)] <- imgdata.o[ i, j]
		}
	}
	imgdata.o <- imgdata.n

	## create an upscaled underlay for this slice
	imgdata.u <- array( 0, dim = dim.o*scale)
	for (i in 1:(dim.o[1]*scale)) {
		for (j in 1:(dim.o[2]*scale)) {
			if (view == "axial") {
				pos <- ind2pos.func( c(x$roixa+(2*i-1)/(2*scale)-0.5,
								x$roiya+(2*j-1)/(2*scale)-0.5, x$roiza + slice - 1) )
			} else if (view == "coronal") {
				pos <- ind2pos.func( c(x$roixa+(2*i-1)/(2*scale)-0.5, 
								x$roiya + slice - 1, x$roiza+(2*j-1)/(2*scale)-0.5) )
			} else if (view == "sagittal") {
				pos <- ind2pos.func( c(x$roixa + slice -1, 
								x$roiya+(2*i-1)/(2*scale)-0.5, x$roiza+(2*j-1)/(2*scale)-0.5) )
			}
			ind.ana <- pos2ind.ana(pos)# this is real(!) index for anatomic image
			ii <- ind.ana[1]
			jj <- ind.ana[2]
			kk <- ind.ana[3]
			iint <- ceiling(ind.ana[1])# these are the integer indices
			jint <- ceiling(ind.ana[2])
			kint <- ceiling(ind.ana[3])
			if ((iint >= 1) & (jint >= 1) & (kint >= 1) &
					(iint <= ddim.ana[1]) & (jint <= ddim.ana[2]) & (kint <= ddim.ana[3])) {
				imgdata.u[i,j] <- ttt.ana[iint,jint,kint] * (ii - iint + 1) * 
						(jj - jint + 1) * (kk - kint + 1)
				if (kint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint,kint-1] * 
							(ii - iint + 1) * (jj - jint + 1) * (kint - kk)
				if (jint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint-1,kint] * 
							(ii - iint + 1) * (jint - jj) * (kk - kint + 1)
				if (iint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint-1,jint,kint] *
							(iint - ii) * (jj - jint + 1) * (kk - kint + 1)
				if ((iint > 1) & (jint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] +
							ttt.ana[iint-1,jint-1,kint] * (iint - ii) * (jint - jj) * (kk - kint + 1)
				if ((iint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] +
							ttt.ana[iint-1,jint,kint-1] * (iint - ii) * (jj - jint + 1) * (kint - kk)
				if ((jint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] + 
							ttt.ana[iint,jint-1,kint-1] * (ii - iint + 1) * (jint - jj) * (kint - kk)
				if ((iint > 1) & (jint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] +  
							ttt.ana[iint-1,jint-1,kint-1] * (iint - ii) * (jint - jj) * (kint - kk)
			}
		}
	}
	
	## user defined data limits to scale the image contrast
	## not sure whether this is what the user wants
	if (any(!is.na(imgdata.o))) {
		if (is.null(zlim.o)) {
			zlim.o <- range( abs(imgdata.o), na.rm = TRUE)
		} else {
			if (length(zlim.o) != 2) stop("zlim.o not length 2")
			if (zlim.o[2] < zlim.o[1]) stop("zlim.o[2] < zlim.o[1]")
			imgdata.o[imgdata.o > zlim.o[2]] <- zlim.o[2]
			imgdata.o[imgdata.o < zlim.o[1]] <- zlim.o[1]
			imgdata.o[imgdata.o < -zlim.o[2]] <- -zlim.o[2]
			imgdata.o[imgdata.o > -zlim.o[1]] <- -zlim.o[1]
		}
	}
	if (is.null(zlim.u)) {
		zlim.u <- range(imgdata.u, na.rm = TRUE)
	} else {
		if (length(zlim.u) != 2) stop("zlim.u not length 2")
		if (zlim.u[2] < zlim.u[1]) stop("zlim.u[2] < zlim.u[1]")
		imgdata.u[imgdata.u > zlim.u[2]] <- zlim.u[2]
		imgdata.u[imgdata.u < zlim.u[1]] <- zlim.u[1]
	}
	
	## create the break points for the color scale
	if (any(!is.na(imgdata.o))) {
		zlim.o <- quantile( abs(imgdata.o), c( 0, 0.9, 1), na.rm = TRUE)
		breaks.o <- c( -zlim.o[3], seq( -zlim.o[2], -zlim.o[1], length = 63), 0, seq( zlim.o[1], zlim.o[2], length = 63), zlim.o[3])
	}
	breaks.u <- seq( zlim.u[1], zlim.u[2], length = 129)
	
	## plot the image 
	graphics::image(1:dim(imgdata.u)[1], 1:dim(imgdata.u)[2], imgdata.u, col = col.u, asp = aspect, axes = FALSE, xlab = "",	ylab = "")
	if (any(!is.na(imgdata.o))) {
		graphics::image(1:dim(imgdata.o)[1], 1:dim(imgdata.o)[2], imgdata.o, asp = aspect, col = col.o, breaks = breaks.o, add = TRUE)
	}
	
	## finally create img for adimpro
	img <- array(0, dim = c( dim(imgdata.u), 3))
	for (i in 1:dim(imgdata.u)[1]) {
		for (j in 1:dim(imgdata.u)[2]) {
			if (!is.na(imgdata.o[ i, j])) { # use overlay
				ind <- (0:128)[imgdata.o[ i, j] < breaks.o]
				level <- ifelse(length(ind) == 0, 128, min(ind))
				img[ i, j, ] <- as.integer( col2rgb( col.o[level])) * 256
			} else { # use underlay
				ind <- (0:128)[imgdata.u[ i, j] < breaks.u]
				level <- ifelse(length(ind) == 0, 128, min(ind))
				img[ i, j, ] <- as.integer( col2rgb( col.u[level])) * 256
			}
		}
	}
	
	img
}

