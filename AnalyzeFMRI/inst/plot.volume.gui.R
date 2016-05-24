  require(tcltk) || stop("tcltk support is absent")
  require("tkrplot") || stop("tkrplot support is absent")


local({


  orientation <- function(L) {
    # Pour déterminer si on est en Radiological ou en Neurological
    # Inspiré de FSL
    # Voir : http://www.fmrib.ox.ac.uk/fslfaq/#general_lr
    
    sform.code <- L$sform.code
    qform.code <- L$qform.code
    
    if (is.null(sform.code)) { # Analyze format
      
      order <- -1*sign(L$pixdim[2])  # Radiological 
      
    } else { # NIFTI format
      
      sform.mat <- rbind(L$srow.x,L$srow.y,L$srow.z,c(0,0,0,1))
      qform.mat <- Q2R(c(L$quatern.b,L$quatern.c,L$quatern.d),if (L$pixdim[1]==0) 1 else L$pixdim[1])
      
      determinant <- -1
      
      if (sform.code != 0) {
        
        determinant <- det(sform.mat)
        
      } else if (qform.code != 0) {
        
        determinant <- det(qform.mat)
      }
      
      if (determinant<0) order <- -1 # Radiological 
      else order <- 1 # Neurological
      
      
    }
    
    return(order)
    
  }


maColorBar <- function (x, horizontal = TRUE, col = heat.colors(50), scale = 1:length(x),
    k = 10, ...)
{
    if (is.numeric(x)) {
        x <- x
        colmap <- col
    }
    else {
        colmap <- x
        low <- range(scale)[1]
        high <- range(scale)[2]
        x <- seq(low, high, length = length(x))
    }
    if (length(x) > k)
        x.small <- seq(x[1], x[length(x)], length = k)
    else x.small <- x
    if (horizontal) {
        image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "",
            ylab = "", col = colmap, ...)
        axis(1, at = rev(x.small), labels = signif(rev(x.small),
            2), srt = 270)
    }
    if (!horizontal) {
        image(1, x, matrix(x, 1, length(x)), axes = FALSE, xlab = "",
            ylab = "", col = colmap, ...)
        par(las = 1)
        axis(4, at = rev(x.small), labels = signif(rev(x.small),
            2))
        par(las = 0)
    }
    box()
}



    gui.file.fonc <- function(){
        tclvalue(file.name.fonc) <- tcl("tk_getOpenFile")   }
    
    
    gui.file.anat <- function(){
        tclvalue(file.name.anat) <- tcl("tk_getOpenFile")  }    
    
    gui.file.time.series <- function(){
        tclvalue(file.name.time.series) <- tcl("tk_getOpenFile")  }    



plot.volume <- function(vol.fonc="",vol.anat="",time.series="") {

## !!! J'ai pris method=3 (lecture de la matrice affine) pour faire conincider l'anat et la fonctionnelle
  ### si sform.code=0 cela ne marchera pas ...

hscaletmp <- as.numeric(tclvalue(hscaletmp))
vscaletmp <- as.numeric(tclvalue(vscaletmp))
  
# Permet d'afficher un volume cérébral anatomique et/ou fonctionnel
# Paramètres d'entrée:
#	vol.fonc: nom du fichier (.img ou .nii) ou bien une array 4D contenant les images volumiques fonctionnelles
#	vol.anat: nom du fichier (.img ou .nii) ou bien une array 4D contenant les images volumiques anatomiques
#       time.series : fichier .dat contenant les composantes temporelles extraites par l'ICA  


  # J'ai rajouté cela à la fonction plot.volume pour qu'elle soit appelable par le tkwidget
  vol.fonc <- tclvalue(file.name.fonc)
  vol.anat <- tclvalue(file.name.anat)
  time.series <- tclvalue(file.name.time.series)

  is.vol.fonc <- FALSE
  is.vol.anat <- FALSE
  is.time.series <- FALSE
  
  if (nzchar(vol.fonc)[1]) is.vol.fonc <- TRUE
  if (nzchar(vol.anat)[1]) is.vol.anat <- TRUE
  if (nzchar(time.series)[1]) is.time.series <- TRUE
  
  my.title <- "FMRI visualisation:"
  
  if (is.vol.fonc) {

    if (exists(vol.fonc)) vol.fonc <- get(vol.fonc)

    
    col.fonc <- heat.colors(256)
    
    
    if (is.character(vol.fonc)) { # reading of the image
      my.title <- paste(my.title,vol.fonc)
      hdr.fonc <- f.read.header(vol.fonc)
      vol.fonc <- f.read.volume(vol.fonc)

      # if flip == - 1 :Radiological else if flip == 1: Neurological
      flip <- orientation(hdr.fonc)
            
    } else {
      if (!exists("flip")) print("You must create variable flip (-1 for Radiological and 1 for Neurological)")
      if (length(dim(vol.fonc)) == 3) {
        vol.fonc <- array(vol.fonc,dim=c(dim(vol.fonc),1))
      } else {vol.fonc <- vol.fonc}
    }
    
    
    dimensions.fonc <- dim(vol.fonc)
    dim.fonc.sagit <- dimensions.fonc[1]
    dim.fonc.coron <- dimensions.fonc[2]
    dim.fonc.axia <- dimensions.fonc[3]
    dim.fonc.time <- dimensions.fonc[4]
      
    
# Variable initialisation
    nn.fonc.sagit <- ceiling(dim.fonc.sagit/2)
    nn.fonc.coron <- ceiling(dim.fonc.coron/2)
    nn.fonc.axia <- ceiling(dim.fonc.axia/2)
    nn.fonc.time <- 1
    
    mini.fonc <- min(vol.fonc[,,,nn.fonc.time])
    maxi.fonc <- max(vol.fonc[,,,nn.fonc.time])
    if ((mini.fonc-maxi.fonc) == 0) maxi.fonc <- maxi.fonc+0.00001
    breaks.fonc <- c(seq(from=mini.fonc,to=0,len=length(col.fonc)/2),seq(from=0,to=maxi.fonc,len=length(col.fonc)/2+1)) 
    
    SliderSagit.fonc <- tclVar(dim.fonc.sagit)
    SliderCoron.fonc <- tclVar(dim.fonc.coron)
    SliderAxia.fonc <- tclVar(dim.fonc.axia)
    SliderTime.fonc <- tclVar(dim.fonc.time)

  }


    if (is.vol.anat) {

      if (exists(vol.anat)) vol.anat <- get(vol.anat)
      
      col.anat <- gray(seq(from=0.2,to=1,len=256))

    
    if (is.character(vol.anat)) { # reading of the image
      my.title <- paste(my.title,vol.anat)
      hdr.anat <- f.read.header(vol.anat)
      vol.anat <- f.read.volume(vol.anat)

      # if flip == - 1 :Radiological else if flip == 1: Neurological
      flip <- orientation(hdr.anat)

    }
      if (!exists("flip")) print("You must create variable flip (-1 for Radiological and 1 for Neurological)")

    if (length(dim(vol.anat)) == 4) vol.anat <- vol.anat[,,,1]
    if (length(dim(vol.anat)) == 2) dim(vol.anat) <- c(dim(vol.anat),1)
    

    dimensions.anat <- dim(vol.anat)
    dim.anat.sagit <- dimensions.anat[1]
    dim.anat.coron <- dimensions.anat[2]
    dim.anat.axia <- dimensions.anat[3]
    
    
    mini.anat <- min(vol.anat[,,])
    maxi.anat <- max(vol.anat[,,])
    if ((mini.anat-maxi.anat) == 0) maxi.anat <- maxi.anat+0.00001
    breaks.anat <- seq(from=mini.anat,to=maxi.anat,len=length(col.anat)+1)
  
  
# Variable initialisation
    nn.anat.sagit <- ceiling(dim.anat.sagit/2)
    nn.anat.coron <- ceiling(dim.anat.coron/2)
    nn.anat.axia <- ceiling(dim.anat.axia/2)

    SliderSagit.anat <- tclVar(dim.anat.sagit)
    SliderCoron.anat <- tclVar(dim.anat.coron)
    SliderAxia.anat <- tclVar(dim.anat.axia)


  }

    
  if (is.vol.fonc) {
    
# Useful functions
    f.fonc.sagit <- function(...) {
      n.fonc.sagit <- as.numeric(tclvalue("nn.fonc.sagit"))
      n.fonc.time <- as.numeric(tclvalue("nn.fonc.time"))
      if (n.fonc.sagit != nn.fonc.sagit) {
        nn.fonc.sagit <<- n.fonc.sagit
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
      }
      if (n.fonc.time != nn.fonc.time) {
        nn.fonc.time <<- n.fonc.time
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
      }
    }
    
    f.fonc.coron <- function(...) {
      n.fonc.coron <- as.numeric(tclvalue("nn.fonc.coron"))
      n.fonc.time <- as.numeric(tclvalue("nn.fonc.time"))
      if (n.fonc.coron != nn.fonc.coron) {
        nn.fonc.coron <<- n.fonc.coron
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
      }
      if (n.fonc.time != nn.fonc.time) {
        nn.fonc.time <<- n.fonc.time
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
      }
    }
    
    
    f.fonc.axia <- function(...) {
      n.fonc.axia <- as.numeric(tclvalue("nn.fonc.axia"))
      n.fonc.time <- as.numeric(tclvalue("nn.fonc.time"))
      if (n.fonc.axia != nn.fonc.axia) {
        nn.fonc.axia <<- n.fonc.axia
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
      }
      if (n.fonc.time != nn.fonc.time) {
        nn.fonc.time <<- n.fonc.time
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
      }
    }
    
    
    f.fonc.time <- function(...) {
      n.fonc.time <- as.numeric(tclvalue("nn.fonc.time"))
      if (n.fonc.time != nn.fonc.time) {
        nn.fonc.time <<- n.fonc.time
        mini.fonc <<- min(vol.fonc[,,,nn.fonc.time])
        maxi.fonc <<- max(vol.fonc[,,,nn.fonc.time])
        if ((mini.fonc-maxi.fonc) == 0) maxi.fonc <<- maxi.fonc+0.00001
        breaks.fonc <<- c(seq(from=mini.fonc,to=0,len=length(col.fonc)/2),seq(from=0,to=maxi.fonc,len=length(col.fonc)/2+1)) 
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)
        if (is.time.series) tkrreplot(img.time.series)
      }
    }
    
    
    OnOK.fonc <- function() {tkdestroy(tt.fonc);if (is.time.series) tkdestroy(tt.fonc.time.series)}
    
    

  #set up base GUI window
    if(.Platform$OS.type == "windows") flush.console()
    
    tt.fonc <- tktoplevel(bg="#555555")

    if (is.vol.fonc) {
      tkwm.title(tt.fonc, my.title)}
    
    if (is.vol.anat) {
      tkwm.title(tt.fonc, my.title)}
    
    
    
  # frame 1 to contain images
    frame1.fonc <- tkframe(tt.fonc, relief = "groove", borderwidth = 0, bg = "#555555")
    label.fonc.sagit <- tklabel(frame1.fonc,text="Sagittal ", bg = "#aaaaaa")
    label.fonc.coron <- tklabel(frame1.fonc,text="Coronal ", bg = "#aaaaaa")
    slider.fonc.sagit <- tkscale(frame1.fonc, command=f.fonc.sagit, from=as.numeric(tclvalue(SliderSagit.fonc)), to=1, variable="nn.fonc.sagit",
                                 showvalue=TRUE, resolution=1, orient="verti") 
    img.fonc.sagit <- tkrplot(parent=frame1.fonc, fun=function() {
      par(mar=c(0,0,0,0), bg = "#555555")
    # coupe sagittale



      image(1:dim.fonc.coron,1:dim.fonc.axia,as.matrix(vol.fonc[nn.fonc.sagit,,,nn.fonc.time]),col=col.fonc,breaks=breaks.fonc,axes=FALSE,xlab="",ylab="",asp=2)
      abline(h=nn.fonc.axia,v=nn.fonc.coron,col="black")
    }
                              , hscale=hscaletmp, vscale=vscaletmp)
    tkbind(img.fonc.sagit, "<Button-1>", function(x, y) {
      asp <- 2
      wid.fonc <- as.integer(tkwinfo("width", img.fonc.sagit))
      hei.fonc <- as.integer(tkwinfo("height", img.fonc.sagit))      
      if (as.numeric(y)>(hei.fonc/2 - dim.fonc.axia*wid.fonc/dim.fonc.coron) & as.numeric(y)<(hei.fonc/2 + dim.fonc.axia*wid.fonc/dim.fonc.coron)) {
        nn.fonc.coron <<- ceiling(dim.fonc.coron * as.numeric(x) /wid.fonc)
        nn.fonc.axia <<- ceiling(dim.fonc.axia/2 + dim.fonc.coron*(hei.fonc/2-as.numeric(y))/(wid.fonc*asp))
        tkset(slider.fonc.coron,nn.fonc.coron)
        tkset(slider.fonc.axia,nn.fonc.axia)
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)

        if (is.vol.anat) {
          ijk <- c(nn.fonc.sagit,nn.fonc.coron,nn.fonc.axia)
          tmp <- round(xyz2ijk(ijk2xyz(ijk,method=3,hdr.fonc)$xyz,method=3,hdr.anat)$ijk)
          nn.anat.sagit <<- tmp[1]
          nn.anat.coron <<- tmp[2]
          nn.anat.axia <<- tmp[3]
          if (exists("img.anat.sagit")) {tkrreplot(img.anat.sagit);tkset(slider.anat.sagit,nn.anat.sagit)}
          if (exists("img.anat.coron")) {tkrreplot(img.anat.coron);tkset(slider.anat.coron,nn.anat.coron)}
          if (exists("img.anat.axia")) {tkrreplot(img.anat.axia);tkset(slider.anat.axia,nn.anat.axia)}
        }

      }
    })
    slider.fonc.coron <- tkscale(frame1.fonc, command=f.fonc.coron, from=as.numeric(tclvalue(SliderCoron.fonc)), to=1, variable="nn.fonc.coron",
                                 showvalue=TRUE, resolution=1, orient="verti")
    img.fonc.coron <- tkrplot(parent=frame1.fonc, fun=function() {
      par(mar=c(0,0,0,0), bg = "#555555")
    # coupe coronale
      image(1:dim.fonc.sagit,1:dim.fonc.axia,as.matrix(vol.fonc[,nn.fonc.coron,,nn.fonc.time]),col=col.fonc,breaks=breaks.fonc,axes=FALSE,xlab="",ylab="",asp=2)
      abline(h=nn.fonc.axia,v=nn.fonc.sagit,col="black")
    }
                              , hscale=hscaletmp, vscale=vscaletmp)
    tkbind(img.fonc.coron, "<Button-1>", function(x, y) {
      asp <- 2
      wid.fonc <- as.integer(tkwinfo("width", img.fonc.coron))
      hei.fonc <- as.integer(tkwinfo("height", img.fonc.coron))
      if (as.numeric(y)>(hei.fonc/2 - dim.fonc.axia*wid.fonc/dim.fonc.sagit) & as.numeric(y)<(hei.fonc/2 + dim.fonc.axia*wid.fonc/dim.fonc.sagit)) {  
        nn.fonc.sagit <<- ceiling(dim.fonc.sagit * as.numeric(x) /wid.fonc)
        nn.fonc.axia <<- ceiling(dim.fonc.axia/2 + dim.fonc.sagit*(hei.fonc/2-as.numeric(y))/(wid.fonc*asp))
        tkset(slider.fonc.sagit,nn.fonc.sagit)
        tkset(slider.fonc.axia,nn.fonc.axia)
        tkrreplot(img.fonc.sagit)
        tkrreplot(img.fonc.coron)
        tkrreplot(img.fonc.axia)
        tkrreplot(img.fonc.palette)

        if (is.vol.anat) {
          ijk <- c(nn.fonc.sagit,nn.fonc.coron,nn.fonc.axia)
          tmp <- round(xyz2ijk(ijk2xyz(ijk,method=3,hdr.fonc)$xyz,method=3,hdr.anat)$ijk)
          nn.anat.sagit <<- tmp[1]
          nn.anat.coron <<- tmp[2]
          nn.anat.axia <<- tmp[3]
          if (exists("img.anat.sagit")) {tkrreplot(img.anat.sagit);tkset(slider.anat.sagit,nn.anat.sagit)}
          if (exists("img.anat.coron")) {tkrreplot(img.anat.coron);tkset(slider.anat.coron,nn.anat.coron)}
          if (exists("img.anat.axia")) {tkrreplot(img.anat.axia);tkset(slider.anat.axia,nn.anat.axia)}
        }

      }
    })

  
  # frame 2 to contain Axial image and Palette
    frame2.fonc <- tkframe(tt.fonc, relief = "groove", borderwidth = 0, bg = "#555555")
    label.fonc.axia <- tklabel(frame2.fonc,text="Axial ", bg = "#aaaaaa")  
    label.fonc.palette <- tklabel(frame2.fonc,text="Palette values ", bg = "#aaaaaa")
    slider.fonc.axia <- tkscale(frame2.fonc, command=f.fonc.axia, from=as.numeric(tclvalue(SliderAxia.fonc)), to=1, variable="nn.fonc.axia",
                                showvalue=TRUE, resolution=1, orient="verti")
    img.fonc.axia <- tkrplot(parent=frame2.fonc, fun=function() {
      par(mar=c(0,0,0,0), bg = "#555555")
    # coupe axiale     
#      image(1:dim.fonc.sagit,1:dim.fonc.coron,as.matrix(vol.fonc[if (flip == 1) rev(1:dim.fonc.sagit) else 1:dim.fonc.sagit,,nn.fonc.axia,nn.fonc.time]),col=col.fonc,breaks=breaks.fonc,axes=FALSE,xlab="",ylab="",asp=1)
      image(1:dim.fonc.sagit,1:dim.fonc.coron,as.matrix(vol.fonc[1:dim.fonc.sagit,,nn.fonc.axia,nn.fonc.time]),col=col.fonc,breaks=breaks.fonc,axes=FALSE,xlab="",ylab="",asp=1)
      abline(h=nn.fonc.coron,v=nn.fonc.sagit,col="black")    
    }
                             , hscale=hscaletmp, vscale=vscaletmp)
    tkbind(img.fonc.axia, "<Button-1>", function(x, y) {
      wid.fonc <- as.integer(tkwinfo("width", img.fonc.axia))
      hei.fonc <- as.integer(tkwinfo("height", img.fonc.axia))
      nn.fonc.sagit <<- ceiling(dim.fonc.sagit * as.numeric(x) /wid.fonc)
      nn.fonc.coron <<- ceiling(dim.fonc.coron - dim.fonc.coron * as.numeric(y) /hei.fonc)
      tkset(slider.fonc.sagit,nn.fonc.sagit)
      tkset(slider.fonc.coron,nn.fonc.coron)
      tkrreplot(img.fonc.sagit)
      tkrreplot(img.fonc.coron)
      tkrreplot(img.fonc.axia)
      tkrreplot(img.fonc.palette)
 
      if (is.vol.anat) {
        ijk <- c(nn.fonc.sagit,nn.fonc.coron,nn.fonc.axia)
        tmp <- round(xyz2ijk(ijk2xyz(ijk,method=3,hdr.fonc)$xyz,method=3,hdr.anat)$ijk)
        nn.anat.sagit <<- tmp[1]
        nn.anat.coron <<- tmp[2]
        nn.anat.axia <<- tmp[3]
        if (exists("img.anat.sagit")) {tkrreplot(img.anat.sagit);tkset(slider.anat.sagit,nn.anat.sagit)}
        if (exists("img.anat.coron")) {tkrreplot(img.anat.coron);tkset(slider.anat.coron,nn.anat.coron)}
        if (exists("img.anat.axia")) {tkrreplot(img.anat.axia);tkset(slider.anat.axia,nn.anat.axia)}
      }

    })
    img.fonc.palette <- tkrplot(parent=frame2.fonc, fun=function() {
      par(mfrow=c(2,1), bg = "#555555")
      nf.fonc <- layout(matrix(c(1,2), nrow=2, ncol=1, byrow = TRUE),respect = FALSE,widths=c(0.8,0.8),height=c(0.3,0.8))
      par(mar=c(5.1,4.1,4.1,2.1), mai=c(0.5,0.82,0.2,0.42), bg = "#555555")
      maColorBar(x=seq(from=mini.fonc,to=maxi.fonc,len=length(col.fonc)+1),col=col.fonc,horizontal=TRUE)
      par(mar=c(1.1,2.1,2.1,0.1), mai=c(0.5,0.42,0.4,0.22))
      plot(vol.fonc[nn.fonc.sagit,nn.fonc.coron,nn.fonc.axia,],type="l",main=paste("Time course of the selected voxel.\nValue at that time: ",round(vol.fonc[nn.fonc.sagit,nn.fonc.coron,nn.fonc.axia,nn.fonc.time],4),sep=""),xlab="",ylab="",cex.main=0.8)
      points(nn.fonc.time,vol.fonc[nn.fonc.sagit,nn.fonc.coron,nn.fonc.axia,nn.fonc.time],col="blue",pch=1,cex=1.2)
    }
                                , hscale=hscaletmp, vscale=vscaletmp)
    
 
  # frame 3 to contain time slider and Quit button
    frame3.fonc <- tkframe(tt.fonc, relief = "groove", borderwidth = 2, bg = "#555555")
    neuro.fonc <- tklabel(frame3.fonc,text=if (flip == 1) "Neurological convention" else if (flip == -1) "Radiological convention")
    timenb.fonc <- tklabel(frame3.fonc,text="Time", bg = "#aaaaaa")
    slider.fonc.time <- tkscale(frame3.fonc, command=f.fonc.time, from=1, to=as.numeric(tclvalue(SliderTime.fonc)), variable="nn.fonc.time",
                                showvalue=TRUE, resolution=1,bigincrement=1,sliderlength=10,digits=0,length=400, orient="horiz")
    OK.but.fonc <- tkbutton(frame3.fonc,text="Quit",command=OnOK.fonc, bg = "#aaaaaa")
  
  
# We build the tk window with its widgets
    tkgrid(tklabel(frame1.fonc,text="             ", bg = "#555555"),label.fonc.sagit,tklabel(frame1.fonc,text="             ", bg = "#555555"),label.fonc.coron ,padx = 1, pady = 1)
    tkgrid(slider.fonc.sagit,img.fonc.sagit,slider.fonc.coron,img.fonc.coron ,padx = 1, pady = 1)
    tkgrid(frame1.fonc)
    tkgrid(tklabel(frame2.fonc,text="             ", bg = "#555555"),label.fonc.axia,tklabel(frame2.fonc,text="             ", bg = "#555555"),label.fonc.palette ,padx = 1, pady = 1)
    tkgrid(slider.fonc.axia,img.fonc.axia,tklabel(frame2.fonc,text="             ", bg = "#555555"), img.fonc.palette ,padx = 1, pady = 1)
    tkgrid(frame2.fonc)
    tkgrid(neuro.fonc,timenb.fonc,slider.fonc.time,OK.but.fonc, padx = 10, pady = 10)
    tkgrid(frame3.fonc, sticky = "ew")
    

    
    if (is.time.series) {

      if (exists(time.series)) time.series <- get(time.series)
      
      
      if (is.character(time.series)) { # reading of the image
        time.series <- as.matrix(read.table(time.series,header=TRUE))
      } else {
        time.series <- as.matrix(time.series)
      }
      
      tt.fonc.time.series <- tktoplevel(bg="#555555")
      tkwm.title(tt.fonc.time.series, "Associated time course")
      
      img.time.series <- tkrplot(parent=tt.fonc.time.series, fun=function() {
        par(mar=c(2,2,2,2), bg = "#555555")
        if (dim(time.series)[1] > dim(time.series)[2]) {
          plot(time.series[,nn.fonc.time],type="l",main=paste("Component number ",nn.fonc.time,sep=""))
        } else {
          plot(time.series[nn.fonc.time,],type="l",main=paste("Component number ",nn.fonc.time,sep=""))
        }
      }
                                 , hscale=1.5, vscale=0.5)                   
      
      tkgrid(img.time.series)
      
    }
    
  }


  if (is.vol.anat) {

# Useful functions
    f.anat.sagit <- function(...) {
      n.anat.sagit <- as.numeric(tclvalue("nn.anat.sagit"))
      if (n.anat.sagit != nn.anat.sagit) {
        nn.anat.sagit <<- n.anat.sagit
        tkrreplot(img.anat.sagit)
        tkrreplot(img.anat.coron)
        tkrreplot(img.anat.axia)
      }
    }
    
    f.anat.coron <- function(...) {
      n.anat.coron <- as.numeric(tclvalue("nn.anat.coron"))
      if (n.anat.coron != nn.anat.coron) {
        nn.anat.coron <<- n.anat.coron
        tkrreplot(img.anat.sagit)
        tkrreplot(img.anat.coron)
        tkrreplot(img.anat.axia)
      }
    }
    
    
    f.anat.axia <- function(...) {
      n.anat.axia <- as.numeric(tclvalue("nn.anat.axia"))
      if (n.anat.axia != nn.anat.axia) {
        nn.anat.axia <<- n.anat.axia
        tkrreplot(img.anat.sagit)
        tkrreplot(img.anat.coron)
        tkrreplot(img.anat.axia)
      }
    }
    
    
    OnOK.anat <- function() {tkdestroy(tt.anat)}
 

  #set up base GUI window
    if(.Platform$OS.type == "windows") flush.console()
    
    tt.anat <- tktoplevel(bg="#555555")
    tkwm.title(tt.anat, "MRI visualisation")


  # frame 1 to contain images
    frame1.anat <- tkframe(tt.anat, relief = "groove", borderwidth = 0, bg = "#555555")
    label.anat.sagit <- tklabel(frame1.anat,text="Sagittal ", bg = "#aaaaaa")
    label.anat.coron <- tklabel(frame1.anat,text="Coronal ", bg = "#aaaaaa")
    slider.anat.sagit <- tkscale(frame1.anat, command=f.anat.sagit, from=as.numeric(tclvalue(SliderSagit.anat)), to=1, variable="nn.anat.sagit",
                                 showvalue=TRUE, resolution=1, orient="verti") 
    img.anat.sagit <- tkrplot(parent=frame1.anat, fun=function() {
      par(mar=c(0,0,0,0), bg = "#555555")
    # coupe sagittale
      image((1:dim.anat.coron)/hdr.anat$pixdim[3],(1:dim.anat.axia)/hdr.anat$pixdim[4],as.matrix(vol.anat[nn.anat.sagit,,]),col=col.anat,breaks=breaks.anat,axes=FALSE,xlab="",ylab="",asp=dim.anat.coron/dim.anat.axia)
      abline(h=nn.anat.axia/hdr.anat$pixdim[4],v=nn.anat.coron/hdr.anat$pixdim[3],col="red")
    }
                              , hscale=hscaletmp, vscale=vscaletmp)
    tkbind(img.anat.sagit, "<Button-1>", function(x, y) {
      wid.anat <- as.integer(tkwinfo("width", img.anat.sagit))
      hei.anat <- as.integer(tkwinfo("height", img.anat.sagit))
      nn.anat.coron <<- ceiling(dim.anat.coron * as.numeric(x) /wid.anat)
      nn.anat.axia <<- ceiling( dim.anat.axia - dim.anat.axia * as.numeric(y) /hei.anat )
      tkset(slider.anat.coron,nn.anat.coron)	
      tkset(slider.anat.axia,nn.anat.axia)
      tkrreplot(img.anat.sagit)
      tkrreplot(img.anat.coron)
      tkrreplot(img.anat.axia)

      if (is.vol.fonc) {
        ijk <- c(nn.anat.sagit,nn.anat.coron,nn.anat.axia)
        tmp <- round(xyz2ijk(ijk2xyz(ijk,method=3,hdr.anat)$xyz,method=3,hdr.fonc)$ijk)
        nn.fonc.sagit <<- tmp[1]
        nn.fonc.coron <<- tmp[2]
        nn.fonc.axia <<- tmp[3]
        if ((1 <= nn.fonc.sagit) & (nn.fonc.sagit <= dim.fonc.sagit) &(1 <= nn.fonc.coron) & (nn.fonc.coron <= dim.fonc.coron) & (1 <= nn.fonc.axia) & (nn.fonc.axia <= dim.fonc.axia)) {
          tkrreplot(img.fonc.sagit)
          tkrreplot(img.fonc.coron)
          tkrreplot(img.fonc.axia)
          tkrreplot(img.fonc.palette)
          tkset(slider.fonc.sagit,nn.fonc.sagit)
          tkset(slider.fonc.coron,nn.fonc.coron)
          tkset(slider.fonc.axia,nn.fonc.axia)
        }
      }

      
    })
    slider.anat.coron <- tkscale(frame1.anat, command=f.anat.coron, from=as.numeric(tclvalue(SliderCoron.anat)), to=1, variable="nn.anat.coron",
                                 showvalue=TRUE, resolution=1, orient="verti")
    img.anat.coron <- tkrplot(parent=frame1.anat, fun=function() {
      par(mar=c(0,0,0,0), bg = "#555555")
    # coupe coronale
      image((1:dim.anat.sagit)/abs(hdr.anat$pixdim[2]),(1:dim.anat.axia)/hdr.anat$pixdim[4],as.matrix(vol.anat[,nn.anat.coron,]),col=col.anat,breaks=breaks.anat,axes=FALSE,xlab="",ylab="",asp=dim.anat.sagit/dim.anat.axia)
      abline(h=nn.anat.axia/hdr.anat$pixdim[4],v=nn.anat.sagit/abs(hdr.anat$pixdim[2]),col="red")
    }
                              , hscale=hscaletmp, vscale=vscaletmp)
    tkbind(img.anat.coron, "<Button-1>", function(x, y) {
      wid.anat <- as.integer(tkwinfo("width", img.anat.coron))
      hei.anat <- as.integer(tkwinfo("height", img.anat.coron))
      nn.anat.sagit <<- ceiling(dim.anat.sagit * as.numeric(x) /wid.anat)
      nn.anat.axia <<- ceiling( dim.anat.axia - dim.anat.axia * as.numeric(y) /hei.anat )
      tkset(slider.anat.sagit,nn.anat.sagit)
      tkset(slider.anat.axia,nn.anat.axia)
      tkrreplot(img.anat.sagit)
      tkrreplot(img.anat.coron)
      tkrreplot(img.anat.axia)

      if (is.vol.fonc) {
        ijk <- c(nn.anat.sagit,nn.anat.coron,nn.anat.axia)
        tmp <- round(xyz2ijk(ijk2xyz(ijk,method=3,hdr.anat)$xyz,method=3,hdr.fonc)$ijk)
        nn.fonc.sagit <<- tmp[1]
        nn.fonc.coron <<- tmp[2]
        nn.fonc.axia <<- tmp[3]
        if ((1 <= nn.fonc.sagit) & (nn.fonc.sagit <= dim.fonc.sagit) &(1 <= nn.fonc.coron) & (nn.fonc.coron <= dim.fonc.coron) & (1 <= nn.fonc.axia) & (nn.fonc.axia <= dim.fonc.axia)) {
          tkrreplot(img.fonc.sagit)
          tkrreplot(img.fonc.coron)
          tkrreplot(img.fonc.axia)
          tkrreplot(img.fonc.palette)
          tkset(slider.fonc.sagit,nn.fonc.sagit)
          tkset(slider.fonc.coron,nn.fonc.coron)
          tkset(slider.fonc.axia,nn.fonc.axia)
        }
      }

      
      
    })

  
  # frame 2 to contain Axial image and Palette
    frame2.anat <- tkframe(tt.anat, relief = "groove", borderwidth = 0, bg = "#555555")
    label.anat.axia <- tklabel(frame2.anat,text="Axial ", bg = "#aaaaaa")  
    label.anat.palette <- tklabel(frame2.anat,text="Palette values ", bg = "#aaaaaa")
    slider.anat.axia <- tkscale(frame2.anat, command=f.anat.axia, from=as.numeric(tclvalue(SliderAxia.anat)), to=1, variable="nn.anat.axia",
                                showvalue=TRUE, resolution=1, orient="verti")
    img.anat.axia <- tkrplot(parent=frame2.anat, fun=function() {
      par(mar=c(0,0,0,0), bg = "#555555")
    # coupe axiale     
      image((1:dim.anat.sagit)/abs(hdr.anat$pixdim[2]),(1:dim.anat.coron)/hdr.anat$pixdim[3],as.matrix(vol.anat[if (flip == 1) rev(1:dim.anat.sagit) else 1:dim.anat.sagit,,nn.anat.axia]),col=col.anat,breaks=breaks.anat,axes=FALSE,xlab="",ylab="",asp=dim.anat.sagit/dim.anat.coron)
      abline(h=nn.anat.coron/hdr.anat$pixdim[3],v=nn.anat.sagit/abs(hdr.anat$pixdim[2]),col="red")    
    }
                             , hscale=hscaletmp, vscale=vscaletmp)
    tkbind(img.anat.axia, "<Button-1>", function(x, y) {
      wid.anat <- as.integer(tkwinfo("width", img.anat.axia))
      hei.anat <- as.integer(tkwinfo("height", img.anat.axia))
      nn.anat.sagit <<- ceiling(dim.anat.sagit * as.numeric(x) /wid.anat)
      nn.anat.coron <<- ceiling(dim.anat.coron - dim.anat.coron * as.numeric(y) /hei.anat)
      tkset(slider.anat.sagit,nn.anat.sagit)
      tkset(slider.anat.coron,nn.anat.coron)
      tkrreplot(img.anat.sagit)
      tkrreplot(img.anat.coron)
      tkrreplot(img.anat.axia)

      if (is.vol.fonc) {
        ijk <- c(nn.anat.sagit,nn.anat.coron,nn.anat.axia)
        tmp <- round(xyz2ijk(ijk2xyz(ijk,method=3,hdr.anat)$xyz,method=3,hdr.fonc)$ijk)
        nn.fonc.sagit <<- tmp[1]
        nn.fonc.coron <<- tmp[2]
        nn.fonc.axia <<- tmp[3]
        if ((1 <= nn.fonc.sagit) & (nn.fonc.sagit <= dim.fonc.sagit) &(1 <= nn.fonc.coron) & (nn.fonc.coron <= dim.fonc.coron) & (1 <= nn.fonc.axia) & (nn.fonc.axia <= dim.fonc.axia)) {
          tkrreplot(img.fonc.sagit)
          tkrreplot(img.fonc.coron)
          tkrreplot(img.fonc.axia)
          tkrreplot(img.fonc.palette)
          tkset(slider.fonc.sagit,nn.fonc.sagit)
          tkset(slider.fonc.coron,nn.fonc.coron)
          tkset(slider.fonc.axia,nn.fonc.axia)
        }
      }

      
    })
    img.anat.palette <- tkrplot(parent=frame2.anat, fun=function() {
      nf.anat <- layout(matrix(c(0,1,0), 1, 3, byrow = TRUE),respect = FALSE,widths=c(4,1,4),height=0.8)
      par(mar=c(0.4,0.1,0.4,0.1), bg = "#555555")
      maColorBar(x=seq(from=mini.anat,to=maxi.anat,len=length(col.anat)+1),col=col.anat,horizontal=FALSE)
    }
                                , hscale=hscaletmp, vscale=vscaletmp)
  
 
  # frame 3 to contain Quit button
    frame3.anat <- tkframe(tt.anat, relief = "groove", borderwidth = 2, bg = "#555555") 
    neuro.anat <- tklabel(frame3.anat,text=if (flip == 1) "Neurological convention" else if (flip == -1) "Radiological convention")
    Quit.but.anat <- tkbutton(frame3.anat,text="Quit",command=OnOK.anat, bg = "#aaaaaa")
  
  
# We build the tk window with its widgets
    tkgrid(tklabel(frame1.anat,text="             ", bg = "#555555"),label.anat.sagit,tklabel(frame1.anat,text="             ", bg = "#555555"),label.anat.coron ,padx = 1, pady = 1)
    tkgrid(slider.anat.sagit,img.anat.sagit,slider.anat.coron,img.anat.coron ,padx = 1, pady = 1)
    tkgrid(frame1.anat)
    tkgrid(tklabel(frame2.anat,text="             ", bg = "#555555"),label.anat.axia,tklabel(frame2.anat,text="             ", bg = "#555555"),label.anat.palette ,padx = 1, pady = 1)
    tkgrid(slider.anat.axia,img.anat.axia,tklabel(frame2.anat,text="             ", bg = "#555555"), img.anat.palette ,padx = 1, pady = 1)
    tkgrid(frame2.anat)
    tkgrid(neuro.anat,Quit.but.anat, padx = 10, pady = 10)
    tkgrid(frame3.anat, sticky = "ew")


    
  }

  
}

 

    gui.end<-function(...){
        tkdestroy(base.plot)
    }






























    ##set tcl variables
    file.name.fonc <- tclVar()
    file.name.anat <- tclVar()
    file.name.time.series <- tclVar()

   
    #set up base GUI window
    if(.Platform$OS.type == "windows") flush.console()
    
    base.plot <- tktoplevel(bg="#555555")
    tkwm.title(base.plot, "(F)MRI visualisation")
    
    
    #frame to contain file selection
    frame1 <- tkframe(base.plot, relief = "groove", borderwidth = 2, bg = "#555555")
    
    file.fonc.entry <- tkentry(frame1, textvariable = file.name.fonc, width = 50, bg = "#ffffff")
    file.fonc.find.but <- tkbutton(frame1, text = "Select functionnal NIFTI file (or give R object name)", width = 50, command = gui.file.fonc, bg = "#aaaaaa", anchor = "c")
    tkgrid(file.fonc.find.but, file.fonc.entry, pady = 10, padx = 10)
    
    file.time.series.entry <- tkentry(frame1, textvariable = file.name.time.series, width = 50, bg = "#ffffff")
    file.time.series.find.but <- tkbutton(frame1, text = "Select time series component file (or give R object name)", width = 50, command = gui.file.time.series, bg = "#aaaaaa", anchor = "c")
    tkgrid(file.time.series.find.but, file.time.series.entry, pady = 10, padx = 10)

    file.anat.entry <- tkentry(frame1, textvariable = file.name.anat, width = 50, bg = "#ffffff")
    file.anat.find.but <- tkbutton(frame1, text = "Select anatomical NIFTI file (or give R object name)", width = 50, command = gui.file.anat, bg = "#aaaaaa", anchor = "c")
    tkgrid(file.anat.find.but, file.anat.entry, pady = 10, padx = 10)
    

    tkgrid(frame1)
    
    hscaletmp <- tclVar(0.6)
    vscaletmp <- tclVar(0.6)
  
    #frame for start and end buttons     
    frame2 <- tkframe(base.plot, borderwidth = 2, bg = "#555555")
    go.but <- tkbutton(frame2, text = "Start", bg = "#aaaaaa", command = plot.volume)
    q.but <- tkbutton(frame2, text = "Quit",command = gui.end, bg = "#aaaaaa")
 
hscale.txt <- tklabel(frame2,text="hscale factor", bg = "#aaaaaa")
    hscale.entry <- tkentry(frame2, textvariable = hscaletmp, width = 5, bg = "#ffffff")
vscale.txt <- tklabel(frame2,text="vscale factor", bg = "#aaaaaa")
    vscale.entry <- tkentry(frame2, textvariable = vscaletmp, width = 5, bg = "#ffffff")


   tkgrid(go.but, q.but, hscale.txt,hscale.entry, vscale.txt,vscale.entry, padx = 30, pady = 20)
    tkgrid(frame2)



















})
