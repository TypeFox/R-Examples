histnbmap<- function(sp.obj, nb.obj, longlat = NULL, nbcol=10, type = c("count","percent", "density"),
sup=FALSE, criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16, col="lightblue3",
xlab="", ylab="count", axes=FALSE, lablong="", lablat="")
{
# Verification of the Spatial Object sp.obj
class.obj<-class(sp.obj)[1]
if(substr(class.obj,1,7)!="Spatial") stop("sp.obj may be a Spatial object")
spdf<-(class.obj=="SpatialPolygonsDataFrame")

# we propose to refind the same arguments used in first version of GeoXp
coords<-coordinates(sp.obj)

if(substr(class.obj,nchar(class.obj)-8,nchar(class.obj))=="DataFrame")
{listvar<-sp.obj@data
listnomvar<-names(sp.obj@data)
}
else
{listvar<-NULL
listnomvar<-NULL
}

# Code which was necessary in the previous version
object<-nb.obj
# if(is.null(carte) & substr(class.obj,1,15)=="SpatialPolygons") carte<-spdf2list(sp.obj)$poly

 # for identifyng the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")

# initialisation
  nointer<-FALSE
  nocart<-FALSE
  buble<-FALSE
  labvar=c(xlab,ylab)
  z<-NULL
  legmap<-NULL
  legends<-list(FALSE,FALSE,"","")

# Transformation data.frame en matrix
  if((length(listvar)>0)&&(dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)

# Récupération des objets spdep 
   nb <- object
   W<-nb2mat(nb,zero.policy=TRUE)
   if (!inherits(nb, "nb"))
        stop("Not a neighbours list")
   c.nb <- card(nb)
   n.nb <- length(nb)
   regids <- attr(nb, "region.id")
   if (is.null(regids))
       regids <- as.character(1:n.nb)

   long <- coords[, 1]
   lat <- coords[, 2]
    
   obs <- matrix(FALSE, nrow=n.nb, ncol=n.nb)
    
   graf<-"Neighbourplot1"

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()


# this code comes from a spdep function 
    if (!is.null(coords)) 
    {
        if (inherits(coords, "SpatialPoints")) 
        {
            if ((is.null(longlat) || !is.logical(longlat)) &&
                !is.na(is.projected(coords)) && !is.projected(coords)) 
            {
            longlat <- TRUE
            }
            else longlat <- FALSE
            coords <- coordinates(coords)
        }
        else if (is.null(longlat) || !is.logical(longlat))
            longlat <- FALSE
        if (!is.matrix(coords))
            stop("Data not in matrix form")
        if (any(is.na(coords)))
            stop("Data include NAs")
        np <- nrow(coords)
        if (np != n.nb)
            stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
        # dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np),
        #    as.integer(dimension), as.integer(longlat), PACKAGE = "spdep")[[1]]
             dlist <- nbdists(nb=nb, coords=coords, longlat=longlat)
    }

# this code also ....
 if(sup==TRUE)
 {
  for (i in 1:n.nb)
   { nb[[i]]=as.integer(nb[[i]][which(dlist[[i]]==max(dlist[[i]]))])
   }
   
    if (!is.null(coords)) 
    {
        if (inherits(coords, "SpatialPoints")) 
        {
            if ((is.null(longlat) || !is.logical(longlat)) &&
                !is.na(is.projected(coords)) && !is.projected(coords)) 
            {
             longlat <- TRUE
            }
            else longlat <- FALSE
            coords <- coordinates(coords)
        }
        else if (is.null(longlat) || !is.logical(longlat))
            longlat <- FALSE
        if (!is.matrix(coords))
            stop("Data not in matrix form")
        if (any(is.na(coords)))
            stop("Data include NAs")
        np <- nrow(coords)
        if (np != n.nb)
            stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
       # dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np),
       #     as.integer(dimension), as.integer(longlat), PACKAGE = "spdep")[[1]]
        dlist <- nbdists(nb=nb, coords=coords, longlat=longlat)
    }   
   
    W<-nb2mat(nb)   
   
 }                                


####################################################
# sélection d'un point sur la carte
####################################################

pointfunc<-function()
{
   if (graf=="Neighbourplot2") SGfunc()

    quit <- FALSE
    graf<<-"Neighbourplot1"
   
    dev.set(2)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    while(!quit)
    {
      dev.set(2)
        if(spdf & nrow(sp.obj)>75 & !buble) 
         {points(long,lat,pch=16,col='royalblue')}
       
      loc<-locator(1)
        if (is.null(loc))
        {
          quit<-TRUE
          carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
          method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
          criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
          next
        }
        
       if(!spdf|nrow(sp.obj)>75)
       { 
       obs<<-selectmap(var1=long,var2=lat,obs=obs,Xpoly=loc[1], Ypoly=loc[2], method="point")}
       else
       {if(gContains(sp.obj,SpatialPoints(cbind(loc$x,loc$y),proj4string=CRS(proj4string(sp.obj)))))
        {for (i in 1:nrow(sp.obj))
          {if(gContains(sp.obj[i,],SpatialPoints(cbind(loc$x,loc$y),proj4string=CRS(proj4string(sp.obj)))))
           {obs[i,]<<-!obs[i,]
            break}  
          }
         } 
       }
        # graphiques
        carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  lablong=lablong, lablat=lablat, label=label, symbol=pch,
        method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
        criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

        graphique(var1=nb,var2=dlist, obs=obs, num=3, bin=type, graph="histo.nb",nbcol=nbcol, W=W,
        labvar=labvar, symbol=pch,couleurs=col)

    }
  }
  

####################################################
# sélection d'un polygone
####################################################

polyfunc<-function()
{
   if (graf=="Neighbourplot2") SGfunc()

    graf<<-"Neighbourplot1"
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    dev.set(2)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
      if(spdf) 
       {points(long,lat,pch=16,col='royalblue')}
         
    while(!quit)
    {
        dev.set(2)
        loc<-locator(1)
        if(is.null(loc))
        {
          quit<-TRUE
          next
        }

        polyX <- c(polyX, loc[1])
        polyY <- c(polyY, loc[2])
        lines(polyX,polyY)
    }

    polyX <- c(polyX, polyX[1])
    polyY <- c(polyY, polyY[1])
if (length(polyX)>0)
{
    lines(polyX,polyY)

    obs <<- selectmap(var1=long, var2=lat, obs=obs, Xpoly=polyX, Ypoly=polyY, method="poly")

    # graphiques
    carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
    method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
    criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
      
    graphique(var1=nb,var2=dlist, obs=obs, num=3,bin=type, graph="histo.nb",nbcol=nbcol, W=W,
    labvar=labvar, symbol=pch,couleurs=col); #   obs <<- matrix(FALSE, nrow=length(long), ncol=length(long));

}
  }




####################################################
# sélection d'une barre de l'histogramme
####################################################

barfunc<-function()
{
   if (graf=="Neighbourplot1") SGfunc()

    graf<<-"Neighbourplot2"
    quit <- FALSE
    dev.set(3)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    while(!quit)
    {
      dev.set(3)
      loc<-locator(1)
        if(is.null(loc))
        {
          quit<-TRUE
          graphique(var1=nb,var2=dlist, obs=obs, num=3, graph="histo.nb",nbcol=nbcol, W=W,
          labvar=labvar, symbol=pch,couleurs=col)
          next
        }
        obs<<-selectstat(var1=nb,var2=dlist,obs=obs,Xpoly=loc[1], Ypoly=loc[2],method="nbhist", nbcol=nbcol)
 

   carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
   method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
   criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
      
   graphique(var1=nb,var2=dlist, obs=obs, num=3,bin=type, graph="histo.nb",nbcol=nbcol, W=W,
   labvar=labvar, symbol=pch,couleurs=col)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    }
  }


####################################################
# contour des unités spatiales
####################################################

cartfunc <- function()
{ 
 if (length(carte) != 0)
   {
    ifelse(!nocart,nocart<<-TRUE,nocart<<-FALSE)
    carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
    method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
    criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
    }
   else
   {
    tkmessageBox(message="Spatial contours have not been given",icon="warning",type="ok")    
   }
}

####################################################
# Open a no interactive selection
####################################################

fnointer<-function() 
{
 if (length(criteria) != 0)
 {
   ifelse(!nointer,nointer<<-TRUE,nointer<<-FALSE)
    carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  lablong=lablong, lablat=lablat, label=label, symbol=pch,
    method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
    criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
 }
 else
 {
  tkmessageBox(message="Criteria has not been given",icon="warning",type="ok")
 }
 
}
  
####################################################
# Bubble
####################################################

fbubble<-function()
{
  res2<-choix.bubble(buble,listvar,listnomvar,legends)
  
  buble <<- res2$buble
  legends <<- res2$legends
  z <<- res2$z
  legmap <<- res2$legmap
  
  carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
  method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
  criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)

}


    
####################################################
# rafraichissement des graphiques
####################################################

SGfunc<-function()
{
    obs <<- matrix(FALSE, nrow=length(long), ncol=length(long));

    # graphiques
    carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
    method="Neighbourplot3", W=W,axis=axes,cex.lab=cex.lab,legmap=legmap,legends=legends,buble=buble,
    criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart)
      
    graphique(var1=nb,var2=dlist, obs=obs, num=3,bin=type, graph="histo.nb",nbcol=nbcol, W=W,
    labvar=labvar, symbol=pch,couleurs=col)
}
  
####################################################
# quitter l'application
####################################################

quitfunc<-function()
{
    assign("GeoXp.open", FALSE, envir = baseenv())
    tkdestroy(tt)
}


####################################################
# Graphique de base
####################################################
 # Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2)
 {
  carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
  method="Neighbourplot3", W=W,axis=axes,legmap=legmap,legends=legends,buble=buble,
  criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart,cex.lab=cex.lab)
     
  graphique(var1=nb,var2=dlist, obs=obs, num=3,bin=type, graph="histo.nb",nbcol=nbcol, W=W,
  labvar=labvar, symbol=pch,couleurs=col)
  assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, lablong=lablong, lablat=lablat, label=label, symbol=pch,
   method="Neighbourplot3", W=W,axis=axes,legmap=legmap,legends=legends,buble=buble,
   criteria=criteria,nointer=nointer,cbuble=z,carte=carte,nocart=nocart,cex.lab=cex.lab)
     
   graphique(var1=nb,var2=dlist, obs=obs, num=3,bin=type, graph="histo.nb",nbcol=nbcol, W=W,
   labvar=labvar, symbol=pch,couleurs=col)
   assign("GeoXp.open", TRUE, envir = baseenv())}
 }
}

####################################################
# Boîte de Dialogue
####################################################

if(interactive())
{
fontheading<-tkfont.create(family="times",size=14,weight="bold")

tt <- tktoplevel()
tkwm.title(tt, "histnbmap")

frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
foreground = "darkred", background = "white"))

point.but <- tkbutton(frame1a, text="Selection by point", command=pointfunc);
poly.but <- tkbutton(frame1a, text="Selection by polygon ", command=polyfunc);
tkpack(point.but, poly.but, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame1a, expand = "TRUE", fill = "x")

frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1c, text = "Work on the histogram", font = "Times 12",
foreground = "darkred", background = "white"))
barre.but <- tkbutton(frame1c, text="Select a bar  ", command=barfunc);
tkpack(barre.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1c, expand = "TRUE", fill = "x")


frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nettoy.but <- tkbutton(frame1b, text="     Reset selection     " , command=SGfunc);
tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1b, expand = "TRUE", fill = "x")


frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame2, text = "Options", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame2, text = "Spatial contours  ", font = "Times 11",
foreground = "darkred", background = "white"),tklabel(frame2, text = "Preselected sites  ", font = "Times 11",
foreground = "darkred", background = "white"),tklabel(frame2, text = "  Bubbles    ", font = "Times 11",
foreground = "darkred", background = "white"),side = "left", fill="x",expand = "TRUE")
tkpack(frame2, expand = "TRUE", fill = "x")

frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nocou1.but <- tkbutton(frame2b, text="On/Off", command=cartfunc)
noint1.but <- tkbutton(frame2b, text="On/Off", command=fnointer)
bubble.but <- tkbutton(frame2b, text="On/Off", command=fbubble)
tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame2b, expand = "TRUE", fill = "x")



frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
foreground = "blue", background = "white"))

quit.but2 <- tkbutton(frame3, text="Exit", command=quitfunc);

tkpack(quit.but2, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame3, expand = "TRUE", fill = "x")


}
####################################################

return(invisible())
}
