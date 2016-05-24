`dblehistomap` <- function(sp.obj, names.var, nbcol=c(10,10), type = c("count","percent", "density"),
names.attr=names(sp.obj), criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16,
col=c("grey","lightblue3"), xlab=c("",""), ylab=c("count","count"), axes=FALSE, lablong="", lablat="")
{
envir = as.environment(1)
# Verification of the Spatial Object sp.obj
class.obj<-class(sp.obj)[1]
spdf<-(class.obj=="SpatialPolygonsDataFrame")
if(substr(class.obj,1,7)!="Spatial") stop("sp.obj may be a Spatial object")
if(substr(class.obj,nchar(class.obj)-8,nchar(class.obj))!="DataFrame") stop("sp.obj should contain a data.frame")
if(!is.numeric(names.var) & length(match(names.var,names(sp.obj)))!=length(names.var) ) stop("At least one component of names.var is not included in the data.frame of sp.obj")
if(length(names.attr)!=length(names(sp.obj))) stop("names.attr should be a vector of character with a length equal to the number of variable")



# we propose to refind the same arguments used in first version of GeoXp
long<-coordinates(sp.obj)[,1]
lat<-coordinates(sp.obj)[,2]

var1<-sp.obj@data[,names.var[1]]
var2<-sp.obj@data[,names.var[2]]

listvar<-sp.obj@data
listnomvar<-names.attr

# Code which was necessary in the previous version
# if(is.null(carte) & class.obj=="SpatialPolygonsDataFrame") carte<-spdf2list(sp.obj)$poly

 # for identifyng the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")

# initialisation
  nointer<-FALSE
  nocart<-FALSE
  buble<-FALSE
  legends<-list(FALSE,FALSE,"","")
  z<-NULL
  legmap<-NULL
  graphChoice <- ""
  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix<-""
  listgraph <- c("Histogram","Barplot","Scatterplot")
  labvar1<-c(xlab[1],ylab[1])
  labvar2<-c(xlab[2],ylab[2])
  obs <- vector(mode = "logical", length = length(long))

# options for adding a graphic with colors
  polyX2 <- NULL
  method <- ""
  col2 <- "blue"
  col3 <- "lightblue3"
  pch2<-pch[1]
  labmod<-""
  labvar1<-c(xlab[1],ylab[1])
  labvar2<-c(xlab[2],ylab[2])

# transformation data.frame en matrix
  if((length(listvar)>0) && (dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()
if(!(4%in%dev.list())) dev.new()

####################################################
# sélection d'un point
####################################################

pointfunc <- function()
 {
  quit <- FALSE

  dev.set(2)
  title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
  title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

  while (!quit)
   {
     dev.set(2)
     loc <- locator(1)
     if(spdf & nrow(sp.obj)>75 & !buble) 
       {points(long,lat,pch=16,col='royalblue')}
         
     if (is.null(loc))
      {
       quit <- TRUE
      carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
      symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
      lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])
          next
      }
    
     if(!spdf|nrow(sp.obj)>75)
       { 
       obs<<-selectmap(var1=long,var2=lat,obs=obs,Xpoly=loc[1], Ypoly=loc[2], method="point")}
       else
       {if(gContains(sp.obj,SpatialPoints(cbind(loc$x,loc$y),proj4string=CRS(proj4string(sp.obj)))))
        {for (i in 1:nrow(sp.obj))
          {if(gContains(sp.obj[i,],SpatialPoints(cbind(loc$x,loc$y),proj4string=CRS(proj4string(sp.obj)))))
           {obs[i]<<-!obs[i]
            break}  
          }
         } 
       }
   
    graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1],bin=type, labvar = labvar1, couleurs=col[1])
    graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2],bin=type, labvar = labvar2, couleurs=col[2])

    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
    symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
    lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

     if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
       {
        graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
        obs=obs, num=5, graph=graphChoice, couleurs=col3,symbol=pch, labvar=c(varChoice1,varChoice2));
       }
    }
 }


####################################################
# sélection d'un polygone
####################################################


    polyfunc <- function() {
        polyX <- NULL
        polyY <- NULL
        quit <- FALSE

        dev.set(2)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
         if(spdf & nrow(sp.obj)>75 & !buble) 
         {points(long,lat,pch=16,col='royalblue')}
         
        while (!quit) {
            dev.set(2)
            loc <- locator(1)
            if (is.null(loc))
            {
              quit <- TRUE
              next
            }
            polyX <- c(polyX, loc[1])
            polyY <- c(polyY, loc[2])
            lines(polyX, polyY)
       }
        polyX <- c(polyX, polyX[1])
        polyY <- c(polyY, polyY[1])

      if (length(polyX)>0)
      {
        lines(polyX, polyY)
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs,
            Xpoly = polyX, Ypoly = polyY, method = "poly")

  graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1],bin=type, labvar = labvar1, couleurs=col[1])
    graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2],bin=type, labvar = labvar2, couleurs=col[2])

    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
    symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
    lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

     if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
       {
        graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
        obs=obs, num=5, graph=graphChoice, couleurs=col3,symbol=pch, labvar=c(varChoice1,varChoice2));
       }
      }
  }


####################################################
# sélection d'une barre sur le premier histogramme
####################################################


    bar1func <- function()
    {
        SGfunc()
        quit <- FALSE

        dev.set(3)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

        while (!quit) {
            dev.set(3)
            loc <- locator(1)
            if (is.null(loc))
            {
              quit <- TRUE
              graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1],bin=type, labvar = labvar1, couleurs=col[1])
              next
            }
            obs <<- selectstat(var1 = var1, obs = obs, Xpoly = loc[1],Ypoly = loc[2], method = "Histogram", nbcol = nbcol[1])

     graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1],bin=type, labvar = labvar1, couleurs=col[1])
     title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
     title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

     graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2],bin=type, labvar = labvar2, couleurs=col[2])

     carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
     symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
     lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
       {
        graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
        obs=obs, num=5, graph=graphChoice, couleurs=col3,symbol=pch, labvar=c(varChoice1,varChoice2));
       }
      }
  }


####################################################
# sélection d'une barre sur le second histogramme
####################################################


    bar2func <- function()
    {
        SGfunc()
        quit <- FALSE

       dev.set(4)
       title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
       title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

        while (!quit)
         {
            dev.set(4)
            loc <- locator(1)
            if (is.null(loc))
            {
              quit <- TRUE
              graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2],bin=type, labvar = labvar2, couleurs=col[2])
              next
            }
            obs <<- selectstat(var1 = var2, obs = obs, Xpoly = loc[1],Ypoly = loc[2], method = "Histogram", nbcol = nbcol[2])

    graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1],bin=type, labvar = labvar1, couleurs=col[1])
    graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2],bin=type, labvar = labvar2, couleurs=col[2])
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
    symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
    lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

     if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
       {
        graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
        obs=obs, num=5, graph=graphChoice, couleurs=col3,symbol=pch, labvar=c(varChoice1,varChoice2));
       }
     }
  }




####################################################
# choix d'un autre graphique
####################################################

graphfunc <- function()
{
   if ((length(listvar) != 0) && (length(listnomvar) != 0))
    {
        choix <<- selectgraph(listnomvar,listgraph)
        varChoice1 <<- choix$varChoice1
        varChoice2 <<- choix$varChoice2
        graphChoice <<- choix$graphChoice

        if ((graphChoice != "") && (varChoice1 != ""))
        {
          if (((graphChoice == "Histogram")&&(!is.numeric(listvar[,which(listnomvar == varChoice1)])))||((graphChoice == "Scatterplot")&&((!is.numeric(listvar[,which(listnomvar == varChoice1)]))||(!is.numeric(listvar[,which(listnomvar == varChoice2)])))))
           {
            tkmessageBox(message="Variables choosed are not in a good format",icon="warning",type="ok")
           }
          else
           {
            res1<-choix.couleur(graphChoice,listvar,listnomvar,varChoice1,legends,col,pch,spdf=spdf)

            method <<- res1$method
            col2 <<- res1$col2
            col3 <<- res1$col3
            pch2 <<- res1$pch2
            legends <<- res1$legends
            labmod <<- res1$labmod

            if(!(5%in%dev.list())) dev.new()
            graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
            obs=obs, num=5, graph=graphChoice, couleurs=col3,symbol=pch, labvar=c(varChoice1,varChoice2))

            carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
            symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
            lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])
           }
       }
   }
   else
   {
    tkmessageBox(message="Variables (listvar) and their names (listnomvar) must have been given",icon="warning",type="ok")
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
    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
    symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
    lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])
   }
   else
   {
    tkmessageBox(message="Spatial contours have not been given",icon="warning",type="ok")
   }
}




####################################################
# rafraichissement des graphiques
####################################################
SGfunc<-function()
{
    obs<<-vector(mode = "logical", length = length(long));

 # graphiques

     graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1],bin=type, labvar = labvar1, couleurs=col[1])
     graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2],bin=type, labvar = labvar2, couleurs=col[2])

     carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
     symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
     lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
       {
        graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
        obs=obs, num=5, graph=graphChoice, couleurs=col3,symbol=pch, labvar=c(varChoice1,varChoice2))
       }

}



####################################################
# quitter l'application
####################################################

quitfunc<-function()
{
    #tclvalue(fin)<<-TRUE
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = baseenv())
   # assign("obs", row.names(sp.obj)[obs], envir = .GlobalEnv)
}

quitfunc2<-function()
{
    #tclvalue(fin)<<-TRUE
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = baseenv())
    print("Results have been saved in last.select object")
    assign("last.select", which(obs), envir = envir)
}


####################################################
# Open a no interactive selection
####################################################

fnointer<-function()
{
 if (length(criteria) != 0)
 {
  ifelse(!nointer,nointer<<-TRUE,nointer<<-FALSE)
  carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
  symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
  lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])
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

  carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
  symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
  lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

}


####################################################
# graphiques
####################################################

# Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2)
 {
  carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
  symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
  lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

  graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1], bin=type, labvar = labvar1, couleurs=col[1])
  graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2], bin=type, labvar = labvar2, couleurs=col[2])
  assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,  label=label,
  symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes, labmod=labmod,
  lablong=lablong,lablat=lablat,cex.lab=cex.lab,method=method,classe=listvar[,which(listnomvar == varChoice1)])

  graphique(var1 = var1, obs = obs, num = 3, graph = "Histogram",nbcol = nbcol[1], bin=type, labvar = labvar1, couleurs=col[1])
  graphique(var1 = var2, obs = obs, num = 4, graph = "Histogram",nbcol = nbcol[2], bin=type, labvar = labvar2, couleurs=col[2])
  assign("GeoXp.open", TRUE, envir = baseenv())}
 }
}

####################################################
# création de la boite de dialogue
####################################################

if(interactive())
{
fontheading<-tkfont.create(family="times",size=14,weight="bold")

tt <- tktoplevel()
tkwm.title(tt, "dblehistomap")

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
tkpack(tklabel(frame1c, text = "Select bar(s)", font = "Times 12",
foreground = "darkred", background = "white"))
barre1.but <- tkbutton(frame1c, text="on the 1st histogram", command=bar1func);
barre2.but <- tkbutton(frame1c, text="on the 2nd histogram", command=bar2func);
tkpack(barre1.but, barre2.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1c, expand = "TRUE", fill = "x")


frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nettoy.but <- tkbutton(frame1b, text="     Reset selection     " , command=SGfunc);
tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1b, expand = "TRUE", fill = "x")


frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame2, text = "Options", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame2, text = "Spatial contours", font = "Times 11",
foreground = "darkred", background = "white"),tklabel(frame2, text = "Preselected sites", font = "Times 11",
foreground = "darkred", background = "white"), side = "left", fill="x",expand = "TRUE")
tkpack(frame2, expand = "TRUE", fill = "x")

frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nocou1.but <- tkbutton(frame2b, text="On/Off", command=cartfunc)
noint1.but <- tkbutton(frame2b, text="On/Off", command=fnointer)
tkpack(nocou1.but,noint1.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame2b, expand = "TRUE", fill = "x")

frame2c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame2c, text = "Bubbles", font = "Times 11",
foreground = "darkred", background = "white"),tklabel(frame2c, text = "Additional graph", font = "Times 11",
foreground = "darkred", background = "white"), side = "left", fill="x",expand = "TRUE")
tkpack(frame2c, expand = "TRUE", fill = "x")

frame2d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
bubble.but <- tkbutton(frame2d, text="On/Off", command=fbubble)
autre.but <- tkbutton(frame2d, text="     OK     " , command=graphfunc)
tkpack(bubble.but,autre.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame2d, expand = "TRUE", fill = "x")

frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
foreground = "blue", background = "white"))

quit.but <- tkbutton(frame3, text="Save results", command=quitfunc2);
quit.but2 <- tkbutton(frame3, text="Exit without saving", command=quitfunc);

tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame3, expand = "TRUE", fill = "x")

}
#########################################
return(invisible())
}

