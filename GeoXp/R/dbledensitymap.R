`dbledensitymap` <- function(sp.obj, names.var, kernel='triweight',
names.attr=names(sp.obj), criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16,
col=c("grey","lightblue3"), xlab=c("",""), ylab="", axes=FALSE, lablong="", lablat="")
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
  obs<-vector(mode = "logical", length = length(long))
  graph1<-"Densityplot2"
  graph2<-"Densityplot2"
  nointer<-FALSE
  nocart<-FALSE
  z<-NULL
  legmap<-NULL
  interv<-NULL
  buble<-FALSE
  legends<-list(FALSE,FALSE,"","")

  # for the slider
  alpha11<-20
  alpha21<-20
  names.slide<-c("Alpha (1st graph)","Alpha (2nd graph)")
  
  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix<-""

  listgraph <- c("Histogram","Barplot","Scatterplot")

  # options for adding a graphic with colors
  polyX2 <- NULL
  method <- ""
  col2 <- "blue"
  col3<-"lightblue3"
  pch2<-pch[1]
  labmod<-""
  labvar1<-c(xlab[1],ylab[1])
  labvar2<-c(xlab[2],ylab[2])


# Transformation de data.frame en matrix
if((length(listvar)>0)&&(dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()
if(!(4%in%dev.list())) dev.new()

####################################################
# sélection d'un point
####################################################

pointfunc<-function()
{
   if((graph1=="Densityplot2")||(graph2=="Densityplot2"))
    {
     #SGfunc()
      graph1<<-"Densityplot1"
      graph2<<-"Densityplot1"
      
       if (length(var1[obs]) > 1)
        {
         graphique(var1=var1, obs=obs, alpha1=alpha11,  num=3, graph=graph1, labvar=labvar1,
         couleurs=col[1],symbol=pch,kernel=kernel)

         graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2, labvar=labvar2,
         couleurs=col[2],symbol=pch,kernel=kernel)
        }
       else
       {dev.set(3)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
       dev.set(4)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }
    }

    quit <- FALSE
    loc <- NULL

    dev.set(2)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
      if(spdf & nrow(sp.obj)>75 & !buble) 
       {points(long,lat,pch=16,col='royalblue')}
    
    while(!quit)
    {
        dev.set(2)
        loc<-locator(1)

    graph1<<-"Densityplot1"
    graph2<<-"Densityplot1"

        if(is.null(loc))
        {
          quit<-TRUE
        carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
        label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
        lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
        classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
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

        carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
        label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
        lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
        classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
        if(spdf & nrow(sp.obj)>75 & !buble) 
        {points(long,lat,pch=16,col='royalblue')}
        # graphiques
      if (length(var1[obs]) > 1)
        {
         graphique(var1=var1, obs=obs, alpha1=alpha11,  num=3, graph=graph1, labvar=labvar1,
         couleurs=col[1],symbol=pch,kernel=kernel)

         graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2, labvar=labvar2,
         couleurs=col[2],symbol=pch,kernel=kernel)
        }
       else
       {dev.set(3)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
       dev.set(4)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }


         if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
          {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
            obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
          }
    }
  }

####################################################
# sélection d'un polygone
####################################################

polyfunc<-function()
{

   if((graph1=="Densityplot2")||(graph2=="Densityplot2"))
    {
     #SGfunc()
      graph1<<-"Densityplot1"
      graph2<<-"Densityplot1"

      graphique(var1=var1, obs=obs, alpha1=alpha11,  num=3, graph=graph1, labvar=labvar1,
      couleurs=col[1],symbol=pch,kernel=kernel)

      graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2, labvar=labvar2,
      couleurs=col[2],symbol=pch,kernel=kernel)

    }

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


      if (length(var1[obs]) > 1)
        {
         graphique(var1=var1, obs=obs, alpha1=alpha11,  num=3, graph=graph1, labvar=labvar1,
         couleurs=col[1],symbol=pch,kernel=kernel)

         graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2, labvar=labvar2,
         couleurs=col[2],symbol=pch,kernel=kernel)
        }
       else
       {dev.set(3)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
       dev.set(4)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }


    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
    label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
    lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
    classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)


    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
     {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
      obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
     }
}

 }

####################################################
# sélection d'un intervalle sous la courbe de densité
####################################################

inter1func<-function()
{
  # SGfunc();

#    quit <- FALSE;

  if(graph1=="Densityplot1"||(graph1==graph2))
   {
    SGfunc()
    graph1<<-"Densityplot2"
    graph2<<-"Densityplot1"
   }

    polyX <- NULL
    n.inter<-length(polyX2)

    while (length(polyX)<2)
    {
      dev.set(3)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
      title(sub = "Click two times to select an interval", cex.sub = 0.8, font.sub = 3,col.sub='red')
      loc<-locator(1)
      polyX <- c(polyX, loc[1])
    }

    polyX2[[n.inter+1]]<<-polyX

    obs<<-selectstat(var1=var1,obs=obs,Xpoly=polyX[1], Ypoly=polyX[2],method="Densityplot")

    # graphiques

         graphique(var1=var1, obs=obs, alpha1=alpha11,   num=3, graph=graph1, Xpoly=polyX2,
         labvar=labvar1, couleurs=col[1], kernel=kernel)

       if (length(var1[obs]) > 1)
        {
        graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2,
        labvar=labvar2, couleurs=col[2], kernel=kernel)
         
        }
      else
       {dev.set(4)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
       }
     
        carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
        label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
        lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
        classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)


    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
       obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
      }

  }

####################################################
# sélection d'un intervalle sous la courbe de densité
####################################################

inter2func<-function()
{
  if(graph2=="Densityplot1"||(graph1==graph2))
   {
    SGfunc()
    graph1<<-"Densityplot1"
    graph2<<-"Densityplot2"
   }

    polyX <- NULL
    n.inter<-length(polyX2)

    while (length(polyX)<2)
    {
      dev.set(4)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
      title(sub = "Click two times to select an interval", cex.sub = 0.8, font.sub = 3,col.sub='red')
      loc<-locator(1)
      polyX <- c(polyX, loc[1])
    }

    polyX2[[n.inter+1]]<<-polyX

    obs<<-selectstat(var1=var2,obs=obs,Xpoly=polyX[1], Ypoly=polyX[2],method="Densityplot")

    # graphiques
     graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2, Xpoly=polyX2,
     labvar=labvar2, couleurs=col[2], kernel=kernel)

     if (length(var2[obs]) > 1)
     {
     graphique(var1=var1, obs=obs, alpha1=alpha11,  num=3, graph=graph1,
     labvar=labvar1, couleurs=col[1],kernel=kernel)
     }
     else
     {dev.set(3)
      title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }
     
     carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
     label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
     lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
     classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)

    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
       obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
      }

 }


####################################################
# Choisir une valeur sur le 1er graphique
####################################################

choixvalue1 <- function()
{
  if(graph1=="Densityplot1"||(graph1==graph2))
   {
    SGfunc()
    graph1<<-"Densityplot2"
    graph2<<-"Densityplot1"
   }

   dev.set(3)
   title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
      
      
  tt1<-tktoplevel()
  Name <- tclVar("1st value")
  Name2 <- tclVar("2nd value")
  entry.Name <-tkentry(tt1,width="8",textvariable=Name)
  entry.Name2 <-tkentry(tt1,width="8",textvariable=Name2)
  tkgrid(tklabel(tt1,text="Please enter values"),entry.Name,entry.Name2)

  OnOK <- function()
  {
	 value1 <- tclvalue(Name)
	 value2 <- tclvalue(Name2)
	 n.inter<-length(polyX2)
	 tkdestroy(tt1)

   if(is.na(as.numeric(value1))||is.na(as.numeric(value2)))
    {
        tkmessageBox(message="Sorry, but you have to choose decimal values",icon="warning",type="ok");
    }
      else
    {

    polyX2[[n.inter+1]]<<- c(as.numeric(value1),as.numeric(value2))

    obs<<-selectstat(var1=var1,obs=obs,Xpoly=as.numeric(value1), Ypoly=as.numeric(value2),method="Densityplot")

    # graphiques
    graphique(var1=var1, obs=obs, alpha1=alpha11,   num=3, graph=graph1, Xpoly=polyX2,
    labvar=labvar1, couleurs=col[1], kernel=kernel)

    if (length(var1[obs]) > 1)
     {
    graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2,
    labvar=labvar2, couleurs=col[2], kernel=kernel)
    }
    else
    {dev.set(4)
     title(sub = "You have to choose one more site to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
    }

    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
    label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
    lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
    classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)

    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
       obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
      }
    }
  }

  OK.but <-tkbutton(tt1,text="   OK   ",command=OnOK)
  #tkbind(entry.Name, "<Return>",OnOK)
  tkgrid(OK.but)
  tkfocus(tt1)

}

####################################################
# Choisir une valeur sur le deuxième graphique
####################################################

choixvalue2 <- function()
{
  if(graph2=="Densityplot1"||(graph1==graph2))
   {
    SGfunc()
    graph1<<-"Densityplot1"
    graph2<<-"Densityplot2"
   }

  dev.set(4)
  title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
      
  tt1<-tktoplevel()
  Name <- tclVar("1st value")
  Name2 <- tclVar("2nd value")
  entry.Name <-tkentry(tt1,width="8",textvariable=Name)
  entry.Name2 <-tkentry(tt1,width="8",textvariable=Name2)
  tkgrid(tklabel(tt1,text="Please enter values"),entry.Name,entry.Name2)


  OnOK <- function()
  {
	 value1 <- tclvalue(Name)
	 value2 <- tclvalue(Name2)
	 n.inter<-length(polyX2)
	 tkdestroy(tt1)

   if(is.na(as.numeric(value1))||is.na(as.numeric(value2)))
    {
        tkmessageBox(message="Sorry, but you have to choose decimal values",icon="warning",type="ok")
    }
   else
    {

    polyX2[[n.inter+1]]<<- c(as.numeric(value1),as.numeric(value2))

    obs<<-selectstat(var1=var2,obs=obs,Xpoly=as.numeric(value1), Ypoly=as.numeric(value2),method="Densityplot")

    # graphiques
     graphique(var1=var2, obs=obs, alpha1=alpha21,  num=4, graph=graph2, Xpoly=polyX2,
     labvar=labvar2, couleurs=col[2], kernel=kernel)

    if (length(var1[obs]) > 1)
     {
     graphique(var1=var1, obs=obs, alpha1=alpha11,  num=3, graph=graph1,
     labvar=labvar1, couleurs=col[1],kernel=kernel)
     }
     else
     {dev.set(3)
      title(sub = "You have to choose one more site to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }

     carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
     label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
     lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
     classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)

    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
       obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
      }
    }
  }

OK.but <-tkbutton(tt1,text="   OK   ",command=OnOK)
#tkbind(entry.Name, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt1)

}
####################################################
# modification du alpha pour la courbe de densité
####################################################


refresh1.code<-function(...)
{
 res<-slider1(names.slide=names.slide,no=1)
 alpha11<<-res$alpha11
 alpha21<<-res$alpha21
 
 if(graph1=="Densityplot1")
 {   if (length(var2[obs]) > 1)
       {graphique(var1=var1, obs=obs, alpha1=alpha11, num=3, graph=graph1, labvar=labvar1, couleurs=col[1],kernel=kernel)
       }
       else
       {dev.set(3)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
      }
  }
 else
 {graphique(var1=var1, obs=obs,alpha1=alpha11,  num=3, graph="Densityplot2", Xpoly=polyX2,
  labvar=labvar1, couleurs=col[1],kernel=kernel)
 }
 
 if(graph2=="Densityplot1")
  {   if (length(var2[obs]) > 1)
      {graphique(var1=var2, obs=obs, alpha1=alpha21, num=4, graph=graph2, labvar=labvar2,
       couleurs=col[2],kernel=kernel)
      }
      else
      {dev.set(4)
       title(sub = "You have to choose at least two sites to represent the sub-density", cex.sub = 0.8, font.sub = 3,col.sub='red')
      }
  }
 else
  {graphique(var1=var2, obs=obs,alpha1=alpha21, num=4, graph="Densityplot2", Xpoly=polyX2,
   labvar=labvar2, couleurs=col[2],kernel=kernel)}
 }
 


####################################################
# rafraichissement des graphiques
####################################################

SGfunc<-function()
{
    obs<<-vector(mode = "logical", length = length(long));
    polyX2 <<- NULL

    graphique(var1=var2, obs=obs, alpha1=alpha21,   num=4, graph="Densityplot1",
    labvar=labvar2, couleurs=col[2], kernel=kernel)

    graphique(var1=var1, obs=obs, alpha1=alpha11,   num=3, graph="Densityplot1", labvar=labvar1,
    couleurs=col[1], kernel=kernel)

    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
    label=label,cex.lab=cex.lab, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
    lablong=lablong, lablat=lablat,symbol=pch2, couleurs=col2,method=method,
    classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)

    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      {graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
       obs=obs, num=5, graph=graphChoice, symbol=pch, labvar=c(varChoice1,varChoice2),couleurs=col3)
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
# contour des unités spatiales
####################################################
cartfunc <- function()
{
 if (length(carte) != 0)
   {
    ifelse(!nocart,nocart<<-TRUE,nocart<<-FALSE)
    carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
    label=label,cex.lab=cex.lab, symbol=pch2, couleurs=col2, carte=carte,nocart=nocart,legmap=legmap,legends=legends,
    axis=axes,lablong=lablong, lablat=lablat,method=method,classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
    }
   else
   {
    tkmessageBox(message="Spatial contours have not been given",icon="warning",type="ok")
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
            tkmessageBox(message="Variables choosed are not in a good format",icon="warning",type="ok");
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
            obs=obs, num=5, graph=graphChoice, couleurs=col3, symbol=pch, labvar=c(varChoice1,varChoice2));

            carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
            label=label,cex.lab=cex.lab, symbol=pch2, couleurs=col2,carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
            lablong=lablong, lablat=lablat,method=method,classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
            }
       }
   }
   else
   {
    tkmessageBox(message="Variables (listvar) and their names (listnomvar) must have been given",icon="warning",type="ok");
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
   carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
   label=label,cex.lab=cex.lab, symbol=pch2, couleurs=col2, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
   lablong=lablong, lablat=lablat,method=method,classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
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

  carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj, buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
  label=label,cex.lab=cex.lab, symbol=pch2, couleurs=col2, carte=carte,nocart=nocart,legmap=legmap,legends=legends,axis=axes,
  lablong=lablong, lablat=lablat,method=method,classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)

}

####################################################
# Représentation graphique
####################################################

# Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2)
 {
  graphique(var1=var2, obs=obs, alpha1=alpha21, num=4, graph=graph1, labvar=labvar2,
  couleurs=col[2],kernel=kernel,Xpoly=NULL)

  graphique(var1=var1, obs=obs, alpha1=alpha11, num=3, graph=graph2, labvar=labvar1,
  couleurs=col[1],kernel=kernel,Xpoly=NULL)

  carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
  label=label,cex.lab=cex.lab, symbol=pch2, couleurs=col2, carte=carte,nocart=nocart,legmap=legmap,legends=legends,
  axis=axes, lablong=lablong, lablat=lablat,method=method,classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
  assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {graphique(var1=var2, obs=obs, alpha1=alpha21, num=4, graph=graph1, labvar=labvar2,
   couleurs=col[2],kernel=kernel,Xpoly=NULL)

   graphique(var1=var1, obs=obs, alpha1=alpha11, num=3, graph=graph2, labvar=labvar1,
   couleurs=col[1],kernel=kernel,Xpoly=NULL)

   carte(long=long, lat=lat,obs=obs,  sp.obj=sp.obj,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,
   label=label,cex.lab=cex.lab, symbol=pch2, couleurs=col2, carte=carte,nocart=nocart,legmap=legmap,legends=legends,
   axis=axes, lablong=lablong, lablat=lablat,method=method,classe=listvar[,which(listnomvar == varChoice1)],labmod=labmod)
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
tkwm.title(tt, "dbledensitymap")

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
tkpack(tklabel(frame1c, text = "Select an interval on the 1st graphic", font = "Times 12",
foreground = "darkred", background = "white"))
intervalle1.but <- tkbutton(frame1c, text="by selecting on graph", command=inter1func);
intervalle11.but <- tkbutton(frame1c, text="by specifying bounds", command=choixvalue1);
tkpack(intervalle1.but,intervalle11.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1c, expand = "TRUE", fill = "x")

frame1d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1d, text = "Select an interval on the 2nd graphic", font = "Times 12",
foreground = "darkred", background = "white"))
intervalle2.but <- tkbutton(frame1d, text="by selecting on graph", command=inter2func);
intervalle22.but <- tkbutton(frame1d, text="by specifying bounds", command=choixvalue2);
tkpack(intervalle2.but,intervalle22.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1d, expand = "TRUE", fill = "x")

frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nettoy.but <- tkbutton(frame1b, text="     Reset selection     " , command=SGfunc);
tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1b, expand = "TRUE", fill = "x")


frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame2, text = "Options", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame2, text = "Spatial contours", font = "Times 10",
foreground = "darkred", background = "white"),tklabel(frame2, text = "Preselected sites", font = "Times 10",
foreground = "darkred", background = "white"),tklabel(frame2, text = "Bubbles", font = "Times 10",
foreground = "darkred", background = "white"),tklabel(frame2, text = "Additional graph", font = "Times 10",
foreground = "darkred", background = "white"), side = "left", fill="x",expand = "TRUE")
tkpack(frame2, expand = "TRUE", fill = "x")

frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nocou1.but <- tkbutton(frame2b, text="On/Off", command=cartfunc)
noint1.but <- tkbutton(frame2b, text="On/Off", command=fnointer)
bubble.but <- tkbutton(frame2b, text="On/Off", command=fbubble)
autre.but <- tkbutton(frame2b, text="     OK     " , command=graphfunc)
tkpack(nocou1.but,noint1.but,bubble.but,autre.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame2b, expand = "TRUE", fill = "x")


frame2e <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
slider1(frame2e,refresh1.code,names.slide,3,100,1,c(alpha11,alpha21))
tkpack(frame2e, expand = "TRUE", fill = "x")

frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
foreground = "blue", background = "white"))

quit.but <- tkbutton(frame3, text="Save results", command=quitfunc2);
quit.but2 <- tkbutton(frame3, text="Exit without saving", command=quitfunc)

tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame3, expand = "TRUE", fill = "x")
}
####################################################

return(invisible())
}

