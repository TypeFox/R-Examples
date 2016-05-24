`clustermap` <- function(sp.obj, names.var, clustnum, method=c("kmeans","hclust"),type=NULL, centers=NULL, scale=FALSE,
names.arg="", names.attr=names(sp.obj), criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16,
col="lightblue3", xlab="Cluster", ylab="Number", axes=FALSE, lablong="", lablat="")
{
envir = as.environment(1)
# Verification of the Spatial Object sp.obj
class.obj<-class(sp.obj)[1]

if(substr(class.obj,1,7)!="Spatial") stop("sp.obj may be a Spatial object")
if(substr(class.obj,nchar(class.obj)-8,nchar(class.obj))!="DataFrame") stop("sp.obj should contain a data.frame")
if(!is.numeric(names.var) & length(match(names.var,names(sp.obj)))!=length(names.var) ) stop("At least one component of names.var is not included in the data.frame of sp.obj")
if(length(names.attr)!=length(names(sp.obj))) stop("names.attr should be a vector of character with a length equal to the number of variable")

# Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2)
 {
  assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {assign("GeoXp.open", TRUE, envir = baseenv())}
 }
}

# we propose to refind the same arguments used in first version of GeoXp
long<-coordinates(sp.obj)[,1]
lat<-coordinates(sp.obj)[,2]

dataset <- sp.obj@data[,names.var]

listvar<-sp.obj@data
listnomvar<-names.attr

# Code which was necessary in the previous version
 if(is.null(carte) & class.obj=="SpatialPolygonsDataFrame") carte<-spdf2list(sp.obj)$poly

 # for colors in map and new grahics
 ifelse(length(col)==1, col2<-"blue", col2<-col)
 col3<-"lightblue3"
 
 # for identifyng the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")


  # initialisation
  nointer<-FALSE
  nocart<-FALSE
  buble<-FALSE
  z<-NULL
  legmap<-NULL
  legends<-list(FALSE,FALSE,"","")
  labvar=c(xlab,ylab)
  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix <- ""
  listgraph <- c("Histogram","Barplot","Scatterplot")
  if(names.arg[1]=="") names.arg <- 1:clustnum
  obs<-vector(mode = "logical", length = length(long))

  # Méthodes de classification
  if(length(type)==0)
  {ifelse(method[1]=="kmeans",type<-"Hartigan-Wong", type <-"ward")}

  # Réduction de la matrice des données
  if(scale && class(dataset)!="dist") dataset<-scale(dataset)

  # Etude des différentes possibilités
  if(class(dataset)!="dist")
  {if(method[1]=="hclust") dataset<-dist(dataset)}
  else
  {if(method[1]=="kmeans") dataset<-as.matrix(dataset)}

  # classification
  if(method[1]=="kmeans")
  {
    num.graph <- 4 
    ifelse(length(centers)==0,res <- kmeans(dataset,clustnum,algorithm=type),
    res <- kmeans(dataset,centers,algorithm=type))
    vectclass <- res$cluster
  }
  else
  {
    num.graph <- 5
    res <- hclust(dataset,method=type)
    vectclass <- as.vector(cutree(res,k=clustnum))
  }


  # Transformation de data.frame en matrix
  if((length(listvar)>0) && (dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()

  if(method[1]=="hclust")
  {if(!(4%in%dev.list())) dev.new()}


####################################################
# sélection d'un point
####################################################

pointfunc<-function() 
{
    quit <- FALSE
   
     dev.set(2)
     title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
     title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    while(!quit)
    {
        #sélection des points
        dev.set(2)
        loc<-locator(1)
        if(is.null(loc)) 
        {
            quit<-TRUE
        carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
        lablong=lablong, lablat=lablat, label=label, symbol=pch,carte=carte,nocart=nocart,
        method="Cluster",classe=vectclass,couleurs=col2,legmap=legmap,legends=legends,labmod=names.arg,
        cex.lab=cex.lab,axis=axes)
            next
        }   
        obs<<-selectmap(var1=long,var2=lat,obs=obs,Xpoly=loc[1], Ypoly=loc[2], method="point")         

        # graphiques

        carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
        lablong=lablong, lablat=lablat, label=label, symbol=pch,carte=carte,nocart=nocart,
        method="Cluster",classe=vectclass,couleurs=col2,legmap=legmap,legends=legends,labmod=names.arg,
        cex.lab=cex.lab,axis=axes)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

        graphique(var1=vectclass, obs=obs, num=3,graph="Barplot", labvar=labvar, symbol=pch,labmod=names.arg,couleurs=col);
        
        if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        {
         graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
         obs=obs, num=num.graph, graph=graphChoice, couleurs=col3, labvar=c(varChoice1,varChoice2))
        }


    }
  }

####################################################
# sélection d'un polygone
####################################################

polyfunc<-function() 
{
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
 
    dev.set(2)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

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

    #graphiques
   carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
   lablong=lablong, lablat=lablat, label=label, symbol=pch,couleurs=col2,carte=carte,nocart=nocart,
   method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes,cex.lab=cex.lab)
   
   graphique(var1=vectclass, obs=obs, num=3,graph="Barplot", labvar=labvar, symbol=pch,labmod=names.arg,
   couleurs=col)
   
        if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        {
         graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
         obs=obs, num=num.graph, graph=graphChoice, couleurs=col3, labvar=c(varChoice1,varChoice2))
        }

}

  }

####################################################
# sélection d'une barre sur le diagramme
####################################################

barfunc<-function()
{
    SGfunc()
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
        graphique(var1=vectclass, obs=obs, num=3,graph="Barplot", labvar=labvar, symbol=pch,labmod=names.arg,
        couleurs=col)
          next
        }           
        obs<<-selectstat(var1=as.vector(vectclass),obs=obs,Xpoly=loc[1], Ypoly=loc[2],method="Barplot");   

        # graphiques
        graphique(var1=vectclass, obs=obs, num=3,graph="Barplot", labvar=labvar, symbol=pch,labmod=names.arg,
        couleurs=col)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

        carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
        lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,couleurs=col2,
        method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)

        if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        {
         graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
         obs=obs, num=num.graph, graph=graphChoice, couleurs=col3, labvar=c(varChoice1,varChoice2))
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
         if(method[1]=="kmeans")
         {if(!(4%in%dev.list())) dev.new()}
         else
         {if(!(5%in%dev.list())) dev.new()}
         graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
         obs=obs, num=num.graph, graph=graphChoice, couleurs=col3, labvar=c(varChoice1,varChoice2))
        }
        else
        {
        tkmessageBox(message="You must choose a variable",icon="warning",type="ok")
        }
        
    }
    else
    {
        tkmessageBox(message="You must give listvar and listnomvar in input",icon="warning",type="ok")
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
   carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
   lablong=lablong, lablat=lablat, label=label, cex.lab=cex.lab, symbol=pch,couleurs=col2,carte=carte,
   nocart=nocart, method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)    
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
   
    carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
    lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,
    couleurs=col2,method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)
    
    graphique(var1=vectclass, obs=obs, num=3, graph="Barplot", labvar=labvar, symbol=pch,labmod=names.arg,couleurs=col)
    
     if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        {
         graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
         obs=obs, num=num.graph, graph=graphChoice, couleurs=col3, labvar=c(varChoice1,varChoice2))
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
    print("Results have saved in last.select object")
    res<-list(obs=which(obs),vectclass=vectclass)
    assign("last.select", res, envir = envir)
    
}
####################################################
# Open a no interactive selection
####################################################

fnointer<-function() 
{
 if (length(criteria) != 0)
 {
   ifelse(!nointer,nointer<<-TRUE,nointer<<-FALSE)
   carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
   lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,
   couleurs=col2,method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)
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
  
  carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
  lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,
  couleurs=col2,method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)
 
}

####################################################
# Représentation graphique
####################################################

carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,
couleurs=col2,method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)
    
graphique(var1=as.vector(vectclass), obs=obs, num=3, graph="Barplot", labvar=labvar, symbol=pch,
labmod=names.arg,couleurs=col)

if(method[1]=="hclust")
{dev.set(4)
plot(res)
}


####################################################
# création de la boite de dialogue to create legens
####################################################

if(interactive())
{
 OnOK <- function()
 { 
  tkdestroy(tt1)	
  msg <- paste("Click on the map to indicate the location of the upper left corner of the legend box")
	tkmessageBox(message=msg)

  dev.set(2)
  loc <- locator(1)
  loc$name <- "Cluster"
  legends<<-list(legends[[1]],TRUE,legends[[3]],loc)

  carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
  lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,
  couleurs=col2,method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)

  graphique(var1=as.vector(vectclass), obs=obs, num=3, graph="Barplot", labvar=labvar, symbol=pch,
  labmod=names.arg,couleurs=col)
 }

OnOK2 <- function()
 { 
  legends<<-list(legends[[1]],FALSE,legends[[3]],"")
  tkdestroy(tt1)	

  carte(long=long, lat=lat,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,obs=obs,
  lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, symbol=pch,carte=carte,nocart=nocart,
  couleurs=col2,method="Cluster",classe=vectclass,legmap=legmap,legends=legends,labmod=names.arg,axis=axes)
    
  graphique(var1=as.vector(vectclass), obs=obs, num=3, graph="Barplot", labvar=labvar, 
  symbol=pch,labmod=names.arg,couleurs=col)
  }


if(length(col)==length(levels(as.factor(vectclass)))||length(pch)==length(levels(as.factor(vectclass))))
 {
  tt1<-tktoplevel()

  labelText12 <- tclVar("Do you want a legend for factors")
  label12 <- tklabel(tt1,justify = "center", wraplength = "3i", text=tclvalue(labelText12))
  tkconfigure(label12, textvariable=labelText12)
  tkgrid(label12,columnspan=2)

  point.but <- tkbutton(tt1, text="  Yes  ", command=OnOK)
  poly.but <- tkbutton(tt1, text=" No ", command=OnOK2)
  tkgrid(point.but, poly.but)
  tkgrid(tklabel(tt1,text="    "))

  tkfocus(tt1)
 }
}
####################################################
# création de la boite de dialogue
####################################################
if(interactive())
{
fontheading<-tkfont.create(family="times",size=14,weight="bold")

tt <- tktoplevel()
tkwm.title(tt, "clustermap")

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
tkpack(tklabel(frame1c, text = "Work on the bar plot", font = "Times 12",
foreground = "darkred", background = "white"))
barre.but <- tkbutton(frame1c, text="  Bar  ", command=barfunc);
tkpack(barre.but, side = "left", expand = "TRUE", fill = "x")
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
quit.but2 <- tkbutton(frame3, text="Exit without saving", command=quitfunc)

tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame3, expand = "TRUE", fill = "x")

}


####################################################
# Fin
####################################################

return(invisible())
  }

