`neighbourmap` <-function(sp.obj, name.var, nb.obj, lin.reg=TRUE,
names.attr=names(sp.obj), criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16, col="lightblue3",
xlab="", ylab="", axes=FALSE, lablong="", lablat="")
{
envir = as.environment(1)
# Verification of the Spatial Object sp.obj
class.obj<-class(sp.obj)[1]
spdf<-(class.obj=="SpatialPolygonsDataFrame")
if(substr(class.obj,1,7)!="Spatial") stop("sp.obj may be a Spatial object")
if(substr(class.obj,nchar(class.obj)-8,nchar(class.obj))!="DataFrame") stop("sp.obj should contain a data.frame")
if(!is.numeric(name.var) & is.na(match(as.character(name.var),names(sp.obj)))) stop("name.var is not included in the data.frame of sp.obj")
if(length(names.attr)!=length(names(sp.obj))) stop("names.attr should be a vector of character with a length equal to the number of variable")

# we propose to refind the same arguments used in first version of GeoXp
long<-coordinates(sp.obj)[,1]
lat<-coordinates(sp.obj)[,2]

var<-sp.obj@data[,name.var]

# verify the type of the main variable
if(!(is.integer(var) || is.double(var))) stop("the variable name.var should be a numeric variable")


listvar<-sp.obj@data
listnomvar<-names.attr

# spatial weight matrix
W<-nb2mat(nb.obj)

# Code which was necessary in the previous version
# if(is.null(carte) & class.obj=="SpatialPolygonsDataFrame") carte<-spdf2list(sp.obj)$poly

 # for identifyng the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")



#initialisation
 obs <- matrix(FALSE, nrow=length(long), ncol=length(long))
 obs2 <- matrix(FALSE, nrow=length(long), ncol=length(long));
 nointer<-FALSE
 nocart<-FALSE
 buble<-FALSE
 legends<-list(FALSE,FALSE,"","")
 z<-NULL
 legmap<-NULL
 inout<-NULL
 graf<-"Neighbourplot1"
 labvar<-c(xlab,ylab)
 classe<-rep(1,length(long))

# Transformation data.frame en matrix
if((length(listvar)>0) && (dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()


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
    if(spdf & nrow(sp.obj)>75 & !buble) 
           {points(long,lat,pch=16,col='royalblue')}
   
    while(!quit)
    {
     dev.set(2)
     loc<-locator(1)
    
      if (is.null(loc)) 
       {
         quit<-TRUE
         carte(long=long, lat=lat, obs=obs, sp.obj=sp.obj, carte=carte,nocart=nocart, classe=classe,
         symbol=c(pch[1],16),W=W, method="Neighbourplot1", buble=buble, cbuble=z, criteria=criteria,
         nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
         label=label, cex.lab=cex.lab)    
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
        
      carte(long=long, lat=lat, obs=obs,    sp.obj=sp.obj, carte=carte,nocart=nocart, classe=classe,
      symbol=c(pch[1],16),W=W, method="Neighbourplot1", buble=buble, cbuble=z, criteria=criteria,
      nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
      label=label, cex.lab=cex.lab)      
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
      title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
      if(spdf & nrow(sp.obj)>75 & !buble) 
           {points(long,lat,pch=16,col='royalblue')}
           
      graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
      couleurs=col, symbol=pch, opt1=lin.reg, W=W)
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
      carte(long=long, lat=lat, obs=obs,    sp.obj=sp.obj, carte=carte,nocart=nocart,  classe=classe,
      symbol=c(pch[1],16),W=W, method="Neighbourplot1", buble=buble, cbuble=z, criteria=criteria,
      nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
      label=label, cex.lab=cex.lab)      
      
      graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
      couleurs=col, symbol=pch, opt1=lin.reg, W=W)
 #   obs <<- matrix(FALSE, nrow=length(long), ncol=length(long));

}
  }

####################################################
# sélection d'un point sur le scatterplot
####################################################

voisfunc <- function()
{
   if (graf=="Neighbourplot1")  SGfunc()
   
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
            graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
            couleurs=col, symbol=pch, opt1=lin.reg , W=W)
            next
        }
 
 
        obs <<- selectstat(var1=var,obs=obs, Xpoly=loc[1], Ypoly=loc[2], method="Neighbourplot", W=W)

        # graphiques
      carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=classe,
      symbol=c(pch[1],16),W=W, method="Neighbourplot2", buble=buble, cbuble=z, criteria=criteria,
      nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
      label=label, cex.lab=cex.lab)      
      
      graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
      couleurs=col, symbol=pch, opt1=lin.reg , W=W)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
      title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    }
  }



####################################################
# sélection d'un polygone sur le scattermap
####################################################


polyscatfunc <- function() 
{ 
  obs2 <<- matrix(FALSE, nrow=length(long), ncol=length(long))
  
  if (graf=="Neighbourplot1")  SGfunc()

  graf<<-"Neighbourplot2"
  quit <- FALSE
  polyX <- NULL
  polyY <- NULL
 
    dev.set(3)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

   while (!quit) 
    {
     dev.set(3)
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

  
  obs2[which(W!=0,arr.ind=TRUE)] <<- inout(cbind(var[which(W!=0,arr.ind=TRUE)[,1]],var[which(W!=0,arr.ind=TRUE)[,2]]),cbind(polyX, polyY), bound = TRUE)

  obs3<-obs+obs2
  obs[which(obs3==1,arr.ind=TRUE)]<<-TRUE
  obs[which(obs3!=1,arr.ind=TRUE)]<<-FALSE

  carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=classe,
  symbol=c(pch[1],16),W=W, method="Neighbourplot2", buble=buble, cbuble=z, criteria=criteria,
  nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
  label=label, cex.lab=cex.lab)      
      
  graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
  couleurs=col, symbol=pch, opt1=lin.reg , W=W)
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
   carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=classe,
   symbol=c(pch[1],16),W=W, method=graf, buble=buble, cbuble=z, criteria=criteria,
   nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
   label=label, cex.lab=cex.lab)       
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
    obs <<- matrix(FALSE, nrow=length(long), ncol=length(long))

  # graphiques
  carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=classe,
  symbol=c(pch[1],16),W=W, method=graf, buble=buble, cbuble=z, criteria=criteria,
  nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
  label=label, cex.lab=cex.lab)  
  
  graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
  couleurs=col, symbol=pch, opt1=lin.reg , W=W)
}

####################################################
# Open a no interactive selection
####################################################

fnointer<-function() 
{
 if (length(criteria) != 0)
 {
  ifelse(!nointer,nointer<<-TRUE,nointer<<-FALSE)
  carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart,  classe=classe,
  symbol=c(pch[1],16),W=W, method=graf, buble=buble, cbuble=z, criteria=criteria,
  nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
  label=label, cex.lab=cex.lab)  
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
  
  carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart,  classe=classe,
  symbol=c(pch[1],16),W=W, method=graf, buble=buble, cbuble=z, criteria=criteria,
  nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
  label=label, cex.lab=cex.lab)  
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
    if(graf=="Neighbourplot1")
    {obs<-unique(which(obs,arr.ind=TRUE)[,1])
     p<-length(obs)
     res<-NULL
     for(k in 1:p)
     {res<-rbind(res,cbind(obs[k],nb.obj[[obs[k]]]))
     }
    }
    else
    {res<-which(obs,arr.ind=TRUE)}
    assign("last.select", res, envir = envir)
}

####################################################
# Représentation des graphiques
####################################################

# Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2)  # new environment
 {
   carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=classe,
   symbol=c(pch[1],16),W=W, method="Neighbourplot1", buble=buble, cbuble=z, criteria=criteria,
   nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
   label=label, cex.lab=cex.lab)      
   
   graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
   couleurs=col, symbol=pch, opt1=lin.reg, W=W)
   assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {carte(long=long, lat=lat, obs=obs,   sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=classe,
   symbol=c(pch[1],16),W=W, method="Neighbourplot1", buble=buble, cbuble=z, criteria=criteria,
   nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
   label=label, cex.lab=cex.lab)      
   
   graphique(var1=var, obs=obs, num=3, graph="Neighbourplot", labvar=labvar,
   couleurs=col, symbol=pch, opt1=lin.reg, W=W)
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
tkwm.title(tt, "neighbourmap")

frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
foreground = "darkred", background = "white"))
point.but <- tkbutton(frame1a, text="Selection by point", command=pointfunc);
poly.but <- tkbutton(frame1a, text="Selection by polygon ", command=polyfunc);
tkpack(point.but, poly.but, side = "left", expand = "TRUE",fill = "x")
tkpack(frame1a, expand = "TRUE", fill = "x")

frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1c, text = "Work on the graph", font = "Times 12",
foreground = "darkred", background = "white"))
point.but2 <- tkbutton(frame1c, text="Selection by point", command=voisfunc);
poly.but2 <- tkbutton(frame1c, text="Selection by polygon ", command=polyscatfunc);
tkpack(point.but2, poly.but2, side = "left", expand = "TRUE",fill = "x")
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

quit.but <- tkbutton(frame3, text="Save results", command=quitfunc2);
quit.but2 <- tkbutton(frame3, text="Exit without saving", command=quitfunc);

tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame3, expand = "TRUE", fill = "x")
}

#############################################
return(invisible())
}

