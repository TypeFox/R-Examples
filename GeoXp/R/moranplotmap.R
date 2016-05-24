`moranplotmap` <-
function(sp.obj, name.var, listw.obj, flower=FALSE, locmoran=FALSE, names.arg=c("H.-H.","L.-H.","L.-L.","H.-L."),
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

# Code which was necessary in the previous version
# if(is.null(carte) & class.obj=="SpatialPolygonsDataFrame") carte<-spdf2list(sp.obj)$poly

 # for identifyng the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")

# We create a spatial weight matrix by using a matrix object
  n <- nrow(sp.obj)
  W<-matrix(0,n,n)
  W.sn<-listw2sn(listw.obj)
  W[as.matrix(W.sn[,1:2])]<-W.sn[,3]

  # Is W normalized ?
  is.norm<-all(apply(W,1,sum)==rep(1,n))
  
  #initialisation
if(xlab=="") xlab=name.var
if(ylab=="") ylab=paste("spatially lagged",name.var)
  obs<-vector(mode = "logical", length = length(long))
  obsq<-rep(0,length(long))
  nointer<-FALSE
  nocart<-FALSE
  buble<-FALSE
  buble2<-FALSE
  maptest<-FALSE
  classe<-rep(1,length(long))

  legends<-list(FALSE,FALSE,"","")
  legends2<-list(FALSE,FALSE,"","")

  z<-NULL
  z2<-NULL
  legmap<-NULL
  legmap2<-NULL
  num<-NULL
  labvar=c(xlab,ylab)

  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix <- ""
  listgraph <- c("Histogram","Barplot","Scatterplot")

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()


quad <- FALSE

# Transformation data.frame en matrix
if((length(listvar)>0)&&(dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)

# Option sur le moran
ifelse(flower,method <- "Neighbourplot1", method <- "")

choix.col<-FALSE

graph <- "Moran"
col2 <- rep(col[1],4)
col3 <- rep("blue",4)
pch2 <- rep(pch[1],4)



# calcul du I de Moran
#var <- var - mean(var)
wvar <- W%*%var

stdvar <- var/sd(var)
uns <- rep(1, length(var))
result <- nonormmoran(stdvar,uns,W)
MORAN <- result$morani
prob.I <- pnorm(result$istat)
rvar <- qr(var)
beta.I <- qr.coef(rvar,wvar)

# calcul de la variable obsq (pour les quadrants)

obsq[which((var > mean(var)) & (wvar >= mean(wvar)))] <- 1
obsq[which((var <= mean(var)) & (wvar > mean(wvar)))] <- 2
obsq[which((var < mean(var)) & (wvar <= mean(wvar)))] <- 3
obsq[which((var >= mean(var)) & (wvar < mean(wvar)))] <- 4
 
 # i de moran local
  x.centre<-(var-mean(var))
  wx.centre<-(W%*%x.centre)
  ilocal <- (x.centre/var(x.centre)) * (wx.centre)  

####################################################
# sélection d'un point
####################################################

pointfunc<-function() 
{
    quit <- FALSE
    quad <<- FALSE
 
    if(maptest) 
      { dev.set(2)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
        if(spdf & nrow(sp.obj)>75 & !buble) 
           {points(long,lat,pch=16,col='royalblue')}
      }
     else
     { dev.set(3)
       title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
       title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }
      
    while(!quit)
    {
        if(maptest) 
        {
            dev.set(2)
            loc<-locator(1)
            if(is.null(loc)) 
            {
              quit<-TRUE
              carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, carte=carte,nocart=nocart, classe=obsq,
              couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
              nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
              label=label, cex.lab=cex.lab, labmod=names.arg)   
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
        }
        else
        {
            dev.set(3)
            loc<-locator(1)
            if(is.null(loc)) 
            {
              quit<-TRUE
              graphique(var1=var, var2=wvar, var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
              symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
              legmap=legmap2, legends=legends2, bin=is.norm)
              next
            }
            obs<<-selectmap(var1=var,var2=wvar,obs=obs,Xpoly=loc[1], Ypoly=loc[2], method="point")
        }
        
        #graphiques
        carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=obsq,
        couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
        nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
        label=label, cex.lab=cex.lab, labmod=names.arg)      

        graphique(var1=var, var2=wvar, var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
        symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
        legmap=legmap2, legends=legends2, bin=is.norm)
      
      if(maptest) 
      { dev.set(2)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
        if(spdf & nrow(sp.obj)>75 & !buble) 
           {points(long,lat,pch=16,col='royalblue')}
      }
     else
     { dev.set(3)
       title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
       title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
     }
     
        
        if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        {
            graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)], 
            obs=obs, num=4, graph=graphChoice, couleurs=col[1], symbol=pch[1], labvar=c(varChoice1,varChoice2))
        }
    }
  }

####################################################
# sélection d'un point sur la carte
####################################################

pt1func <- function()
{
  ifelse(flower,method <<- "Neighbourplot1", method <<- "Quadrant")
  graph <<- "Moran" 
  maptest <<- TRUE  
  pointfunc()
}

####################################################
# sélection d'un point sur le graphique
####################################################

pt2func <- function()
{    
  ifelse(flower,method <<- "Neighbourplot1", method <<- "Quadrant")
  graph <<- "Moran" 
  maptest <<- FALSE  
  pointfunc()
  }

####################################################
# sélection d'un polygone
####################################################

polyfunc<-function() 
{
   polyX <- NULL
   polyY <- NULL
   quit <- FALSE
   quad <<- FALSE
   
   if(maptest) 
      { dev.set(2)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
        title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
      }
     else
     { dev.set(3)
       title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
       title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')
    }
      
   while(!quit)
    {
        if(maptest)
        {
          dev.set(2)
          if(spdf) 
          {points(long,lat,pch=16,col='royalblue')}
          loc<-locator(1)
          if(is.null(loc)) 
            {
              quit<-TRUE
              next
            }
        }
        else
        {
            dev.set(3);
            loc<-locator(1)
            if(is.null(loc)) 
            {
                quit<-TRUE
                next
            }   
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

    if (maptest)
    {
        obs <<- selectmap(var1=long, var2=lat, obs=obs, Xpoly=polyX, Ypoly=polyY, method="poly")
    }
    else
    {
        obs <<- selectmap(var1=var, var2=wvar, obs=obs, Xpoly=polyX, Ypoly=polyY, method="poly")
    }

    #graphiques
        carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=obsq,
        couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
        nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
        label=label, cex.lab=cex.lab, labmod=names.arg)        

        graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
        symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
        legmap=legmap2, legends=legends2, bin=is.norm)
        
     if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        {
          graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)], 
          obs=obs, num=4, graph=graphChoice, couleurs=col[1], symbol=pch[1], labvar=c(varChoice1,varChoice2))
        }
  }
}

####################################################
# sélection d'un polygone sur la carte
####################################################

poly1func <- function()
{
  ifelse(flower,method <<- "Neighbourplot1", method <<- "")
  graph <<- "Moran" 
       
    if (quad == TRUE)
    {
      SGfunc()
      quad <<- FALSE
    }
    
    maptest <<- TRUE
    polyfunc()
  }

####################################################
# sélection d'un polygone sur le graphique
####################################################

poly2func <- function()
{
 
  ifelse(flower,method <<- "Neighbourplot1", method <<- "")
  graph <<- "Moran" 
       
    if (quad == TRUE)
    {
      SGfunc()
      quad <<- FALSE
    }
    maptest <<- FALSE
    polyfunc()
  }

####################################################
# sélection d'un quadrant
####################################################

quadfunc <- function()
{
    if (quad == FALSE)
    {
      SGfunc()
      quad <<- TRUE
    }

    obs[which(obsq == num)] <<- !obs[which(obsq == num)]

    carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  classe=obsq,  carte=carte,nocart=nocart, labmod=names.arg, W=W,
    couleurs=col3,symbol=pch2, method=method,buble=buble,cbuble=z,criteria=criteria,nointer=nointer,legmap=legmap,
    legends=legends,axis=axes,lablong=lablong, lablat=lablat, label=label, cex.lab=cex.lab) 

    graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
    symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
    legmap=legmap2, legends=legends2, bin=is.norm)
    
        
            
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
     {
      graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
      obs=obs, num=4, graph=graphChoice, symbol=pch[1], couleurs=col[1], labvar=c(varChoice1,varChoice2))
     }
}

####################################################
# sélection d'un des 4 quadrants 
####################################################

quad1func <- function()
{
  num <<- 1 
  method <<- ifelse(flower,"Neighbourplot1","Quadrant") 
#  graph <<- "Quadrant"  
  quadfunc()
}

quad2func <- function()
{
 num <<- 2  
  method <<- ifelse(flower,"Neighbourplot1","Quadrant") 
 #graph <<- "Quadrant"      
 quadfunc()
}

quad3func <- function()
{
 num <<- 3
  method <<- ifelse(flower,"Neighbourplot1","Quadrant") 
# graph <<- "Quadrant"    
 quadfunc()
}

quad4func <- function()
{
 num <<- 4
  method <<- ifelse(flower,"Neighbourplot1","Quadrant") 
# graph <<- "Quadrant"    
 quadfunc()
}

####################################################
# Différentes couleurs selon le quadrant
####################################################

 colfunc <- function()
{
    if(!choix.col)
    {choix.col<<-TRUE
     method <<- "Quadrant"
     res1<-choix.couleur("Moran",col=col[1],pch=pch[1],legends=legends,spdf=spdf)     
     if(length(res1$col2)==4)
      {col2 <<- res1$col2
       col3 <<- col2 }
      else
      {col2 <<- rep(col[1],4)
       col3 <<- rep("blue",4) }
      
     ifelse(length(res1$pch2)==4,pch2 <<- res1$pch2,pch2 <<- rep(pch[1],4))
     legends <<- res1$legends
     }
    else
    {choix.col<<-FALSE
     col2 <<- rep(col[1],4)
     col3 <<- rep("blue",4)
     pch2 <<- rep(pch[1],4)
     legends <<- list(legends[[1]],FALSE,legends[[3]],"")
    }

carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=obsq,
couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
label=label, cex.lab=cex.lab,labmod=names.arg)     

graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
legmap=legmap2, legends=legends2, bin=is.norm)
        
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
            if(!(4%in%dev.list())) dev.new()
            graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
            obs=obs, num=4, graph=graphChoice, couleurs=col[1], symbol=pch[1], labvar=c(varChoice1,varChoice2))
        }
    }
    else
    {
        tkmessageBox(message="List of Variables and list of variables names must have been given",icon="warning",type="ok")
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
    carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=obsq,
    couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
    nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
    label=label, cex.lab=cex.lab,labmod=names.arg)        
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

carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj,  carte=carte,nocart=nocart, classe=obsq,
couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
label=label, cex.lab=cex.lab,labmod=names.arg)     

graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
legmap=legmap2, legends=legends2, bin=is.norm)
    
if ((graphChoice != "") && (varChoice1 != ""))
 {
 graphique(var1=listvar[,which(listnomvar == varChoice1)], var2=listvar[,which(listnomvar == varChoice2)],
 obs=obs, num=4, graph=graphChoice, couleurs=col[1], symbol=pch[1], labvar=c(varChoice1,varChoice2))
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
    res<-list(obs=which(obs),MORAN=MORAN)
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
 carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, carte=carte,nocart=nocart, classe=obsq,
 couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
 nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
 label=label, cex.lab=cex.lab,labmod=names.arg)       
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
  
carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, carte=carte,nocart=nocart, classe=obsq,
couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
label=label, cex.lab=cex.lab,labmod=names.arg)    
 
}

####################################################
# Bubble  of Lisa
####################################################

lisa<-function()
{        
  res3<-choix.bubble(buble2,abs(ilocal),"ilocal",legends2)
  
  buble2 <<- res3$buble
  legends2 <<- res3$legends
  z2 <<- res3$z
  legmap2 <<- res3$legmap
  
graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph=graph, labvar=labvar, couleurs=col2,
symbol=pch2, locmoran=locmoran,obsq=obsq, cex.lab=cex.lab,buble=buble2, cbuble=z2, 
legmap=legmap2, legends=legends2, bin=is.norm)
 
}

####################################################
# Permutation
####################################################
permutation<-function()
{
  tt1<-tktoplevel()
  Name <- tclVar("n")
  entry.Name <-tkentry(tt1,width="5",textvariable=Name)
  tkgrid(tklabel(tt1,text="Number of simulations"),entry.Name)

  OnOK <- function()
  { 
	 value1 <- tclvalue(Name)
	 tkdestroy(tt1)
       if(is.na(as.integer(value1)))
      {
        tkmessageBox(message="Sorry, but you have to choose decimal values",icon="warning",type="ok");
      }
      else
      {
       n=as.integer(value1)
       perm<-NULL
      for (i in 1:n)
       {
       sam<-sample(var,length(var))
       sam <- sam - mean(sam);

       epe <- sam %*% sam
       mi <- (sam %*% W %*% sam)/epe
       morani <- round(mi,4);
       perm <- c(perm,morani);
      }
     msg <- paste("The p-value of the permutation test is :",1-length(which(perm<rep(MORAN,n)))/n)
     tkmessageBox(message=msg,icon="info",type="ok") 
    }   
}


  OK.but <-tkbutton(tt1,text="   OK   ",command=OnOK)
  tkgrid(OK.but)
  tkfocus(tt1)
}

#################################################
###########   Representation

# Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2) # new environment
 {
   carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, carte=carte,nocart=nocart, classe=classe,
   couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
   nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
   label=label, cex.lab=cex.lab)     

   graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph="Moran", labvar=labvar, couleurs=col2, symbol=pch2,
   locmoran=locmoran,obsq=obsq,buble=buble2, cbuble=z2, legmap=legmap2, legends=legends2, bin=is.norm)
   assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {carte(long=long, lat=lat, obs=obs,  sp.obj=sp.obj, carte=carte,nocart=nocart, classe=classe,
    couleurs=col3, symbol=pch2, W=W, method=method, buble=buble, cbuble=z, criteria=criteria,
    nointer=nointer, legmap=legmap, legends=legends,axis=axes,lablong=lablong, lablat=lablat,
    label=label, cex.lab=cex.lab)     

   graphique(var1=var, var2=wvar,  var3=ilocal, obs=obs, num=3, graph="Moran", labvar=labvar, couleurs=col2, symbol=pch2,
   locmoran=locmoran,obsq=obsq,buble=buble2, cbuble=z2, legmap=legmap2, legends=legends2, bin=is.norm)
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
tkwm.title(tt, "moranplotmap")

frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
foreground = "darkred", background = "white"))
point.but <- tkbutton(frame1a, text="Selection by point", command=pt1func);
poly.but <- tkbutton(frame1a, text="Selection by polygon ", command=poly1func);
tkpack(point.but, poly.but, side = "left", expand = "TRUE",fill = "x")
tkpack(frame1a, expand = "TRUE", fill = "x")

frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1c, text = "Work on the graphic", font = "Times 12",
foreground = "darkred", background = "white"))
point2.but <- tkbutton(frame1c, text="Selection by point", command=pt2func);
poly2.but <- tkbutton(frame1c, text="Selection by polygon ", command=poly2func);
tkpack(point2.but, poly2.but, side = "left", expand = "TRUE",fill = "x")
tkpack(frame1c, expand = "TRUE", fill = "x")

frame1d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1d, text = "Select a quadrant", font = "Times 12",
foreground = "darkred", background = "white"))
HH.but <- tkbutton(frame1d, text="H.-H.", command=quad1func);
LH.but <- tkbutton(frame1d, text="H.-L.", command=quad4func);
LL.but <- tkbutton(frame1d, text="L.-L.", command=quad3func);
HL.but <- tkbutton(frame1d, text="L.-H.", command=quad2func);
tkpack(HH.but,LH.but,LL.but,HL.but, side = "left", expand = "TRUE",fill = "x")
tkpack(frame1d, expand = "TRUE", fill = "x")

frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nettoy.but <- tkbutton(frame1b, text="     Reset selection     " , command=SGfunc);
tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame1b, expand = "TRUE", fill = "x")


frame11a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame11a, text = "Moran test", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(frame11a, expand = "TRUE", fill = "x")

frame11b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
msg <- paste("Moran index ",ifelse(is.norm,"(W normalized)","(W not normalized)"),": ", MORAN, " - ","p-value (Gaussian Test) : ", ifelse(round(1-prob.I,4)<0.0001,"<0.0001",round(1-prob.I,4)))
tkgrid(tklabel(frame11b,text=msg),columnspan=2)
tkpack(frame11b, expand = "TRUE", fill = "x")

frame11c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
noint10.but <- tkbutton(frame11c, text="Permutation Test", command=permutation);
tkpack(noint10.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame11c, expand = "TRUE", fill = "x")



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
tkpack(tklabel(frame2e, text = "Print different colors by quadrant  ", font = "Times 10",
foreground = "darkred", background = "white"),tklabel(frame2e, text = "Bubles of LISA", font = "Times 10",
foreground = "darkred", background = "white"), side = "left", fill="x",expand = "TRUE")
tkpack(frame2e, expand = "TRUE", fill = "x")

frame2f <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
colquad.but <- tkbutton(frame2f, text="On/Off", command=colfunc)
lisa.but <- tkbutton(frame2f, text="     OK     " , command=lisa)
tkpack(colquad.but,lisa.but, side = "left", expand = "TRUE", fill = "x")
tkpack(frame2f, expand = "TRUE", fill = "x")

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

