`angleplotmap` <- function(sp.obj, name.var, quantiles=TRUE,
names.attr=names(sp.obj), criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16, col="lightblue3",
xlab="angle", ylab="absolute magnitude", axes=FALSE, lablong="", lablat="")
{
envir = as.environment(1)
# Verification of the Spatial Object sp.obj
class.obj<-class(sp.obj)[1]

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

 if(is.null(carte) & class.obj=="SpatialPolygonsDataFrame") carte<-spdf2list(sp.obj)$poly

# for identifying the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")

 # initialisation
  nointer<-FALSE
  nocart<-FALSE
  buble<-FALSE
  legends<-list(FALSE,FALSE,"","")
  z<-NULL
  legmap<-NULL
  inout=NULL
  labvar=c(xlab,ylab)
  obs <- matrix(FALSE, nrow = length(long), ncol = length(long))
  names.slide="Alpha Quantile Value"
# Transformation d'un data.frame en matrix

if((length(listvar)>0) && (dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)


# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()


# calcul des matrices theta et absvar

    long1 <- matrix(rep(t(long), length(long)), ncol = dim(t(long))[2],byrow = FALSE)
    long2 <- matrix(rep(t(long), length(long)), ncol = dim(t(long))[2],byrow = TRUE)
    lat1 <- matrix(rep(t(lat), length(lat)), ncol = dim(t(lat))[2],byrow = FALSE)
    lat2 <- matrix(rep(t(lat), length(lat)), ncol = dim(t(lat))[2],byrow = TRUE)
    theta <- matrix(0, nrow = length(long), ncol = length(long))
    numer <- lat2 - lat1
    denom <- long2 - long1
   

   for (i in 1:length(long)) 
       {
        for (j in 1:length(long)) 
           {
             if (denom[i, j] == 0) 
               {
                theta[i, j] <- pi/2
               }
             else 
               {
                theta[i, j] <- atan(numer[i, j]/denom[i, j])
               }
           }
       }

    theta[which(theta < 0)] <- theta[which(theta < 0)] + pi
    v1 <- matrix(rep(t(var), length(var)), ncol = dim(t(var))[2],byrow = FALSE)
    v2 <- matrix(rep(t(var), length(var)), ncol = dim(t(var))[2],byrow = TRUE)
    absvar <- abs(v1 - v2)
   

####################################################
#choix des bornes des réglettes (inspiré de la documentation de Matlab)
####################################################


#        u1 <- sort(theta)
#        u1 <- as.vector(u1)
#        z <- seq(1, max(u1), by = (max(u1)/3500))
#        z <- round(z)
#        z1 <- z[2:length(z)] - z[1:(length(z) - 1)]
#        h <- mean(z1)
#        p <- 1/(1 + (h^3/6))
#        p1 <- 1/(1 + (h^3/60))
#        p2 <- 1/(1 + (h^3/0.6))
#        alpha <- (1 - p)/p
#        borne1 <- (1 - p1)/p1
#        borne2 <- (1 - p2)/p2
#

borne1=0.01
borne2=0.99
alpha=0.5
####################################################
# sélection d'un point sur l'angleplot
####################################################

 pointfunc <- function() 
 {
  quit <- FALSE
 
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
        graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "Angleplot", labvar = labvar,
        couleurs=col,symbol = pch, quantiles = quantiles,alpha1 = alpha)
        next
       }

   obs <<- selectstat(var1 = theta, var2 = absvar, obs = obs,Xpoly = loc[1], Ypoly = loc[2],
   method = "AnglePoint",long = long, lat = lat)
   
   diag(obs)<<-FALSE
   
   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "Angleplot", labvar = labvar,
   couleurs=col,symbol = pch, quantiles = quantiles,alpha1 = alpha)
   title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
   title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "Angleplot",
   axis=axes,legmap=legmap,legends=legends) 
 }
}



####################################################
# sélection d'un polygone sur l'angleplot
####################################################


 polyfunc <- function() 
 {
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

      if (length(polyX)>0)
      {
        lines(polyX, polyY)
     }
    }  
     
   polyX <- c(polyX, polyX[1])
   polyY <- c(polyY, polyY[1])

  if (length(polyX)>0)
    {
     lines(polyX, polyY)
        for (i in 1:length(long)) 
         {
          def <- inout(cbind(theta[, i], absvar[, i]),cbind(polyX, polyY), bound = TRUE)       
          obs[def,i ] <<- !obs[def,i ]    
         }

   graphique(var1 = theta, var2 = absvar, obs = obs, num = 3,graph = "Angleplot", labvar = labvar,couleurs=col,
   symbol = pch,quantiles = quantiles, alpha1 = alpha)
   
   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,
   method = "Angleplot",axis=axes,legmap=legmap,legends=legends) 

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
    carte(long=long, lat=lat,criteria=criteria,buble=buble,cbuble=z,nointer=nointer,obs=obs,
    lablong=lablong, lablat=lablat,method="Angleplot",label=label,cex.lab=cex.lab, symbol=pch,
    carte=carte,nocart=nocart,axis=axes,legmap=legmap,legends=legends) 
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
    obs<<-matrix(FALSE,nrow=length(long), ncol=length(long));
    
    carte(long=long, lat=lat,criteria=criteria,buble=buble,cbuble=z,nointer=nointer,obs=obs,
    lablong=lablong, lablat=lablat,method="Angleplot",label=label,cex.lab=cex.lab, symbol=pch,
    carte=carte,nocart=nocart,axis=axes,legmap=legmap,legends=legends) 
    
    graphique(var1=theta, var2=absvar, obs=obs, num=3, graph="Angleplot", labvar=labvar,
    couleurs=col,symbol=pch,  quantiles=quantiles, alpha1=alpha);
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
    obs[lower.tri(obs)]<-FALSE
    assign("last.select", which(obs,arr.ind=TRUE), envir = envir)
}


####################################################
# Open a no interactive selection
####################################################

fnointer<-function() 
{
 if (length(criteria) != 0)
 {
    ifelse(!nointer,nointer<<-TRUE,nointer<<-FALSE)
    carte(long=long, lat=lat,criteria=criteria,buble=buble,cbuble=z,nointer=nointer,obs=obs,
    lablong=lablong, lablat=lablat,method="Angleplot",label=label,cex.lab=cex.lab, symbol=pch,
    carte=carte,nocart=nocart,axis=axes,legmap=legmap,legends=legends) 
   
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
  
  carte(long=long, lat=lat,criteria=criteria,buble=buble,cbuble=z,nointer=nointer,obs=obs,
  lablong=lablong, lablat=lablat,method="Angleplot",label=label,cex.lab=cex.lab, symbol=pch,
  carte=carte,nocart=nocart,axis=axes,legmap=legmap,legends=legends) 

}


####################################################
# Pour le alpha 
####################################################

    refresh.code <- function(...)
     {
        alpha <<- slider1(names.slide=names.slide,no=1)
        graphique(var1 = theta, var2 = absvar, obs = obs, num = 3, 
        graph = "Angleplot",couleurs=col, labvar = labvar, symbol = pch, 
        quantiles = quantiles, alpha1 = alpha)
    }

####################################################
# Représentation graphique
####################################################
# Is there a Tk window already open ?
if(interactive())
{
 if(!exists("GeoXp.open",envir = baseenv())||length(ls(envir=.TkRoot$env, all.names=TRUE))==2)  # new environment
 {
  carte(long=long, lat=lat, obs=obs, lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, 
  symbol=pch,carte=carte,nocart=nocart, method="Angleplot",axis=axes,legends=legends)

  graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "Angleplot", labvar = labvar,
  couleurs=col,symbol = pch, quantiles = quantiles,alpha1 = alpha)
  assign("GeoXp.open", TRUE, envir = baseenv())
 }
 else
 {if(get("GeoXp.open",envir= baseenv()))
   {stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")}
  else
  {carte(long=long, lat=lat, obs=obs, lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, 
  symbol=pch,carte=carte,nocart=nocart, method="Angleplot",axis=axes,legends=legends)

  graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "Angleplot", labvar = labvar,
  couleurs=col,symbol = pch, quantiles = quantiles,alpha1 = alpha)
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
tkwm.title(tt, "angleplotmap")

frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame1a, text = "Work on the graph", font = "Times 12",
foreground = "darkred", background = "white"))

point.but <- tkbutton(frame1a, text="Selection by point", command=pointfunc);
poly.but <- tkbutton(frame1a, text="Selection by polygon ", command=polyfunc);
tkpack(point.but, poly.but, side = "left", expand = "TRUE",fill = "x")

tkpack(frame1a, expand = "TRUE", fill = "x")

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

if(quantiles)
{
frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")

slider1(frame1c, refresh.code, names.slide=names.slide,
        borne1, borne2, (borne2 - borne1)/100, alpha)

tkpack(frame1c, expand = "TRUE", fill = "x")
}

frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
foreground = "blue", background = "white"))

quit.but <- tkbutton(frame3, text="Save results", command=quitfunc2);
quit.but2 <- tkbutton(frame3, text="Exit without saving", command=quitfunc);

tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame3, expand = "TRUE", fill = "x")
}
####################################################
return(invisible())

# return(obs)
}

