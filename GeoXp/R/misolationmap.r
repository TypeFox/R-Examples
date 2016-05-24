`misolationmap` <- function(sp.obj, nb.obj, names.var, propneighb=0.4,chisqqu=0.975,
names.attr=names(sp.obj), criteria=NULL, carte=NULL, identify=FALSE, cex.lab=0.8, pch=16, col="lightblue3",
xlab="degree of isolation", ylab="Pairwise Mahalanobis distances", axes=FALSE, lablong="", lablat="")
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

 # for identifyng the selected sites
ifelse(identify, label<-row.names(listvar),label<-"")


 # initialisation
  xy<-cbind(long, lat)
  nointer<-FALSE
  nocart<-FALSE
  buble<-FALSE
  legends<-list(FALSE,FALSE,"","")
  z<-NULL
  legmap<-NULL
  inout=NULL
  labvar=c(xlab,ylab)
  obs <- matrix(FALSE, nrow = length(long), ncol = length(long))
    graf<-"Neighbourplot1"
  obsref<-obs
  
  outselect=TRUE
  outselect2=TRUE
  orderselect=TRUE

   W<-nb2mat(nb.obj)
   Wref<-W
   
# Transformation d'un data.frame en matrix
if((length(listvar)>0) && (dim(as.matrix(listvar))[2]==1)) listvar<-as.matrix(listvar)
if((length(dataset)>0) && (dim(as.matrix(dataset))[2]==1)) dataset<-as.matrix(dataset)

# Windows device
if(!(2%in%dev.list())) dev.new()
if(!(3%in%dev.list())) dev.new()

# calcul des matrices theta et absvar

n=nrow(dataset)
p=ncol(dataset)

covr=covMcd(dataset,alpha=0.75)
cinv=solve(covr$cov)
MDglobal=sqrt(mahalanobis(dataset, covr$center, cinv, inverted=TRUE))


# TRUE/FALSE for non-outlying/outlying:
qchi=sqrt(qchisq(chisqqu,p))
MDglobalTF <- (MDglobal<qchi)

idx=matrix(1:n,n,n)
se=as.vector(idx[lower.tri(idx)])
hlp=as.matrix(dataset[rep(1:(n-1),seq((n-1),1)),]-dataset[se,])
MDij=sqrt(rowSums((hlp%*%cinv)*hlp))

MDpair=matrix(0,n,n)
MDpair[lower.tri(MDpair)] <- MDij
MDpair=t(MDpair)
MDpair[lower.tri(MDpair)] <- MDij

MDpairN=vector("list", n)

# boundary that should include required proportion of neighbors:
chibound=rep(NA,n)
theta<-matrix(0,n,n)
absvar<-matrix(0,n,n)

for (i in 1:n){
  MDpairN[[i]] <- MDpair[i,nb.obj[[i]]]
  nn=max(1,round(length(MDpairN[[i]])*propneighb)) # number of neighbors in tolerance ellipse
  nval=MDpairN[[i]][order(MDpairN[[i]])][nn]  # value of largest neighbor to be included
  chibound[i]=pchisq(nval^2,p,MDglobal[i]^2)
  
  theta[i,nb.obj[[i]]]<-chibound[i]
  absvar[i,nb.obj[[i]]]<-MDpair[i,nb.obj[[i]]]

}

thetaref=theta
absvaref=absvar

thetaref2=thetaref
absvaref2=absvaref

# sort according to values of "chibound" - separately for outliers and non-outliers
idx1=order(chibound[MDglobalTF])
idx0=order(chibound[!MDglobalTF])

idxg<-sort(chibound,index.return=TRUE)

# calcul des distances de Mahalanobis par site

rd <- sqrt(mahalanobis(dataset, center = covr$center, cov = covr$cov))

pcrit<-ifelse(p <= 10,(0.24 - 0.003 * p)/sqrt(n), (0.252 - 0.0018 * p)/sqrt(n))
 delta <- qchisq(1 - 0.025, p)
 
 d2 <- mahalanobis(dataset, covr$center, covr$cov)
 d2ord <- sort(d2)
 dif <- pchisq(d2ord, p) - (0.5:n)/n
 i <- (d2ord >= delta) & (dif > 0)
 alfan<-ifelse(sum(i) == 0,0,max(dif[i]))
 if (alfan < pcrit) 
        alfan <- 0
 cn<-ifelse(alfan > 0, max(d2ord[n - ceiling(n * alfan)], delta), Inf)

alphab<-ifelse(cn != Inf, sqrt(c(cn, qchisq(c(0.75, 0.5, 0.25), ncol(dataset)))),sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(dataset))))


chi2.quant<-rep(0,n)
lalpha <- length(alphab)
    for (j in 1:lalpha) {
            if (j == 1) {
                (chi2.quant[which(rd >= alphab[j])]<-lalpha)
            }
            else {
                chi2.quant[which((rd < alphab[j - 1]) & (rd >= alphab[j]))]<-lalpha+1-j

            }
     }



####################################################
# sélection d'un point sur la carte
####################################################

pointfunca<-function()
{
   if (graf=="pairwise") SGfunc()
   graf<<-"Neighbourplot1"
   
    quit <- FALSE
    dev.set(2)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
    title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

    while(!quit)
    {
        dev.set(2)
        loc<-locator(1)
        if (is.null(loc))
        {
            quit<-TRUE
            carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
            nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
            axis=axes,legmap=legmap,legends=legends)
            next
        }
        obs2<-selectmap(var1=long,var2=lat,obs=obs,Xpoly=loc[1], Ypoly=loc[2], method="point");

       obs<<-(W*obs2>0)

       if(!outselect)
        {
         obs[MDglobalTF,]<<-FALSE
        }

    if(!outselect2)
     {
      obs[!MDglobalTF,]<<-FALSE
     }
     
    #    diag(obs)<<-FALSE
        # graphiques
   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
   axis=axes,legmap=legmap,legends=legends)
   title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
   title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

   #     carte(long=long, lat=lat, obs=obs, lablong=lablong, lablat=lablat, label=label, symbol=16,
   #     method="Neighbourplot1", W=W,axis=axes,legmap=legmap,legends=legends,buble=buble,criteria=criteria,
   #     nointer=nointer,cbuble=z,carte=carte,nocart=nocart,couleurs="blue",classe=card(object),cex.lab=cex.lab)


   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)

  #  obs <<- matrix(FALSE, nrow=length(long), ncol=length(long));

    }
  }


####################################################
# sélection d'un polygone
####################################################

polyfunca<-function()
{
   if (graf=="pairwise") SGfunc()
   graf<<-"Neighbourplot1"
   
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

    obs2 <- selectmap(var1=long, var2=lat, obs=obs, Xpoly=polyX, Ypoly=polyY, method="poly")
    obs<<-(W*obs2>0)

    if(!outselect)
     {
      obs[MDglobalTF,]<<-FALSE
     }
     
     if(!outselect2)
     {
      obs[!MDglobalTF,]<<-FALSE
     }

   #    diag(obs)<<-FALSE
       
    # graphiques
   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
   axis=axes,legmap=legmap,legends=legends)

   #     carte(long=long, lat=lat, obs=obs, lablong=lablong, lablat=lablat, label=label, symbol=16,
   #     method="Neighbourplot1", W=W,axis=axes,legmap=legmap,legends=legends,buble=buble,criteria=criteria,
   #     nointer=nointer,cbuble=z,carte=carte,nocart=nocart,couleurs="blue",classe=card(object),cex.lab=cex.lab)


   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)

 #   obs <<- matrix(FALSE, nrow=length(long), ncol=length(long));

}
  }


####################################################
# sélection d'un point sur l'angleplot
####################################################

   
 pointfunc <- function() 
 {
   if (graf=="Neighbourplot1") SGfunc()
   graf<<-"pairwise"
   
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
        graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
        couleurs=col,symbol = pch, direct=propneighb)
        next
       }

   obs <<- selectstat(var1 = theta, var2 = absvar, obs = obs,Xpoly = loc[1], Ypoly = loc[2],
   method = "AnglePoint",long = long, lat = lat)
   
    if(!outselect)
     {
      obs[MDglobalTF,]<<-FALSE
     }

    if(!outselect2)
     {
      obs[!MDglobalTF,]<<-FALSE
     }
     
   diag(obs)<<-FALSE
   
   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)
   title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')
   title(sub = "To stop selection, click on the right button of the mouse and stop (for MAC, ESC)", cex.sub = 0.8, font.sub = 3,col.sub='red')

   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
   axis=axes,legmap=legmap,legends=legends) 
 }
}


####################################################
# sélection des outliers
####################################################


 outlier <- function()
 {

   outselect<<-!outselect
   labvar=c(xlab,ylab)


   if (!outselect) tkmessageBox(message="You have selected only global outliers")
   if (outselect) tkmessageBox(message="You have remited non global outliers")

   if(outselect)
   {
    theta<<-thetaref
    absvar<<-absvaref
    thetaref2<<-thetaref
    absvaref2<<-absvaref
      if(!orderselect)
       {
       labvar=c("Rank",ylab)
       ind0<-which(theta==0,arr.ind=TRUE)
       theta[idxg$ix,]<<-t(matrix(rep(1:n,each=length(idxg$ix)),length(idxg$ix),length(idxg$ix)))
       theta[ind0]<<-0
       theta<<-theta
       }
    }
   else
   {
   SGfunc()
   if(!outselect2)
    {
    theta<<-thetaref
    absvar<<-absvaref
      if(!orderselect)
       {
       labvar=c("Rank",ylab)
       ind0<-which(theta==0,arr.ind=TRUE)
       theta[idxg$ix,]<<-t(matrix(rep(1:n,each=length(idxg$ix)),length(idxg$ix),length(idxg$ix)))
       theta[ind0]<<-0
       theta<<-theta
       }
    }
    #theta<<-thetaref[!MDglobalTF,]
    #absvar<<-absvaref[!MDglobalTF,]
    theta[MDglobalTF,]<<-0
    absvar[MDglobalTF,]<<-0

    thetaref2<<-theta
    absvaref2<<-absvar
   # obs[MDglobalTF,]<<-FALSE
   # W[MDglobalTF,]<<-0
       }

   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)

   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
   axis=axes,legmap=legmap,legends=legends)

}

####################################################
# sélection des non outliers
####################################################


 outlier2 <- function()
 {

   outselect2<<-!outselect2
   labvar=c(xlab,ylab)

   if (!outselect) tkmessageBox(message="You have selected non global outliers")
   if (outselect) tkmessageBox(message="You have remited global outliers")

   if(outselect2)
   {
    theta<<-thetaref
    absvar<<-absvaref
    thetaref2<<-thetaref
    absvaref2<<-absvaref
       if(!orderselect)
       {
       labvar=c("rank",ylab)
       ind0<-which(theta==0,arr.ind=TRUE)
       theta[idxg$ix,]<<-t(matrix(rep(1:n,each=length(idxg$ix)),length(idxg$ix),length(idxg$ix)))
       theta[ind0]<<-0
       theta<<-theta
       }
    }
   else
   {
   SGfunc()
    #theta<<-thetaref[!MDglobalTF,]
    #absvar<<-absvaref[!MDglobalTF,]
    if(!outselect)
    {
    theta<<-thetaref
    absvar<<-absvaref
      if(!orderselect)
       {
       labvar=c("Rank",ylab)
       ind0<-which(theta==0,arr.ind=TRUE)
       theta[idxg$ix,]<<-t(matrix(rep(1:n,each=length(idxg$ix)),length(idxg$ix),length(idxg$ix)))
       theta[ind0]<<-0
       theta<<-theta
       }
    }
    theta[!MDglobalTF,]<<-0
    absvar[!MDglobalTF,]<<-0
    
    thetaref2<<-theta
    absvaref2<<-absvar
   # obs[MDglobalTF,]<<-FALSE
   # W[MDglobalTF,]<<-0
       }

   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)

   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
   axis=axes,legmap=legmap,legends=legends)

}

####################################################
# trier les observations
####################################################


 ordering <- function()
 {
   orderselect<<-!orderselect

   if (orderselect) tkmessageBox(message="Degree of isolation on the x-axis")
   if (!orderselect) tkmessageBox(message="Rank on the x-axis")


   if(orderselect)
   {
    theta<<-thetaref
    absvar<<-absvaref
  #  obs<<-obsref
  #  W<<-Wref
    }
   else
   {
    labvar=c("Rank",ylab)
    ind0<-which(theta==0,arr.ind=TRUE)
    theta[idxg$ix,]<<-t(matrix(rep(1:n,each=length(idxg$ix)),length(idxg$ix),length(idxg$ix)))

    theta[ind0]<<-0
  #   print(idxg$ix)
   #  print(theta[1:10,1:10])
  #   theta <- t(theta)
   #  print(theta[1:10,1:10])
     theta<<-theta

    #print(absvar[1:10,1:10])
    #print(idxg$ix)
    #print(absvar[1:10,1:10])
    
#    long<<-long[idxg$ix]
#    lat<<-lat[idxg$ix]
   # obs[MDglobalTF,]<<-FALSE
   # W[MDglobalTF,]<<-0
       }

   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)

   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,method = "pairwise",
   axis=axes,legmap=legmap,legends=legends)

}
####################################################
# sélection d'un polygone sur l'angleplot
####################################################


   
 polyfunc <- function() 
 {
   if (graf=="Neighbourplot1") SGfunc()
   graf<<-"pairwise"
   
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

    if(!outselect)
     {
      obs[MDglobalTF,]<<-FALSE
     }

    if(!outselect2)
     {
      obs[!MDglobalTF,]<<-FALSE
     }

    #  diag(obs)<<-FALSE
      
   graphique(var1 = theta, var2 = absvar, obs = obs, num = 3,graph = "pairwise", labvar = labvar,couleurs=col,
   symbol = pch, direct=propneighb)
   
   carte(long = long, lat = lat, obs = obs,buble=buble,criteria=criteria,nointer=nointer,cbuble=z,carte=carte,
   nocart=nocart, lablong = lablong,lablat = lablat,label = label,cex.lab=cex.lab, symbol = pch,
   method = "pairwise",axis=axes,legmap=legmap,legends=legends)

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
    lablong=lablong, lablat=lablat,method="pairwise",label=label,cex.lab=cex.lab, symbol=pch,
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
    lablong=lablong, lablat=lablat,method="pairwise",label=label,cex.lab=cex.lab, symbol=pch,
    carte=carte,nocart=nocart,axis=axes,legmap=legmap,legends=legends) 
    
    graphique(var1=theta, var2=absvar, obs=obs, num=3, graph="pairwise", labvar=labvar,
    couleurs=col,symbol=pch, direct=propneighb)
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
    lablong=lablong, lablat=lablat,method="pairwise",label=label,cex.lab=cex.lab, symbol=pch,
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
  res2<-choix.bubble(buble,cbind(chi2.quant,listvar),c("chi2.quant",listnomvar),legends)
  
  buble <<- res2$buble
  legends <<- res2$legends
  z <<- res2$z
  legmap <<- res2$legmap

  if(legends[[1]])
  {if((legmap[length(legmap)]=="chi2.quant"))
   {
    legmap<<-c(legmap,paste(">",round(alphab[1],2)),paste(round(alphab[2],2),"-",round(alphab[1],2)),
    paste(round(alphab[3],2),"-",round(alphab[2],2)),paste(round(alphab[4],2),"-",round(alphab[3],2)),
    paste("<",round(alphab[4],2)),"Mahalanobis")
   }
  }

  carte(long=long, lat=lat,criteria=criteria,buble=buble,cbuble=z,nointer=nointer,obs=obs,
  lablong=lablong, lablat=lablat,method="pairwise",label=label,cex.lab=cex.lab, symbol=pch,
  carte=carte,nocart=nocart,axis=axes,legmap=legmap,legends=legends) 

}



####################################################
# Représentation graphique
####################################################

carte(long=long, lat=lat, obs=obs, lablong=lablong, lablat=lablat, label=label,cex.lab=cex.lab, 
      symbol=pch,method="pairwise",axis=axes,legends=legends)

   graphique(var1 = theta, var2 = absvar, obs = obs,num = 3, graph = "pairwise", labvar = labvar,
   couleurs=col,symbol = pch, direct=propneighb)

####################################################
# création de la boite de dialogue
####################################################

if(interactive())
{
fontheading<-tkfont.create(family="times",size=14,weight="bold")

tt <- tktoplevel()
tkwm.title(tt, "misolationmap")

frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
foreground = "blue", background = "white"))
tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
foreground = "darkred", background = "white"))
point.but <- tkbutton(frame1a, text="Selection by point", command=pointfunca);
poly.but <- tkbutton(frame1a, text="Selection by polygon", command=polyfunca);
tkpack(point.but, poly.but, side = "left", expand = "TRUE",
        fill = "x")

tkpack(frame1a, expand = "TRUE", fill = "x")

frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame1c, text = "Work on the graphic", font = "Times 12",
foreground = "darkred", background = "white"))
intervalle1.but <- tkbutton(frame1c, text="Selection by point", command=pointfunc);
intervalle11.but <- tkbutton(frame1c, text="Selection by polygon", command=polyfunc);
tkpack(intervalle1.but,intervalle11.but, side = "left", expand = "TRUE", fill = "x")
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
foreground = "darkred", background = "white"),tklabel(frame2, text = "Bubbles", font = "Times 11",
foreground = "darkred", background = "white"),side = "left", fill="x",expand = "TRUE")
tkpack(frame2, expand = "TRUE", fill = "x")

frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nocou1.but <- tkbutton(frame2b, text="On/Off", command=cartfunc)
noint1.but <- tkbutton(frame2b, text="On/Off", command=fnointer)
bubble.but <- tkbutton(frame2b, text="On/Off", command=fbubble)
tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE",
        fill = "x")
tkpack(frame2b, expand = "TRUE", fill = "x")


frame2c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
tkpack(tklabel(frame2c, text = "Only global outliers ", font = "Times 11",
foreground = "darkred", background = "white"),tklabel(frame2c, text = "Without global outliers ", font = "Times 11",
foreground = "darkred", background = "white"),tklabel(frame2c, text = "Rank on the x-axis", font = "Times 11",
foreground = "darkred", background = "white"),side = "left", fill="x",expand = "TRUE")
tkpack(frame2c, expand = "TRUE", fill = "x")

frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
nocou1.but <- tkbutton(frame2b, text="On/Off", command=outlier)
noint1.but <- tkbutton(frame2b, text="On/Off", command=outlier2)
bubble.but <- tkbutton(frame2b, text="On/Off", command=ordering)
tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE",
        fill = "x")
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
####################################################

return(invisible())
}

