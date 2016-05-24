##
## Summary Statistics another built in
##
spT.MCMC.stat<-function(x, nBurn=0)
{
  options(warn=-1)
  if(class(x) != "spT"){
    stop("\n# Error: provide valid posterior output \n")
  }
  model<-x$model
  nItr<-x$iterations
  if(!is.null(nBurn)){
    if(nBurn == 0){
    nBurn<-0
    nItr<-nItr-x$nBurn
    }
    if(nBurn != 0){ 
    nBurn<-nBurn
    nItr<-nItr-x$nBurn
    }
  }
  if((nBurn+x$nBurn) >= (nItr+x$nBurn)){
   cat("# Number of Iterations:            ", nItr+x$nBurn, "\n")
   cat("# Number of Burn-in (fitted model):", x$nBurn, "\n")
   cat("# More Burn-in:                    ", nBurn, "\n")
   cat("# Total Number of Burn-in:         ", nBurn+x$nBurn, "\n")
   stop("\n# Error: iteration (",nItr+x$nBurn,") is less than or equal to total burn-in (",nBurn+x$nBurn,") \n")
  }
  cat("\n")
  cat("# Model:", model, "\n")
  if(is.null(model)==TRUE){
   stop("\n# Error: need to define the model")
  }
  else if(model=="AR"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-dim(x$mu_lp)[[1]]
  p<-dim(x$betap)[[1]]
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  }
  #
  para<-spT.Summary.Stat(para)
  dimnames(para)[[1]][1:(1+p+3)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  round(para,4)
  }
  else if(model == "GPP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  }
  #
  para<-spT.Summary.Stat(para)
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  round(para,4)
  }
  else if(model == "GP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  }
  #
  para<-spT.Summary.Stat(para)
  dimnames(para)[[1]]<-c(dimnames(x$X)[[2]],
     "sig2eps","sig2eta","phi")
  round(para,4)
  }
  else{

  }
}
##
## MCMC plot another built in 
##
spT.MCMC.plot<-function(x, nBurn=0, ACF="FALSE", PARTIAL.acf="FALSE")
{
  options(warn=-1)
  if(class(x) != "spT"){
    stop("\n# Error: provide valid posterior output \n")
  }
  model<-x$model
  nItr<-x$iterations
  if(!is.null(nBurn)){
    if(nBurn == 0){
    nBurn<-0
    nItr<-nItr-x$nBurn
    }
    if(nBurn != 0){ 
    nBurn<-nBurn
    nItr<-nItr-x$nBurn
    }
  }
  if((nBurn+x$nBurn) >= (nItr+x$nBurn)){
   cat("# Number of Iterations:            ", nItr+x$nBurn, "\n")
   cat("# Number of Burn-in (fitted model):", x$nBurn, "\n")
   cat("# More Burn-in:                    ", nBurn, "\n")
   cat("# Total Number of Burn-in:         ", nBurn+x$nBurn, "\n")
   stop("\n# Error: iteration (",nItr+x$nBurn,") is less than or equal to total burn-in (",nBurn+x$nBurn,") \n")
  }
  cat("\n")
  cat("# Model:", model, "\n")
  if(is.null(model)==TRUE){
   stop("\n# Error: need to define the model")
  }
  else if(model=="AR"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi","nu")
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  }
  #
  #
  }
  else if(model == "GPP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-dim(x$mu_lp)[[1]]
  p<-dim(x$betap)[[1]]
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi","nu")
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  }
  #
  }
  else if(model == "GP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],
     "sig2eps","sig2eta","phi","nu")
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],
     "sig2eps","sig2eta","phi")
  }
  #
  #
  }
  else{

  }
  #
  ##
   dev.new()
   for(i in 1:dim(para)[[1]]){
    MCMC.plot.obj(para[i,],name=c(dimnames(para)[[1]][i]),ACF=ACF,PARTIAL.acf=PARTIAL.acf)
    par(ask=TRUE)
   }
}
##
## segment plot for upper and lower limit
##  
spT.segment.plot<-function(obs, est, up, low, limit=NULL){
  #
  tmp<-cbind(obs,est,up,low)
  tmp<-na.omit(tmp)
  if(is.null(limit)==TRUE){
  plot(tmp[,1],tmp[,2],xlab="Observations",ylab="Predictions", pch="*",
  xlim=c(min(c(tmp),na.rm=TRUE),max(c(tmp),na.rm=TRUE)),
  ylim=c(min(c(tmp),na.rm=TRUE),max(c(tmp),na.rm=TRUE)))
  }
  #
  else{
  plot(tmp[,1],tmp[,2],xlab="Observations",ylab="Predictions",
   xlim=c(limit[1],limit[2]),ylim=c(limit[1],limit[2]),pch="*")
  }
  #
  segments(tmp[,1],tmp[,2],tmp[,1],tmp[,3])
  segments(tmp[,1],tmp[,2],tmp[,1],tmp[,4])
  #  
}
##
## hit and false alarm function for forecast
##
 spT.hit.false<-function(obs,fore,tol){
   #
   tmp<-cbind(obs,fore)
   tmp<-na.omit(tmp)
   #
   c11<-tmp[tmp[,1]<=tol & tmp[,2]<=tol,]
   c11<-length(c11)/2
   c12<-tmp[tmp[,1]<=tol & tmp[,2]>tol,]
   c12<-length(c12)/2  
   c21<-tmp[tmp[,1]>tol & tmp[,2]<=tol,]
   c21<-length(c21)/2  
   c22<-tmp[tmp[,1]>tol & tmp[,2]>tol,]
   c22<-length(c22)/2
   mat<-matrix(c(c11,c21,c12,c22),2,2)
   dimnames(mat)[[1]]<-c(paste("[Obs:<=",tol,"]"),paste("[Obs:> ",tol,"]")) 
   dimnames(mat)[[2]]<-c(paste("[Forecast:<=",tol,"]"),paste("[Forecast:>",tol,"]")) 
   POD<-round(mat[1,1]/sum(diag(mat)),4)
   FAR<-round(mat[2,1]/sum(mat[,1]),4)
   HAR<-round(mat[1,1]/sum(mat[1,]),4)
   #FAR<-round(mat[1,2]/sum(mat[1,]),4)
   #HAR<-round(sum(diag(mat))/sum(mat),4)
   top<-2*(mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1])
   bot<-mat[1,2]^2+mat[2,1]^2+2*mat[1,1]*mat[2,2]+(mat[1,2]+mat[2,1])*sum(diag(mat))
   S<-round(top/bot,4)
   x<-list(False.Alarm=FAR,Hit.Rate=HAR,Probability.of.Detection=POD,
      Heidke.Skill=S,cross.table=mat,tolerance.limit=tol) 
   x
}
##
## For data split
##
spT.subset<-function (data, var.name, s = NULL, reverse = FALSE) 
{
#
# This function is to select and deduct the sites used to
# fit or valid the model
# Input: 	data
#		s = the site numbers to be selected/deselected, e.g., c(4,7,10)

    if (missing(var.name)) {
        stop("Error: need to define the var.name")
    }
    if(!is.data.frame(data)){
        stop("Error: data should be in data.frame")
    }
    if(!var.name %in%names(data)){
        stop("Error: var.name should be in data")
    }
    data<-data[,c(var.name,names(data)[!names(data)%in%var.name])]
    s.index <- unique(data[, var.name[1]])
    if (reverse == FALSE) {
        data <- data[data[,1] %in% (s.index[s.index %in% s]),]
        data
    }
    else {
        data <- data[data[,1] %in% (s.index[!s.index %in% s]),]
        data
    }
}
# spT.subset<-function(data, var.name, s = NULL, reverse=FALSE) 
#{
#
# This function is to select and deduct the sites used to
# fit or valid the model
# Input: 	data
#		s = the site numbers to be selected/deselected, e.g., c(4,7,10)
#	if(missing(var.name)){
#		stop("Error: need to define the var.name")
#	}
#   s.index <- unique(data[,var.name[1]])
#    if(reverse==FALSE){
#		dat <- subset(data, s.index %in% s)
#		dat
#	}
#	else{
#		dat <- subset(data, !(s.index %in% s))
#		dat
#	}
#}
##
# spT.data.selection<-function(data, rs = NULL, s = NULL, reverse=FALSE) 
#{
#
# This function is to select and deduct the sites used to
# fit or valid the model
# Input: 	data
#		s = the site numbers to be selected/deselected, e.g., c(4,7,10)
#		rs = the total number of random sites to be selected , e.g., 3
#  if(reverse==FALSE){
#	if(is.null(rs) & !is.null(s)){
	# not random
#		rs<-NULL
#		dat<-NULL
#			for(i in 1:length(s)){
#			dat<-rbind(dat,data[data[,var.name[1]]==s[i], ])
#			}	
#		dat	
#	}	
#	else if (!is.null(rs) & is.null(s)){
	# for randomly selected sites
#		s<-NULL
#		a.s<-unique(data[,1])
#		a.s<-sort(sample(a.s, rs))
#		cat('Randomly selected sites =', a.s, '\n')
#		dat<-NULL
#			for(i in 1:rs){
#			dat<-rbind(dat,data[data[,1]==a.s[i], ])
#			}	
#		dat
#	}
#  }
#  else{
#	if(is.null(rs) & !is.null(s)){
#		rs<-NULL
#			for(i in 1:length(s)){
#			data<-data[data[,1] != s[i], ]
#			}	
#		data	
#	}	
#	else{
#		s<-NULL
#		a.s<-unique(data[,1])
#		a.s<-sort(sample(a.s, rs))
#		cat('Randomly deducted sites =', a.s, '\n')
#			for(i in 1:rs){
#			data<-data[data[,1] != s[i], ]
#			}	
#		data
#	}
#  }
# }
##
## grid coordinates
##
 spT.grid.coords<-function(Longitude=c(max,min),Latitude=c(max,min),by=c(NA,NA))
{
      max.lon <- Longitude[1]
      min.lon <- Longitude[2]
      max.lat <- Latitude[1]
      min.lat <- Latitude[2]
      if(is.na(by[1]) || is.na(by[2])){
        stop('Error: need to specify grid dimension, n x m')
      }
      knots.lon <- seq(max.lon,min.lon,length=by[1])
      knots.lat <- seq(max.lat,min.lat,length=by[2])
      #knots.coords <- cbind(rep(knots.lon,by),c(outer(rep(1,by),knots.lat)))
      knots.coords <- cbind(c(outer(rep(1,by[2]),knots.lon)),rep(knots.lat,by[1]))
      as.matrix(knots.coords)
}
##
## Percentage Coverage
##
 spT.pCOVER<-function(z=NULL,zup=NULL,zlow=NULL,zsample=NULL,level=95)
{
    if(is.null(z)){
      stop('\n # Provide observed values\n')
    }
    if(!is.null(zsample)){
      z<-c(z)
      low<-(1-level/100)/2
      up<-1-low
      x<-apply(zsample,1,quantile,prob=c(low,up))
      x<-rbind(x,z)
      x<-t(x)
      x<-na.exclude(x)
      y<-x[x[, 1] <= x[, 3] & x[, 3] <= x[, 2], ]
      round(length(y[,1])/length(x[,1])*100,2)
    }
    else if(is.null(zsample)){
      z<-as.matrix(z)
      zlow<-as.matrix(zlow)
      zup<-as.matrix(zup)
      x<-cbind(zlow,z,zup)
      x<-na.exclude(x)
      y<- x[x[,1] <= x[,2] & x[,2] <= x[,3],]
      round(length(y[,1])/length(x[,1])*100,2)
    }
}
##
## Geodetic distance in K.M. or Miles
##
 spT.geodist<-function(Lon, Lat, KM=TRUE){
 #
 # This function is for geodistance using C/C++
 #
   n<-length(Lat)
   if(KM==TRUE){	
   ds<-.C("GeoDist_km", as.integer(n),as.double(Lat),
    as.double(Lon),out=matrix(double(n*n),n,n))$out
   }
   else {
   ds<-.C("GeoDist_miles", as.integer(n),as.double(Lat),
    as.double(Lon),out=matrix(double(n*n),n,n))$out
   }
   #dis<-matrix(ds,length(Lat),length(Lat))
   #return(dis)
   ds
 }
##
## Geodetic distance: another approach: two points
##
 spT.geo.dist <- function(point1, point2){
 #
 # The following program computes the distance on the surface 
 # of the earth between two points point1 and point2. 
 # Both the points are of the form (Longitude, Latitude)
 #
   R <- 6371
   p1rad <- point1 * pi/180
   p2rad <- point2 * pi/180
   d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))	
   d <- acos(d)
   R*d
 }
##
## Geodetic distance: using points 1:(lon, lat) 2:(lon, lat)
##
spT.geo_dist <- function(points)
{
        point1 <- c(points[1], points[2])
        point2 <- c(points[3], points[4])
        #The following program computes the distance on the surface 
        #The argument points is of the form (Long, Lat, Long, Lat)
        dist <- 0
        if(sum((point1 - point2)^2) > 1e-05) {
                R <- 6371
                p1rad <- (point1 * pi)/180
                p2rad <- (point2 * pi)/180
                d <- sin(p1rad[2]) * sin(p2rad[2]) + cos(p1rad[2]) * cos(p2rad[
                        2]) * cos(abs(p1rad[1] - p2rad[1]))
                d <- acos(d)
                dist <- R * d
        }
        dist
}
##
## To keep the long lat positions based on tol.dist
##
 spT.keep.morethan.dist <- function(coords, tol.dist=100)  
{
 #
 # coords must have two columns named long and lat 
 #
   a <- as.data.frame(coords)
   names(a) <- c("long","lat")
   n <- nrow(coords)
   c1 <- rep(1:n, each=n)
   c2 <- rep(1:n, n)
   b <- matrix(NA, nrow=n, ncol=n)
   w <- as.vector(upper.tri(b))
   bigmat <- matrix(0, nrow=n*n, ncol=7)
   bigmat[, 1] <- c1
   bigmat[, 2] <- c2
   bigmat[, 3] <- a$long[c1]
   bigmat[, 4] <- a$lat[c1]
   bigmat[, 5] <- a$long[c2]
   bigmat[, 6] <- a$lat[c2]
   ubmat <- bigmat[w, ]
   ubmat[,7] <- as.vector(apply(ubmat[,3:6], 1, spT.geo_dist)) 

   v <- ubmat[,7]
   w <- ubmat[v<tol.dist, ]

   z <- unique(w[,1])
   a <- coords[-z,  ]
   a
}
##
## This function extracts the data using formula
##
 Formula.matrix <- function(formula, data, na.rm=FALSE)
{
      if(na.rm == FALSE){
           mt <- terms.formula(formula, data=data, specials=c("sp","tp"))
           if(missing(data)) data <- sys.frame(sys.parent())
           spcheck<-attr(mt, "specials")$sp
           tpcheck<-attr(mt, "specials")$tp
           mf <- model.frame(mt, data=data, na.action=na.pass)
      }
      else{
           mt <- terms.formula(formula, data=data, specials=c("sp","tp"))
           if(missing(data)) data <- sys.frame(sys.parent())
           spcheck<-attr(mt, "specials")$sp
           tpcheck<-attr(mt, "specials")$tp
           mf <- model.frame(mt, data=data)
      }
      # null model support
           "%w/o%" <- function(x, y) x[!x %in% y]
           x.names <- attr(mt, "term.labels") # X variable names
           X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
           if((!is.null(spcheck)) & (is.null(tpcheck))){
           # only spatially varying  
                x.names.sp <- x.names[attr(mt, "specials")$sp-1]
                X.sp <- X[,c(x.names.sp)]
                x.names <- dimnames(X)[[2]]
                x.names <- x.names %w/o% x.names.sp
                X <- X[,c(x.names)]
                X <- as.matrix(X)
                dimnames(X)[[2]]<-as.list(x.names) 
                X.sp <- as.matrix(X.sp) 
                dimnames(X.sp)[[2]]<-as.list(x.names.sp) 
                X.tp <- NULL; x.names.tp<-NULL
                if(dim(X)[[2]]==0){
                  X<-as.matrix(rep(0, dim(mf)[[1]])); 
                  x.names <- "(Intercept)"
                }  
           }
           else if((!is.null(tpcheck)) & (is.null(spcheck))) {
           # only temporally varying
                x.names.tp <- x.names[attr(mt, "specials")$tp-1]
                X.tp <- X[,c(x.names.tp)]
                x.names <- dimnames(X)[[2]]
                x.names <- x.names %w/o% x.names.tp
                X <- X[,c(x.names)]
                X <- as.matrix(X)
                dimnames(X)[[2]]<-as.list(x.names) 
                X.tp <- as.matrix(X.tp) 
                dimnames(X.tp)[[2]]<-as.list(x.names.tp) 
                X.sp <- NULL; x.names.sp<-NULL 
                if(dim(X)[[2]]==0){
                  X<-as.matrix(rep(0, dim(mf)[[1]])); 
                  x.names <- "(Intercept)"
                }  
           }
           else if((!is.null(spcheck)) & (!is.null(tpcheck))){
           # both spatially and temporally varying
                x.names.sp <- x.names[attr(mt, "specials")$sp-1]
                x.names.tp <- x.names[attr(mt, "specials")$tp-1]
                X.sp <- X[,c(x.names.sp)]
                X.tp <- X[,c(x.names.tp)]
                x.names <- dimnames(X)[[2]]
                x.names <- x.names %w/o% x.names.sp
                x.names <- x.names %w/o% x.names.tp
                X <- X[,c(x.names)]
                X <- as.matrix(X)
                dimnames(X)[[2]]<-as.list(x.names) 
                X.sp <- as.matrix(X.sp) 
                dimnames(X.sp)[[2]]<-as.list(x.names.sp) 
                X.tp <- as.matrix(X.tp) 
                dimnames(X.tp)[[2]]<-as.list(x.names.tp) 
                if(dim(X)[[2]]==0){
                  X<-as.matrix(rep(0, dim(mf)[[1]])); 
                  x.names <- "(Intercept)"
                }  
          }
           else {
                x.names <- dimnames(X)[[2]] 
                X.sp <- NULL; x.names.sp<-NULL; X.tp <- NULL; x.names.tp<-NULL 
           }
           Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
       #
       return(list(Y=Y, X=X, x.names=x.names, X.sp=X.sp, x.names.sp=x.names.sp, X.tp=X.tp, x.names.tp=x.names.tp))
}
##
## This function extracts the coords using formula
##
 Formula.coords <- function(formula, data, na.rm=TRUE)
{
     mt <- terms(formula, data=data)
     if(missing(data)) data <- sys.frame(sys.parent())
     mf <- model.frame(mt, data=data)
	 # null model support
     X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
     X <- as.matrix(X[,2:3])         # X matrix
	 X <- unique(X); dimnames(X)<-NULL
     X
}
##
## Summary Statistics
##
 spT.Summary.Stat <- function(y) 
{

     ## define Upper and lower limit function
       up.low.limit<-function(y,limit)
       {
         #
          y<-sort(y)
          y<-y[limit]
          y
         #
       }
       #
         if(is.vector(y)==TRUE){
            y <- matrix(y[!is.na(y)])
         }
         else{
            y <- t(y)
         } 
         N <- length(y[1,  ])
         nItr <- length(y[, 1])
         z <- matrix(nrow = N, ncol = 5)
         dimnames(z) <- list(dimnames(y)[[2]], c("Mean","Median","SD","Low2.5p","Up97.5p"))
        #
         if (nItr < 40) {
           stop("\n##\n# Error: number of samples must be >= 40\n##")
         }
         z[, 1] <- apply(y,2,mean)
         z[, 2] <- apply(y,2,median)
         z[, 3] <- apply(y,2,sd)
         #z[, 4] <- apply(y,2,quantile,0.025)
         #z[, 5] <- apply(y,2,quantile,0.975)
        nl <- as.integer(nItr * 0.025)
        nu <- as.integer(nItr * 0.975)
         z[, 4] <- apply(y,2,up.low.limit,limit=nl)
         z[, 5] <- apply(y,2,up.low.limit,limit=nu)
         
        as.data.frame(z)
}
##
## Predictive Model Choice Criteria
##
 PMCC<-function(z=NULL, z.mean=NULL, z.sd=NULL, z.samples=NULL)
{
    #
    # Predictive Model Choice Criteria
    #
    if(is.null(z)){
      stop("Error: need to provide z values.")
    }
    #
    if(is.null(z.mean) | is.null(z.sd)){
      if(!is.null(z.mean)){
       stop("Error: need to provide z.sd value.")
      }   
      if(!is.null(z.sd)){
       stop("Error: need to provide z.mean value.")
      }   
      if(is.null(z.samples)){
       stop("Error: need to provide z.samples value.")
      }
    }
    #
    if(!is.null(z.samples)){
     if ( !is.matrix(z.samples)) {
         stop("Error: z.samples must be a (N x nItr) matrix")
     }
     if (dim(z.samples)[1] != length(z)) {
         stop("Error: observations in z.samples in each iteration must be equal to length of z")
     }
     if ( dim(z.samples)[2] < 40) {
         stop("Error: samples are too small to obtain summary statistics")
     }
     sum.stat = matrix(NA,length(c(z)),6)
     sum.stat[,1:5] = as.matrix(spT.Summary.Stat(z.samples))
     sum.stat[,6] = c(z)
     sum.stat = sum.stat[!is.na(sum.stat[,6]),]
     goodness.of.fit = round(sum((sum.stat[,1]-sum.stat[,6])^2),2)
     penalty = round(sum(sum.stat[,3]^2),2)
     pmcc =  round(goodness.of.fit + penalty,2)
     out = NULL
     out$pmcc = pmcc; 
     out$goodness.of.fit = goodness.of.fit
     out$penalty = penalty
     #class(out) <- "PMCC"
     out
    }
    else{
     if(is.null(z.mean) | is.null(z.sd)){
       stop("Error: need to provide z.mean and/or z.sd values.")
     }
     if(length(c(z)) != length(c(z.mean))){
       stop("Error: z and z.mean should be in same length.")
     }
     if(length(c(z)) != length(c(z.sd))){
       stop("Error: z and z.sd should be in same length.")
     }
     sum.stat = matrix(NA,length(c(z)),3)
     sum.stat[,1] = c(z)
     sum.stat[,2] = c(z.mean)
     sum.stat[,3] = c(z.sd)
     sum.stat = sum.stat[!is.na(sum.stat[,1]),]
     goodness.of.fit = round(sum((sum.stat[,1]-sum.stat[,2])^2),2)
     penalty = round(sum(sum.stat[,3]^2),2)
     pmcc =  round(goodness.of.fit + penalty,2)
     out = NULL
     out$pmcc = pmcc; 
     out$goodness.of.fit = goodness.of.fit
     out$penalty = penalty
     #class(out) <- "PMCC"
     out
    }
}

##
## To identify the Integer Number
##

#is.wholenumber<-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

##
## To check sites that has <tol km of distances
## 
 spT.check.locations<-function(fit.locations, pred.locations,
              method="geodetic:km", tol=5){
  #
      #
      if(!method %in% c("geodetic:km", "geodetic:mile", "euclidean",
        "maximum", "manhattan", "canberra")){
        stop("\n# Error: correctly define distance.method \n")
      }
      #
           coords.all <- rbind(fit.locations,pred.locations)
           tn.fitsites <- length(fit.locations[, 1])
           nfit.sites <- 1:tn.fitsites
           tn.predsites <- length(coords.all[, 1]) - tn.fitsites
           npred.sites <- (tn.fitsites + 1):(length(coords.all[, 1]))
      #
      if(method=="geodetic:km"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else{
       coords.D <- as.matrix(dist(coords.all, method, diag = T, upper = T))
      }  
           coords.D[is.na(coords.D)]<-0
      #
           diag(coords.D)<-NA
      # 
           fdmis<-coords.D[nfit.sites, npred.sites]
     if(is.matrix(fdmis)==TRUE){
           fdmis<-cbind(c(t(fdmis)),1:dim(fdmis)[[2]],sort(rep(1:dim(fdmis)[[1]],dim(fdmis)[[2]]))) # add pred sites and fitted sites
           fdmis<-fdmis[fdmis[,1] < tol,]
      #
           if(!is.na(fdmis[1])==TRUE){
            cat("#\n# Tolerance Limit:", paste(tol))
            cat("\n# There are some Prediction locations very close to the Fitted locations.\n#\n")
            fdmis<-matrix(fdmis,(length(fdmis)/3),3) 
            for(i in 1:dim(fdmis)[[1]]){
            print(paste("Distance:", round(fdmis[i,1],2)," Predicted location:",fdmis[i,2]," Fitted location:", fdmis[i,3],""))
            }
            cat("#\n# Romove the locations and run again. \n#\n")
            dimnames(fdmis)[[2]]<-c('distance','pred_location','fit_location')
            fdmis
           }
           else{
            cat("#\n# Tolerance Limit (unit):", paste(tol))
            #cat("\n# Fitted and Predicted location distances are alright \n#\n")  
            cat("\n# Location distances are alright \n#\n")  
           }
      }
      else{
           fdmis<-cbind(c(fdmis),1:length(fdmis)) # 
           fdmis<-fdmis[fdmis[,1] < tol,]
           if(!is.na(fdmis[1])==TRUE){
            cat("#\n# Tolerance Limit:", paste(tol))
            cat("\n# There are some Prediction locations very close to the Fitted locations.\n#\n")
            fdmis<-matrix(fdmis) 
            for(i in 1:dim(fdmis)[[1]]){
            print(paste("Distance:", round(fdmis[i,1],2)," Predicted location:",1," Fitted location:", fdmis[i,2],""))
            }
            cat("#\n# Romove the locations and run again. \n#\n")
            dimnames(fdmis)[[2]]<-c('distance','fit_location')
            fdmis
           }
           else{
            cat("#\n# Tolerance Limit (unit):", paste(tol))
            #cat("\n# Fitted and Predicted location distances are alright \n#\n")  
            cat("\n# Location distances are alright \n#\n")  
           }
       }
}
##
## To check sites inside codes

## 
spT.check.sites.inside<-function(coords, method, tol=0.1){
      #
      #
      if(!method %in% c("geodetic:km", "geodetic:mile", "euclidean",
        "maximum", "manhattan", "canberra")){
        stop("\n# Error: correctly define distance.method \n")
      }
      #
      #
      if(method=="geodetic:km"){
          fdm<- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
          fdm<- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=FALSE))
      }
      else{
           fdm<- as.matrix(dist(coords, method, diag = TRUE, upper = TRUE))
      } 
      #
           diag(fdm)<-NA
      # 
           fdm<-cbind(c(fdm),1:dim(fdm)[[2]],sort(rep(1:dim(fdm)[[1]],dim(fdm)[[2]]))) 
      #
           fdm<-fdm[!is.na(fdm[,1]),]
      #
           #tol <- 0.01
           fdmis<-fdm[fdm[,1] < tol,]
      #
           if(!is.na(fdmis[1])==TRUE){
            cat("#\n# There are some locations very close:")
            print(paste("( < ",tol," unit) to each other."))
            fdmis<-matrix(fdmis,(length(fdmis)/3),3) 
            for(i in 1:dim(fdmis)[[1]]){
            print(paste("Distance (unit):", round(fdmis[i,1],2)," site:",fdmis[i,2]," site:", fdmis[i,3],""))
            }
            dimnames(fdmis)[[2]]<-c('dist_km','pred_site','fit_site')
            fdmis
            cat("#\n# Romove the sites and run again. \n#\n")
            stop("Error: Termination.")
           }
       #
           #tol <- 0.5
           #fdmis<-fdm[fdm[,1] < tol,]
      #
           #if(!is.na(fdmis[1])==TRUE){
           # cat("#\n# Warnings: There are some locations very close ( < 0.5 unit) to each other.\n#\n")
           #}
      #
}
##
##
#check.sites.inside.pred<-spT.check.locations
##
## convert seconds into min. hour. and day
##
 fnc.time<-function(t)
{
     #
     if(t < 60){
      t <- round(t,2)
      tt <- paste(t," - Sec.")
      cat(paste("##\n# Elapsed time:",t,"Sec.\n##\n"))
     } 
     #
     if(t < (60*60) && t >= 60){
      t1 <- as.integer(t/60)
      t <- round(t-t1*60,2) 
      tt <- paste(t1," - Mins.",t," - Sec.")
      cat(paste("##\n# Elapsed time:",t1,"Min.",t,"Sec.\n##\n"))
     }
     #
     if(t < (60*60*24) && t >= (60*60)){
      t2 <- as.integer(t/(60*60))
      t <- t-t2*60*60
      t1 <- as.integer(t/60)
      t <- round(t-t1*60,2) 
      tt <- paste(t2," - Hour/s.",t1," - Mins.",t," - Sec.")
      cat(paste("##\n# Elapsed time:",t2,"Hour/s.",t1,"Min.",t,"Sec.\n##\n"))
     }
     #
     if(t >= (60*60*24)){
      t3 <- as.integer(t/(60*60*24))
      t <- t-t3*60*60*24
      t2 <- as.integer(t/(60*60))
      t <- t-t2*60*60
      t1 <- as.integer(t/60)
      t <- round(t-t1*60,2)
      tt <- paste(t3," - Day/s.",t2," - Hour/s.",t1," - Mins.",t," - Sec.")
      cat(paste("##\n# Elapsed time:",t3,"Day/s.",t2,"Hour/s.",t1,"Mins.",t,"Sec.\n##\n"))
     }
     #
     tt
}
##
## validation criteria
##
spT.validation <- function(z, zhat)
{
 ##
 ## Validation Mean Squared Error (VMSE) 
 ##
 VMSE <- function(z, zhat)
 {
       z<-as.matrix(z)
       zhat<-as.matrix(zhat)
       x <- c(z-zhat)
       u <- x[!is.na(x)]
       round(sum(u^2)/length(u), 4)
 }
 ##
 ## Root Mean Squared Error (RMSE) 
 ##
 RMSE<-function (z,zhat) 
 {
       z<-as.matrix(z)
       zhat<-as.matrix(zhat)
       x <- c(z-zhat)
       u <- x[!is.na(x)]
       round(sqrt(sum(u^2)/length(u)), 4)
 }
 ##
 ## Mean Absolute Error (MAE) 
 ##
 MAE<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- abs(c(zhat-z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
 }
 ##
 ## Mean Absolute Percentage Error (MAPE) 
 ##
 MAPE<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- abs(c(zhat-z))/z
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u)*100, 4)
 }
 ##
 ## Bias (BIAS) 
 ##
 BIAS<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- c(zhat-z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
 }
 ##
 ## Relative Bias (rBIAS) 
 ##
 rBIAS<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- c(zhat-z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u)*mean(z,na.rm=TRUE)), 4)
 }
 ##
 ## Relative Mean Separation (rMSEP) 
 ##
 rMSEP<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- c(zhat-z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat,na.rm=TRUE)-z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
 }
 ##
 cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Relative Mean Separation (rMSEP)\n##\n") 
 ##
   out<-NULL
   out$MSE<-VMSE(c(z), c(zhat))
   out$RMSE<-RMSE(c(z), c(zhat))
   out$MAE<-MAE(c(z), c(zhat))
   out$MAPE<-MAPE(c(z), c(zhat))
   out$BIAS<-BIAS(c(z), c(zhat))
   out$rBIAS<-rBIAS(c(z), c(zhat))
   out$rMSEP<-rMSEP(c(z), c(zhat))
   unlist(out)
}
##
## MCMC plots for individual values (trace, density, acf)
##
MCMC.plot.obj<-function(post_val, nItr, nBurn=0, 
                name=c('....'), ACF="FALSE", PARTIAL.acf="FALSE") 
  {
      nBurn=nBurn+1; x<-post_val;
      #if(missing(initial_val)){ 
      #   y<-0
      #}
      #else{
      #   y<-initial_val
      #} 
      if(missing(nItr)){ 
         nItr<-length(post_val)
      }
      if(ACF==FALSE & PARTIAL.acf==FALSE){
        #windows()
        #dev.new()
        par(mfrow=c(1,2))
        #plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
        #   ylab = paste("Values of  (", name, ")", sep=''), type = "l",
        #   ylim=c(min(y,x[nBurn:nItr]),max(y,x[nBurn:nItr])),main='Trace plot')
        #abline(h = y, lty = 2,col="blue")
        #plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(y,x[nBurn:nItr]),max(y,x[nBurn:nItr])))
        #abline(v = y, lty = 2,col="blue")
        plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
           ylab = paste("Values of  (", name, ")", sep=''), type = "l",
           ylim=c(min(x[nBurn:nItr]),max(x[nBurn:nItr])),main='Trace plot')
        plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(x[nBurn:nItr]),max(x[nBurn:nItr])))
     }	
     else if(ACF==TRUE & PARTIAL.acf==TRUE){
      #windows()
      #dev.new() 
      par(mfrow=c(1,4))
      plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
        ylab = paste("Values of  (", name, ")", sep=''), type = "l", 
        ylim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])),main='Trace plot')
      #abline(h = y, lty = 2,col="blue")
      plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])))
      #abline(v = y, lty = 2,col="blue")
      acf(x[nBurn:nItr],main='ACF plot',col="red")
      pacf(x[nBurn:nItr],main='Partial ACF plot',col="red")
     }
     else if(ACF=="TRUE" & PARTIAL.acf=="FALSE"){
      #dev.new() 
      par(mfrow=c(1,3))
      plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
        ylab = paste("Values of  (", name, ")", sep=''), type = "l", 
        ylim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])),main='Trace plot')
      #abline(h = y, lty = 2, col="blue")
      plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])))
      #abline(v = y, lty = 2, col="blue")
      acf(x[nBurn:nItr],main='ACF plot',col="red")
     }
     else {

     }
  }
##
## Multiple imputation for initial values using Amelia with m=1
##
#imputation.z<-function(x){
  #
  # x is the r*T by n matrix
  #
#  x<-amelia(x,m=1)$imputation[[1]]
#  x
#}
##
## To use in coda package
##
as.mcmc.spT<-function(x, ...){

    if (x$combined.fit.pred == TRUE) {
        stop("\n# Error: not available for combined fit-prediction option. Please read the MCMC samples from the *.txt file manually.")
    }
    model <- x$model
    if (is.null(model) == TRUE) {
        stop("\n# Error: need to define the model")
    }
    else if (model == "AR") {
        r <- x$r
        p <- x$p
      #          
      if(x$cov.fnc=="matern"){
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip), t(x$nup))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi", "nu")
      }
      else {
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi")
      }
      #
        para<-t(para)
        para<-mcmc(para)
        para
    }
    else if (model == "GPP") {
        r <- x$r
        p <- x$p
      #          
      if(x$cov.fnc=="matern"){
        if((!is.null(x$sp.covariate.names))){  
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$sig2betap), t(x$phip), t(x$nup))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "sig2beta", "phi", "nu")
        }
        else{
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip), t(x$nup))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi", "nu")
        }
      }
      else {
        if((!is.null(x$sp.covariate.names))){  
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$sig2betap), t(x$phip))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "sig2beta", "phi")
        }
        else{
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi")
        }
      }
      #
        para<-t(para)
        para<-mcmc(para)
        para
    }
    else if (model == "GP") {
        r <- x$r
        p <- x$p
      #          
      if(x$cov.fnc=="matern"){
         if((is.null(x$sp.covariate.names)) & (is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$phip), t(x$nup))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "phi", "nu")
         }
         else if((!is.null(x$sp.covariate.names)) & (is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$sig2betap), t(x$phip), t(x$nup))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "sig2beta", "phi", "nu")
         }
         else if((is.null(x$sp.covariate.names)) & (!is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$sig2deltap), t(x$sig2op), t(x$phip), t(x$nup))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "sig2deltap", "sig2op", "phi", "nu")
         } 
         else if((!is.null(x$sp.covariate.names)) & (!is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$sig2betap), t(x$sig2deltap), t(x$sig2op), t(x$phip), t(x$nup))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "sig2beta", "sig2deltap", "sig2op", "phi", "nu")
         }
         else {
            stop("Error")
         } 
      }
      else {
         if((is.null(x$sp.covariate.names)) & (is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$phip))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "phi")
         }
         else if((!is.null(x$sp.covariate.names)) & (is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$sig2betap), t(x$phip))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "sig2beta", "phi")
         }
         else if((is.null(x$sp.covariate.names)) & (!is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$sig2deltap), t(x$sig2op), t(x$phip))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "sig2deltap", "sig2op", "phi")
         } 
         else if((!is.null(x$sp.covariate.names)) & (!is.null(x$tp.covariate.names))){  
           para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$sig2betap), t(x$sig2deltap), t(x$sig2op), t(x$phip))
           dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "sig2beta", "sig2deltap", "sig2op", "phi")
         }
         else {
            stop("Error")
         } 
      }
      #
        para<-t(para)
        para<-mcmc(para)
        para
    }
    else {
    }
}
##
## use of model.frame
##
model.frame.spT<-function(formula, ...){
#   tmp<-cbind(formula$Y,formula$X[,-1])
#   dimnames(tmp)[[2]]<-dimnames(attr(terms(formula$call),"factors"))[[1]]
#   tmp
   if(formula$combined.fit.pred==TRUE){
      stop("\n# Error: not useful for output with combined fit and predict \n")
    }
   else{
      model.frame(formula$call,formula$data,na.action=na.pass)
   }
}
##
## use of model.matrix
##
model.matrix.spT<-function(object, ...){
   if(object$combined.fit.pred==TRUE){
      stop("\n# Error: not useful for output with combined fit and predict \n")
    }
   else{
      Formula.matrix(object$call,object$data)[[2]]
   }
}
##
## to include dimnames: function used in function summary statistics for spatially varying coef
##
sp.dimname.fnc<-function(x,names,n,q){
  dimnames(x)[[2]][1:(n*q)]<-1:(n*q)
  for(i in 1:q){
    dimnames(x)[[2]][(1+(i-1)*n):(n*i)]<-rep(paste(names[i],"site",1:n,sep=""))
  }
  x
}
##
## to include dimnames: function used in function summary statistics for temporally varying coef
##
tp.dimname.fnc<-function(x,names,u,T){
  dimnames(x)[[2]][1:(u*T)]<-1:(u*T)
  for(i in 1:u){
    dimnames(x)[[2]][(1+(i-1)*T):(T*i)]<-rep(paste(names[i],"time",1:T,sep=""))
  }
  x
}
##
## use of summary
##
summary.spT<-function(object, digits=4, package="spTimer", ...){
   coefficient=NULL
   if(package=="coda"){
    if(object$combined.fit.pred==TRUE){
      stop("\n# Error: coda package is not useful for output with combined fit and predict \n")
    }
    else{
     if(is.null(coefficient)){
       cat("\n#### MCMC summary statistics using coda package ####\n")
       if(!is.null(object$sp.covariate.names)) {cat("\n## \n# Spatially varying parameters are not included.\n##\n ")}
       tmp<-as.mcmc(object)
       summary(tmp, ...)
     }
     else if(coefficient=="spatial"){
       if(is.null(object$sp.covariate.names)) {stop("\n## \n# Spatially varying parameters are not used in the model.\n##\n ")}
       cat("\n#### MCMC summary statistics for the spatially varying coefficients ####\n")
       tmp<-as.mcmc(t(object$betasp))
       if(object$model=="GPP"){n<-object$knots} 
       else{n<-object$n}
       tmp<-sp.dimname.fnc(x=tmp,names=object$sp.covariate.names,n=n,q=object$q) 
       summary(tmp, ...)
     }
     else if(coefficient=="temporal"){
       if(is.null(object$tp.covariate.names)) {stop("\n## \n# Temporally varying dynamic parameters are not used in the model.\n##\n ")}
       cat("\n#### MCMC summary statistics for the temporally varying dynamic coefficients ####\n")
       tmp<-as.mcmc(t(object$betatp))
       tmp<-tp.dimname.fnc(x=tmp,names=object$tp.covariate.names,u=object$u,T=object$T) 
       summary(tmp, ...)
     }
     else if(coefficient=="rho"){
       if(is.null(object$rhotp)) {stop("\n## \n# rho parameter has not been sampled in the temporally varying dynamic model.\n##\n ")}
       cat("\n#### MCMC summary statistics for the rho coefficients ####\n")
       tmp<-as.mcmc(t(object$rhotp))
       dimnames(tmp)[[2]]<-paste("rho",1:object$u,sep="") 
       summary(tmp, ...)
     }
     else{
       stop("Error: the argument coefficient only takes charecter 'spatial' and 'temporal'.")
     }
    }
   }
   else{
     if(is.null(coefficient)){
       if((!is.null(object$sp.covariate.names)) & (!is.null(object$tp.covariate.names))){cat("\n## \n# Spatially and temporally varying parameters are not included.\n##\n ")}
       else if((!is.null(object$sp.covariate.names)) & (is.null(object$tp.covariate.names))) {cat("\n## \n# Spatially varying parameters are not included.\n##\n ")}
       else if((is.null(object$sp.covariate.names)) & (!is.null(object$tp.covariate.names))) {cat("\n## \n# Temporally varying dynamic parameters are not included.\n##\n ")}
       else { }
	   #ans<-NULL
	   #ans$Model=object$model;ans$Call=object$call;ans$Iterations=object$iterations;ans$nBurn=object$nBurn;ans$PMCC=object$pmcc;ans$Parameters=round(object$parameter,digits=digits)
       #ans
       print(object)
       cat("-----------------------------------------------------"); cat('\n');
       cat("Parameters:\n")
       print(round(object$parameter,digits=digits)); #cat("\n");
       cat("-----------------------------------------------------"); cat('\n');
     }
     else if(coefficient=="spatial"){
       if(is.null(object$sp.covariate.names)) {stop("\n## \n# Spatially varying parameters are not used in the model.\n##\n ")}
       cat("\n#### MCMC summary statistics for the spatially varying coefficients ####\n")
       tmp<-as.mcmc(t(object$betasp))
       if(object$model=="GPP"){n<-object$knots} 
       else{n<-object$n}
       tmp<-sp.dimname.fnc(x=tmp,names=object$sp.covariate.names,n=n,q=object$q) 
       summary(tmp, ...)
     }
     else if(coefficient=="temporal"){
       if(is.null(object$tp.covariate.names)) {stop("\n## \n# Temporally varying dynamic parameters are not used in the model.\n##\n ")}
       cat("\n#### MCMC summary statistics for the spatially varying coefficients ####\n")
       tmp<-as.mcmc(t(object$betatp))
       tmp<-tp.dimname.fnc(x=tmp,names=object$tp.covariate.names,u=object$u,T=object$T) 
       summary(tmp, ...)
     }
     else if(coefficient=="rho"){
       if(is.null(object$rhotp)) {stop("\n## \n# rho parameter has not been sampled in the temporally varying dynamic model.\n##\n ")}
       cat("\n#### MCMC summary statistics for the rho coefficients ####\n")
       tmp<-as.mcmc(t(object$rhotp))
       dimnames(tmp)[[2]]<-paste("rho",1:object$u,sep="") 
       summary(tmp, ...)
     }
     else{
       stop("Error: the argument coefficient only takes charecter 'spatial' and 'temporal'.")
     }
   }
}
##
## use of plot
##
plot.spT<-function(x, residuals=FALSE, ...){
   coefficient=NULL 
   if(as.logical(residuals)==FALSE){
     if(x$combined.fit.pred==TRUE){
       if(!is.null(x$sp.covariate.names)) {cat("## \n# Spatially varying parameters are not included.\n##\n ")}
       cat("\n## Only residual plots are available for output with combined fit and predict option \n\n")
       plot(x$fitted[,1],residuals(x),ylab="Residuals",xlab="Fitted values");abline(h=0,lty=2);title("Residuals vs Fitted")
       par(ask=TRUE)
       qqnorm(residuals(x));qqline(residuals(x),lty=2)
     }
     else{
       if(is.null(coefficient)){
         if((!is.null(x$sp.covariate.names)) & (!is.null(x$tp.covariate.names))){cat("\n## \n# Spatially and temporally varying parameters are not included.\n##\n ")}
         else if((!is.null(x$sp.covariate.names)) & (is.null(x$tp.covariate.names))) {cat("\n## \n# Spatially varying parameters are not included.\n##\n ")}
         else if((is.null(x$sp.covariate.names)) & (!is.null(x$tp.covariate.names))) {cat("\n## \n# Temporally varying dynamic parameters are not included.\n##\n ")}
         else {  }
         tmp<-as.mcmc(x)
         plot(tmp, ...)
       }
       else if(coefficient=="spatial"){
         if(is.null(x$sp.covariate.names)) {stop("\n## \n# Spatially varying parameters are not used in the model.\n##\n ")}
         tmp<-as.mcmc(t(x$betasp))
         if(x$model=="GPP"){n<-x$knots} 
         else{n<-x$n}
         tmp<-sp.dimname.fnc(x=tmp,names=x$sp.covariate.names,n=n,q=x$q) 
         plot(tmp, ...)
       }
       else if(coefficient=="temporal"){
         if(is.null(x$tp.covariate.names)) {stop("\n## \n# Temporally varying dynamic parameters are not used in the model.\n##\n ")}
         tmp<-as.mcmc(t(x$betatp))
         tmp<-tp.dimname.fnc(x=tmp,names=x$tp.covariate.names,u=x$u,T=x$T) 
         plot(tmp, ...)
       }
       else if(coefficient=="rho"){
         if(is.null(x$rhotp)) {stop("\n## \n# rho parameter has not been sampled in the temporally varying dynamic model.\n##\n ")}
         tmp<-as.mcmc(t(x$rhotp))
         dimnames(tmp)[[2]]<-paste("rho",1:x$u,sep="") 
         plot(tmp, ...)
       }
       else{
         stop("Error: the argument coefficient only takes charecter 'spatial' and 'temporal'.")
       }
     } 
   }
   else{
     if(!is.null(x$sp.covariate.names)) {cat("## \n# Models fitted with spatially varying parameters.\n## \n ")}
     plot(x$fitted[,1],residuals(x),ylab="Residuals",xlab="Fitted values");abline(h=0,lty=2);title("Residuals vs Fitted")
     par(ask=TRUE)
     qqnorm(residuals(x));qqline(residuals(x),lty=2)
   } # end of 2nd if 
}
##
## use of contour
##
#contour.spT<-function(x, surface="Mean", time=c(1), ...){
#	z<-array(fitted(x)[,paste(surface)],dim=c(x$T*x$r,x$n))
#	xyz<-cbind(x$coords,c(z[time,]))
#	xyz<-interp(x=xyz[,1], y=xyz[,2], z=xyz[,3], xo=seq(min(xyz[,1]), max(xyz[,1]), length = 150),
#		 yo=seq(min(xyz[,2]), max(xyz[,2]), length = 150),linear = TRUE, extrap=FALSE, 
#		 duplicate = "error", dupfun = NULL, ncp = NULL)
#	contour(xyz, ...) 
#}
##
## use of coefficients
##
coef.spT<-function(object, digits=4, ...){
   round(t(object$parameter)[1,],digits=digits)
}
##
## use of residuals
##
residuals.spT<-function(object, ...){
   #if(object$combined.fit.pred==TRUE){
   #   stop("\n# Error: not useful for output with combined fit and predict")
   #}
   #else{
     if(object$scale.transform=="NONE"){
     tmp<-object$Y-object$fitted[,1]
     tmp
     }
     else if(object$scale.transform=="SQRT"){
     tmp<-sqrt(object$Y)-object$fitted[,1]
     tmp
     }
     else if(object$scale.transform=="LOG"){
     tmp<-log(object$Y)-object$fitted[,1]
     tmp
     }
     else{
     }
   #}
}
##
## use of formula
##
formula.spT<-function(x, ...){
  x$call
}
##
## use of terms
##
terms.spT<-function(x, ...){
  terms(x$call)
}
##
## For spatially varying beta
##
sp<-function(x){
  class(x)<-"spBeta"
  x
}
##
## For temporally varying beta
##
tp<-function(x){
  class(x)<-"tpBeta"
  x
}
##
## find close locations
##
spT.ObsGridLoc<-function(obsLoc, gridLoc, distance.method="geodetic:km",
          plot=FALSE)
{
  #
      #
      if(!distance.method %in% c("geodetic:km", "geodetic:mile", "euclidean",
        "maximum", "manhattan", "canberra")){
        stop("\n# Error: correctly define distance.method \n")
      }
      #
           coords.all <- rbind(obsLoc, gridLoc)
           n1.sites <- 1:dim(obsLoc)[[1]]
           n2.sites <- (dim(obsLoc)[[1]] + 1):(dim(coords.all)[[1]])
      #
      if(distance.method=="geodetic:km"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(distance.method=="geodetic:mile"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else{
       coords.D <- as.matrix(dist(coords.all, distance.method, diag = T, upper = T))
      }  
      #
      coords.D[is.na(coords.D)]<-0
      diag(coords.D)<-NA
      # 
      fdmis<-coords.D[n1.sites, n2.sites]
      min.fnc<-function(x){
        x<-cbind(x,1:length(x))
        x<-x[order(x[,1]),] 
        x<-x[1,]
        x  
      }
      fdmis<-t(apply(fdmis,1,min.fnc))
      fdmis<-cbind(obsLoc,gridLoc[fdmis[,2],],fdmis[,2],fdmis[,1])
      dimnames(fdmis)[[2]]<-c("obsLon","obsLat","gridLon","gridLat","gridNum","dist")
      if(plot==TRUE){
        plot(fdmis[,1:2],pch="*")
        points(fdmis[,3:4],pch=19,col=2)
        legend("bottomleft",pch=c(8,19),col=c(1,2),legend=c("Observation locations", "Grid locations"),cex=0.9)
      }
      fdmis
}
##
## combine grid data and coords for spTimer input, gridData should be in array with dim = lon, lat, day, year 
##
spT.gridTodata<-function(gridData, gridLoc=NULL, gridLon=NULL, gridLat=NULL)
{
  # gridData should be in array with dim = lon, lat, day, year
  # for morethan one array use list argument
  # gridLoc should be obtain using expand.grid() argument

  options(warn=-1)
  if(!is.null(gridLoc)){
    lonlat<-gridLoc
    dimnames(lonlat)[[2]]<-c("lon","lat")
  }
  else if((!is.null(gridLon)) & (!is.null(gridLat))){ 
    lonlat<-expand.grid(gridLon, gridLat)
    dimnames(lonlat)[[2]]<-c("lon","lat")
  }
  else{
    stop("Error: check grid location input in the function")
  }
  if(is.list(gridData)){
      n<-length(gridData)
      if(n>1){
        dat<-NULL
        for(i in 1:n){
          grid<-c(gridData[[i]])
          dat<-cbind(dat,grid)
        }
        for(i in 1:n){ dimnames(dat)[[2]][i]<-paste("grid",i,sep="") }
        dat<-cbind(lon=rep(lonlat[,1],length(c(gridData[[1]]))/dim(lonlat)[[1]]),
           lat=rep(lonlat[,2],length(c(gridData[[1]]))/dim(lonlat)[[1]]),dat)
      }
      else{
      dat<-cbind(lon=rep(lonlat[,1],length(c(gridData[[1]]))/dim(lonlat)[[1]]),
           lat=rep(lonlat[,2],length(c(gridData[[1]]))/dim(lonlat)[[1]]),grid=c(gridData))
      }
  }
  else{
      dat<-cbind(lon=rep(lonlat[,1],length(c(gridData))/dim(lonlat)[[1]]),
           lat=rep(lonlat[,2],length(c(gridData))/dim(lonlat)[[1]]),grid=c(gridData))
  }
  dat<-dat[order(dat[,1],dat[,2]),]
  dat
}
##
## combine observation and grid data for model fitting 
##
spT.ObsGridData<-function(obsData, gridData, obsLoc, gridLoc, distance.method="geodetic:km")
{
  # obsData should be in data form
  # gridData should be in array form with dim = lon, lat, day, year

  obsLoc<-as.matrix(obsLoc)
  gridLoc<-as.matrix(gridLoc)
  dimnames(obsLoc)<-NULL
  dimnames(gridLoc)<-NULL

  if((is.list(gridData)) & (is.list(gridLoc))){
      n1<-length(gridData)
      n2<-length(gridLoc)
      if(n1 != n2){ stop("\n Error: list of gridData and gridLoc should be same") }
      datt<-NULL
      for(i in 1:n1){
         gLoc<-gridLoc[[i]]
         gLoc<-gLoc[order(gLoc[,1],gLoc[,2]),]
         comb<-spT.ObsGridLoc(obsLoc=obsLoc, gridLoc=gLoc, distance.method=distance.method)
         grid<-spT.gridTodata(gridData=gridData[[i]], gridLoc=gLoc)
         grid<-cbind(gridNum=sort(rep(1:dim(gLoc)[[1]],dim(grid)[[1]]/dim(gLoc)[[1]])),grid)
         dat<-NULL
         for(j in comb[,5]){ dat<-rbind(dat,grid[grid[,1]==j,]) }
         dimnames(dat)[[2]]<-c(paste("gridNum",i,sep=""),paste("grid.lon",i,sep=""),paste("grid.lat",i,sep=""),paste("grid",i,sep=""))
         if(dim(obsData)[[1]] != dim(dat)[[1]]){ stop("Error: check the day and year for gridData\n obsData and gridData time points mismatch.\n")}
         datt<-cbind(datt,dat)
      }
      dat<-cbind(obsData,datt); rm(datt);
      dat
  }
  else{
      comb<-spT.ObsGridLoc(obsLoc=obsLoc, gridLoc=gridLoc[order(gridLoc[,1],gridLoc[,2]),], distance.method=distance.method)
      grid<-spT.gridTodata(gridData=gridData, gridLoc=gridLoc[order(gridLoc[,1],gridLoc[,2]),])
      grid<-cbind(gridNum=sort(rep(1:dim(gridLoc)[[1]],dim(grid)[[1]]/dim(gridLoc)[[1]])),grid)
      dat<-NULL
      for(i in comb[,5]){ dat<-rbind(dat,grid[grid[,1]==i,]) }
      dimnames(dat)[[2]][2:3]<-c("grid.lon","grid.lat")
      if(dim(obsData)[[1]] != dim(dat)[[1]]){ stop("Error: check the day and year for gridData\n   obsData and gridData time points mismatch.\n")}
      dat<-cbind(obsData,dat)
      dat
  }
}
##
## Print function
##
print.spT<-function(x, ...) {
    cat("-----------------------------------------------------"); cat('\n');
    cat("Model: "); cat(x$model); cat('\n');
    cat("Call: "); print(x$call); #cat('\n')
    cat("Iterations: "); cat(x$iterations); cat("\n")
    cat("nBurn: "); cat(x$nBurn); cat("\n")
    cat("Acceptance rate for phi (%): "); cat(x$accept); cat("\n")
    cat("-----------------------------------------------------"); cat('\n');
    #cat("PMCC: "); cat("\n");
    print(x$PMCC); 
    cat("-----------------------------------------------------"); cat('\n');
    #cat("Parameters:\n")
    #print(x$parameter); #cat("\n");
    #cat("-----------------------------------------------------"); cat('\n');
    cat("Computation time: "); cat(x$computation.time); cat("\n")
}
##
## fitted function
##
fitted.spT<-function(object, ...){
    x<-data.frame(object$fitted)
    x
}
##
## use of package forecast
##
as.forecast.object<-function(object, site=1, level=c(80,95), ...){
   # object is the output from the predict.spT for forecast
   x<-NULL
   x$model<-list(Name = object$model)
   x$method<-paste(object$model,"spatio-temporal model, site:",site)
   x$x<-ts(object$obsData)[,site]
   x$fitted<-ts(object$fittedData)[,site]
   x$residuals<-ts(object$residuals)[,site]
   x$mean<-ts(object$Mean,start=c(dim(object$obsData)[[1]]+1,1))[,site]
   x$level<-level
   x$upper<-array(apply(object$fore.samples,1,quantile,probs=c(level/100)),dim=c(length(level),dim(object$Mean)[[1]],dim(object$Mean)[[2]]))
   x$upper<-t(x$upper[,,site])
   x$lower<-array(apply(object$fore.samples,1,quantile,probs=c(1-level/100)),dim=c(length(level),dim(object$Mean)[[1]],dim(object$Mean)[[2]]))
   x$lower<-t(x$lower[,,site])
   cat(paste("Forecast for site:",site),"\n")
   class(x)<-"forecast"
   x
}
##
## use of confint()
##
confint.spT<-function(object, parm, level=0.95, ...){
   #
   x<-as.mcmc(object)
   up<-level+(1-level)/2
   low<-(1-level)/2
   FUN <- function(x){quantile(x,probs=c(low,up))}
   out<-apply(x,2,FUN=FUN)
   out<-t(out)
   if(missing(parm)){
     out
   }
   else{
    if(length(parm)>1){ 
      out<-as.matrix(out[dimnames(out)[[1]] %in% parm,])
	  out
	}
	else{
      out<-t(as.matrix(out[dimnames(out)[[1]] %in% parm,]))
  	  dimnames(out)[[1]]<-parm
	  out
	}
   }
}
##
## Gamma prior
##
Gamm<-function(a=NA,b=NA){
   out<-matrix(c(a,b),1,2)
   class(out)<-"Gamma"
   out
}
##
## Normal prior
##
Norm<-function(mu=NA,sig=NA){
   out<-matrix(c(mu,sig),1,2)
   class(out)<-"Normal"
   out
}
##
## Uniform prior
##
Unif<-function(low=NA,up=NA){
   out<-matrix(c(low,up),1,2)
   class(out)<-"Uniform"
   out
}
##
##
##