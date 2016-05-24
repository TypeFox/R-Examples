#' Generate random points on the sphere
#' @param n Nnumber of data points on the sphere
#' @return latitudes Latitudes of n random points on the sphere
#' @return longitudes Longitudes of n random points on the sphere
#' @export
#' @seealso \code{\link{sphere_khat},\link{sphere_montekhat}}
sphere_random<-function(n){
  phi<-acos(1-2*runif(n))-pi/2 #radius for latitudes
  theta<-2*pi*runif(n)-pi #radius for longitudes
  latr<-phi*180/pi
  lonr<-theta*180/pi
  list(latitudes = latr, longitudes = lonr)
}
#' Latitude-longitude grids
#' @description Generates widely used Latitude-longitude grids for structuring global-scale data. A 10 * 10 equal-area grid is created in the example for future point-pattern analysis.
#' @param degr A number that is applied to create a degr * degr grid, where 2.5 * 2.5 or 5 * 5 or 10 * 10 grids are commonly used.
#' @return latitudes Latitudes of degr * degr grid points
#' @return longitudes Longitudes of degr * degr grid points
#' @export
#' @seealso \code{\link{sphere_montekhat}, \link{sphere_khat}}
sphere_grid<-function(degr){
  grid_lon<-function(degr){
    nlon<-(180*2/degr)
    lon<-matrix(0,1,nlon)
    for(i in 1:nlon){
      lon[i]<--180+degr/2+(i-1)*degr
    }
    return(lon)
  }
  grid_lat<-function(degr){
    nlat<-(90*2/degr)
    lat<-matrix(0,1,nlat)
    for(i in 1:nlat){
      lat[i]<--90+degr/2+(i-1)*degr
    }
    return(lat)
  }
  nlat<-length(grid_lat(degr))
  nlon<-length(grid_lon(degr))
  lonm<-matrix(0,nlat,nlon)
  for (i in 1:nlon){
    lonm[,i]<-grid_lon(degr)[1,i]
  }
  latm<-matrix(0,nlat,nlon)
  for (i in 1:nlon){
    latm[,i]<-grid_lat(degr)[1,]
  }
  return(list(latitudes=latm,longitudes=lonm))
}
#' Calculate spherical K-function
#' @description Main function to obtain spherical K-function for point-pattern analysis on the sphere
#' @param latitudes Latitudes of observed points on the sphere in degrees
#' @param longitudes Longitudes of observed points on the sphere in degrees
#' @param dis Vector of values for the argument r (from 0 to pi), at which K(r) is evaluated. By default, dis = seq(from=0,to=pi,by=0.1).
#' @return Khats Estimated values of K-function
#' @export
#' @seealso \code{\link{sphere_montekhat}}
#' @examples
#' lat<-sphere_random(100)$latitudes; lon<-sphere_random(100)$longitudes
#' d<-seq(from=0,to=pi,by=0.1)
#' sphere_khat(lat,lon,d)
sphere_khat<-function(latitudes, longitudes, dis){
    if (missing(dis)) dis <- seq(0, pi, length.out = 30)
    r<-length(latitudes)
    ndis<-length(dis)
    dis[ndis]<-pi
    CSR<-2*pi*(1-cos(dis))
    phi<-latitudes*pi/180
    theta<-longitudes*pi/180
    delta<-rep(0,r*r)
    sinPhi <- sin(phi) %*% t(sin(phi))
    cosPhi <- cos(phi) %*% t(cos(phi)) * cos(theta %*% t(rep(1, length(theta))) - t(theta %*% t(rep(1, length(theta)))))
    u <- sinPhi + cosPhi
    u[u <= -1] <- -1
    u[u >= 1] <- 1
    diag(u)=0
    delta<-acos(u)
    diag(delta)<-2*pi
    delta<-matrix(delta,ncol=1)
    khat<-rep(0,ndis)
    for (k in 1:ndis){
        khat[k]<-sum(delta<=dis[k])
    }
    khats<-4*pi*khat/(r*(r-1))
    Khats<-khats-CSR
    #   list(khats = khats,Khats = Khats)
    return(Khats)
}
#' K-functions under complete spatial randomness (CSR) by Monte Carlo tests
#' @description Monte Carlo confidence intervals of K-functions under CSR are provided for point-pattern analysis.
#' @param n Number of observed points 
#' @param nsim Number of simulations for K-function
#' @param dis Vector of values for the argument r (from 0 to pi), at which K(r) is evaluated. By default, dis = seq(from=0,to=pi,by=0.1).
#' @return Kci Simulated K-functions under CSR
#' @export
#' @seealso \code{\link{sphere_khat}}
#' @examples
#' #Spherical K function (minus CSR) with 99% confidence intervals 
#' #for point patterns associated with 100 random points
#'
#' sphererandom<-sphere_random(100)
#' latrd<-sphererandom$latitudes;lonrd<-sphererandom$longitudes
#' d<-seq(from=0,to=pi,by=0.1)
#' nd<-length(d)
#' d[nd]<-pi
#' khatrd<-sphere_khat(latrd,lonrd,d)
#' Kcird<-sphere_montekhat(100,100,d)
#' plot(d,khatrd,type='n', ylim=c(-0.4,0.4),xlim=c(0,pi),xaxt = "n",
#' ylab = expression(K - CSR),xlab = expression("Spherical Angle"))
#' axis(1, at = c(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi),
#' labels = expression(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi))
#' polygon(c(d, rev(d)), c(Kcird[1,], rev(Kcird[99,])),col = "grey79", border = FALSE)
#' lines(d,khatrd,col = 4, lwd=2)
#' lines(y=c(0,0),x=c(0,pi),type='l',lty=2,lwd=2)
#' 
#' #Spherical K function (minus CSR) with 99% confidence intervals 
#' #for point patterns associated with 10 * 10 latitude-longitude grid
#'
#' spheregrid<-sphere_grid(25)
#' latm<-as.vector(spheregrid$latitudes)
#' lonm<-as.vector(spheregrid$longitudes)
#' d<-seq(from=0,to=pi,by=0.2)
#' nd<-length(d)
#' d[nd]<-pi
#' khatsg<-sphere_khat(latm,lonm,d)
#' Kcisg<-sphere_montekhat(98,50,d)
#' plot(d,khatsg,type='n', ylim=c(-0.4,0.4),xlim=c(0,pi),xaxt = "n",
#' ylab = expression(K - CSR),xlab = expression("Spherical Angle"))
#' axis(1, at = c(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi),
#' labels = expression(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi))
#' polygon(c(d, rev(d)), c(Kcisg[1,], rev(Kcisg[99,])),col = "grey79", border = FALSE)
#' lines(d,khatsg,col = 4, lwd=2)
#' lines(y=c(0,0),x=c(0,pi),type='l',lty=2,lwd=2)
#' 
#' #Spherical K function (minus CSR) with 95% confidence intervals 
#' #for point patterns associated with 272 points global hexagonal grid
#' 
#' data(Hex272)
#' lath272<-Hex272[,3]
#' lonh272<-Hex272[,2]
#' d<-seq(from=0,to=pi,by=0.2)
#' nd<-length(d)
#' d[nd]<-pi
#' khatsh272<-sphere_khat(lath272,lonh272,d)
#' Kcih272<-sphere_montekhat(272,100,d)
#' plot(d,khatsh272,type='n', ylim=c(-0.15,0.1),xlim=c(0,pi),xaxt = "n",
#' ylab = expression(K - CSR),xlab = expression("Spherical Angle"))
#' axis(1, at = c(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi),
#' labels = expression(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi))
#' polygon(c(d, rev(d)), c(Kcih272[3,], rev(Kcih272[97,])),col = "grey79", border = FALSE)
#' lines(d,khatsh272,col = 4, lwd=2)
#' lines(y=c(0,0),x=c(0,pi),type='l',lty=2,lwd=2)
#' 
#' #Spherical K function (minus CSR) with 95% confidence intervals 
#' #for point patterns associated with 172 upper-air monitoring stations points
#' 
#' data(GUAN)
#' latg<-GUAN[,4]
#' long<-GUAN[,5]
#' d<-seq(from=0,to=pi,by=0.1)
#' nd<-length(d)
#' d[nd]<-pi
#' khatsg<-sphere_khat(latg,long,d)
#' Kcig<-sphere_montekhat(172,100,d)
#' plot(d,khatsg,type='n', ylim=c(-0.1,0.15),xlim=c(0,pi),xaxt = "n",
#' ylab = expression(K - CSR),xlab = expression("Spherical Angle"))
#' axis(1, at = c(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi),
#' labels = expression(0,pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi))
#' polygon(c(d, rev(d)), c(Kcig[3,], rev(Kcig[97,])),col = "grey79", border = FALSE)
#' lines(d,khatsg,col = 4, lwd=2)
#' lines(y=c(0,0),x=c(0,pi),type='l',lty=2,lwd=2)
sphere_montekhat<-function(n,nsim,dis){
    if (missing(dis)) dis<-seq(from=0,to=pi,by=0.1)
    if (missing(nsim)) nsim<-1000
    ndis<-length(dis)
    dis[ndis]<-pi
    Kci<-matrix(data=0,nrow=ndis,ncol=nsim)
    for (i in 1:nsim){
        randompoints<-sphere_random(n)
        latr<-randompoints$latitudes
        lonr<-randompoints$longitudes
        Kci[,i]<-sphere_khat(latr,lonr,dis)
    }
    Kci<-apply(Kci,1,sort)
    return(Kci)
}