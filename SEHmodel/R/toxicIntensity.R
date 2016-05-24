# toxicIntensity.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

# @import MASS
#' @importFrom pracma trapz
#' @importClassesFrom raster RasterLayer
#' @importFrom raster raster
#' @importFrom raster extent
#' @importFrom raster extent<-
#' @importFrom raster rasterize
#' @importFrom raster as.matrix
# @importMethodsFrom raster as.vector
#' @importFrom raster projectRaster
#' @importFrom raster addLayer
#' @importFrom raster writeRaster
#' @importFrom fftwtools fftw_r2c_2d
#' @importFrom fftwtools fftw2d
#' @importFrom mvtnorm pmvnorm

#' @title toxicIntensity Method
#' 
#' @description Simulate contaminants intensity over the landscape by two steps : dispersal of toxic particules and local intensity of particules after dispersal.
#' 
#' @name toxicIntensity
#' @param objectL A Landscape object
#' @param ... parameters
#' @rdname Landscape-toxicIntensity-method
#' @aliases toxicIntensity,Landscape-method
#' @exportMethod toxicIntensity
setGeneric(name="toxicIntensity",
           def=function(objectL,...)
             standardGeneric("toxicIntensity")
)

#' @name toxicIntensity
#' @details The dispersal of contaminants is implemented by rastering the landscape and by computing the convolution between sources emissions and a dispersal kernel.
#' 
#' The dispersion kernel by default is Normal Inverse Gaussian kernel ("NIG" function). Currently, two others are implemented "geometric" (with parameter \code{a}) and "2Dt" kernels (with parameters \code{a}, \code{b}, \code{c1}, \code{c2}).
#' 
#' Local intensity depends of \code{beta} and \code{alpha} parameters. Beta represents the toxic adherence between [0,1].
#' Alpha represents a list of parameters of the lost of toxic particules due to covariates (precipitation).
#' There are two configurations to integrate the loss in the function : 
#' (i) simulating covariate (simulate=TRUE) or (ii) uploading covariate (simulate=FALSE).
#' The covariate is linked to the loss by a linear regression with paramaters minalpha, maxalpha, covariate_threshold.
#' 
#' @param toxic_emission Matrix of sources emissions, row as sources ID, col as time
#' @param mintime Start simulation time (default=1)
#' @param maxtime End simulation time
#' @param size_raster raster size (default = 2^10)
#' @param kernel dispersion kernel, function name (default = NIG)
#' @param kernel.options parameters list for the kernel function
#' @param beta toxic adherence parameter between 0 and 1 (default = 0.4)
#' @param alpha list of toxic loss options
#' 
#' (default = list(minalpha=0.1,maxalpha=0.95,covariate_threshold=30,simulate=TRUE,covariate=NULL))
#' @return A ToxicIntensityRaster, a 3D array as time matrix dispersion, [t,x,y]
#' @aliases toxicIntensity,Landscape-method
#' @rdname Landscape-toxicIntensity-method
setMethod(f="toxicIntensity",
          signature="Landscape",
          definition=function(objectL,toxic_emission,mintime=1,maxtime=60,size_raster=2^10,kernel="NIG",kernel.options=list("a1"=0.2073 ,"a2"=0.2073 , "b1"=0.3971, "b2"=0.3971, "b3"= 0.0649, "theta"=0),beta=0.4,alpha=list(minalpha=0.1,maxalpha=0.95,covariate_threshold=30,simulate=T,covariate=NULL)) {
            sources<-getSPSources(objectL)
            if(length(sources) != nrow(toxic_emission)) {
              cat("srcdist vector length differ from landscape sources number")
              return(NULL)
            }
            
            if(ncol(toxic_emission) != maxtime) {
              cat("srcdist col number differ from maxtime")
              return(NULL)
            }
            
            if(size_raster >= 2^13) cat("WARNING : raster_size may be too large for memory usage RAM ")
            
            #cell surface
            S<-(objectL@xmax-objectL@xmin)*(objectL@ymax-objectL@ymin)/(size_raster*size_raster)
            
            #--------------------------------------
            # perte sur les feuilles de plante hôte : dépend du temps et des conditions météo
            # fonction simulateur stochastique de climat :(attention chemin vers data de climat dans la fonction)
            ayear<-seq(from=as.Date("01/07/2014","%d/%m/%Y"),to=as.Date("30/06/2015","%d/%m/%Y"),by="day")
                   
            if(alpha$simulate==T) {
              precip=simul.precipitation(starttime="01/07",endtime=format(ayear[maxtime],"%d/%m"),data=alpha$covariate)
            } else {
              precip=alpha$covariate
            }

            loss=numeric(maxtime)
            loss[precip==0]=alpha$minalpha
            pente=(alpha$maxalpha-alpha$minalpha)/alpha$covariate_threshold
            for(p in which(precip>0&precip<=alpha$covariate_threshold)) {
              loss[p]=precip[p]*pente + alpha$minalpha
            }
            loss[precip>alpha$covariate_threshold]=alpha$maxalpha
                       
            #rasterize
            r<-raster(ncol=size_raster,nrow=size_raster,ext=extent(sources),crs=proj4string(objectL@thelandscape))
            #extent(r)<-extent(sources)
            
            # raster avec ID parcelle sources
            rb<-rasterize(sources,r,field=as.numeric(row.names((sources@data)))) 
            rm(r)
            
            convol.env <- new.env(parent = as.environment("package:SEHmodel"))
            #----------------------------------------------
            # flux.convol : fait le produit de convolution pour calculer les flux (dispersion avec Emission de maïs=1) dans le domaine
            # calcul d'une matrice de convolution par parcelle OGM :
            size_domain<-objectL@xmax-objectL@xmin
            
            for(i in row.names(sources)) {
              rb_temp<-rb
              #rb_temp[which(rb_temp@data@values!=i)]<-NA
              rb_temp[!rb_temp@data@values%in%as.numeric(i)] <- NA
              rb_temp[rb_temp@data@values%in%as.numeric(i)] <- 1
              #rb_temp[which(rb_temp@data@values==i)]<-1
              assign(paste0("field",i),flux.convol(kernel,kernel.options,size_domain,size_raster,rb_temp,convol.env),envir = convol.env)
              rm(rb_temp)
            }
            
            res=array(0,c(maxtime,size_raster,size_raster))
            for(t in mintime:maxtime) {
              #all_toxic=array(0,c(size_raster,size_raster))
              all_toxic <- matrix(data = 0,nrow = size_raster,ncol = size_raster)
              for(i in row.names(sources)) {
                if(toxic_emission[as.character(i),paste0("t.",t)] != 0) {
                  all_toxic <- all_toxic + (S * toxic_emission[as.character(i),paste0("t.",t)] * as.matrix(get(paste0("field",i),envir = convol.env)) * beta)
                }
              }
              if(t == mintime) {res[t,,] <- all_toxic
              } else { res[t,,] <- all_toxic + res[t-1,,]*(1-loss[t]) }
            }
            
            rm(convol.env)
            
            attr(res,"class") <- "ToxicIntensityRaster"
            return(res)
          }
)


#' PlotLandscapetoxicIntensity Method
#' 
#' Plot a landscape and the toxic intensity at a time t.
#' 
#' @name plotLandscapetoxicIntensity
#' 
#' @param objectL a Landscape object
#' @param objectT a 3D array from toxicIntensity, [time,x,y]
#' @param time time to plot the toxic dispersion
#' @export
plotLandscapetoxicIntensity<-function(objectL,objectT,time) {
  plot(objectL)
  p<-heat.colors(100, alpha = 0.6)
  p[100]=rgb(0,0,0,alpha=0)
  temp<-objectT[time,,]
  temp[which(temp<=0)] <- NA
  r<-raster(as.matrix(temp),crs=CRS("+init=epsg:2154"))
  extent(r)<-extent(objectL@xmin,objectL@xmax,objectL@ymin,objectL@ymax)
  if( !is.na(proj4string(objectL@thelandscape)) ) { r<-projectRaster(r,crs=proj4string(objectL@thelandscape)) }
  raster::image(r,col=p[length(p):1],useRaster=F,add=T,bg="transparent")
  # fast way but not in coordinates mapping
  #image(x=seq(x@thelandscape@bbox[1,1],x@thelandscape@bbox[1,2],(x@thelandscape@bbox[1,2]-x@thelandscape@bbox[1,1])/nrow(temp)),y=seq(x@thelandscape@bbox[2,1],x@thelandscape@bbox[2,2],(x@thelandscape@bbox[2,2]-x@thelandscape@bbox[2,1])/ncol(temp)),z=t(as.matrix(temp))[1:nrow(temp),ncol(temp):1],col=p[length(p):1],useRaster=T,add=T,bg="transparent")
  
  fields::image.plot(as.matrix(temp),legend.only=T,smallplot=c(0.85,0.88,0.20,0.95),col=p[length(p):1])
  mtext(text ="Toxic intensity",line = 0,side = 3,adj = 1.1,padj = 1)
}

# plottoxicIntensity<-function(td,mintime=1,maxtime) {
#   
#   
#   pdf(file="toxicIntensity.pdf")
#   par(mfrow=c(2,2),mar=c(1,1,1,1))
#   #par(mar=rep(2,4))
#   for(t in mintime:maxtime) {
#     image(as.matrix(td[t,,]))
#   }
#   dev.off()
# }

############## PRIVATE ####################

# Convolution function
flux.convol<-function(kernel,kernel.options,size_domain,size_raster,rb,envir){ 
  
  if(!exists("z",envir=envir,inherits = FALSE) | is.null(envir$z)) {
  
    if(size_domain%%2 == 1) {
      size_domain <- size_domain+1
    }
  
    x=y=seq(from=-(size_domain/2),to=size_domain/2,length.out=size_raster)
    #z<-kernel
    z<-outer(x,y,kernel,kernel.options)
    envir$z<-z/trapz2d(z) #on normalise le noyau
    rm(x);rm(y);rm(z);
  }
  e<-as.matrix(rb)
  e[is.na(e)]<-0
  F<-convolution(envir$z,e,size_raster)
  
  #rm(e);gc()
  return (F)
}

############ kernel for pollen dispersal

#paramètres, distance entre l'émission du pollen (maïs) et la récéption (plantes hôtes),hyp: h=2m
#pas de vent: mu=teta=0
## noyau Normal Inverse Gaussian model (MAPOD, Klein et al 2003)
# NIG<-function(x,y,h=2,mu=0,teta=0){
#   
#   lambda_z=0.027*h/0.831
#   lambda_x=0.165*h*mu*cos(teta)/(2*0.831)
#   lambda_y=0.165*h*mu*sin(teta)/(2*0.831)
#   delta_x=0.499*0.831/h
#   delta_y=delta_x
#   
#   p=lambda_z^2+lambda_x^2+lambda_y^2
#   q=1+delta_x^2*x^2+delta_y^2*y^2
#   
#   A=delta_x*delta_y*exp(lambda_z)/(2*pi)
#   B=(q^-0.5+p^0.5)/q
#   C=exp(-sqrt(p*q))*exp(lambda_x*delta_x*x+delta_y*lambda_y*y)
#   return(A*B*C)
# }

NIG<-function(x,y,kernel.options=list("a1"=0.2073 ,"a2"=0.2073 , "b1"=0.3971 ,"b2"=0.3971 ,"b3"= 0.0649,"theta"=0)){
  if(is.null(kernel.options$theta)) { theta=0 }
  else { theta=kernel.options$theta }
  lambda_z=kernel.options$b3                 #=0.027*h/0.831 avec h=2
  lambda_x=kernel.options$b1*cos(theta)      #0.165*h*mu*cos(theta)/(2*0.831) # avec h=mu=2
  lambda_y=kernel.options$b2*sin(theta)     #0.165*h*mu*sin(theta)/(2*0.831)  avec h=mu=2
  delta_x=kernel.options$a1                  #0.499*0.831/h
  delta_y=kernel.options$a2
  p=lambda_z^2+lambda_x^2+lambda_y^2
  q=1+delta_x^2*x^2+delta_y^2*y^2
  A=delta_x*delta_y*exp(lambda_z)/(2*pi)
  B=(q^-0.5+p^0.5)/q
  C=exp(-sqrt(p*q))*exp(lambda_x*delta_x*x+delta_y*lambda_y*y)
  return(A*B*C)
}

#2Dt Student kernel

student<-function(x,y,kernel.options=list("c1"=1.12,"a"=1.55,"b"=1.45,"c2"=0,"theta"=0)){
  a=kernel.options$a
  b=kernel.options$b
  c2=kernel.options$c2
  c1=kernel.options$c1
  theta=kernel.options$theta
  A=(b-1)/(pi*a*a)
  B=(1+(x*x+y*y)/(a*a))^(-b)
  C=exp(c1*cos(theta-c2))
  return(A*B*C)
}

#Geometric kernel Soubeyrand

geometric<-function(x,y,kernel.options=list("a"=-2.59))   {
  aa=kernel.options$a
  A=(1+sqrt(x*x+y*y))^aa;B=(aa+1)*(aa+2)/(2*pi);return(A*B)}


#### FatTail 
FatTail<-function(x,y) {
  
  return(sqrt(x^2+y^2)^-2)
}

#-------------------------------------
#Hoffman 2014, heavy tail power dispersal kernel
kernel_fat_tail<-function(x,y){res<-1.271e6*sqrt(x*x+y*y)^-0.585;return(res)}

#-------------------------------------
## convolution : convolution with FFT
#kernel: noyau, emission: répartition de la masse,s taille de la matrice (cf package fftwtools)

convolution<-function(kernel,emission,s){
  fp=fftw_r2c_2d(kernel)
  fz=fftw_r2c_2d(emission)
  f1=fz*fp
  p1=fftw2d(f1,inverse=1) #inverse du produit
  p<-Re(p1)/(s*s) #on prend la partie réelle et on renormalise par la taille de la matrice
  p<-shift_fft(p)
  return (p)
}

# Trapz2D, 2D quadrature with trapeze method
trapz2d<-function(p){res=sum(apply(p,2,trapz)); return(res)}

# remet les données en forme pour la convolution en 2D via FFT
shift_fft<-function(m){
  i<-dim(m)[1]
  j<-dim(m)[2]
  n<-matrix(rep(0,i*j),nrow=i,ncol=j)
  i21<-i/2+1
  j21<-j/2+1
  i2<-i/2
  j2<-j/2
  
  n<-m[c(i21:i,1:i2),c(j21:j,1:j2)]
  
  return(n)
}
