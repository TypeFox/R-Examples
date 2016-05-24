ellipsoid3d<-function(rx=1,ry=1,rz=1,n=30,ctr=c(0,0,0),
                        qmesh=FALSE,
                        trans = rgl::par3d("userMatrix"),...) {
							
  if(requireNamespace("rgl", quietly=TRUE)==FALSE){stop("rgl not loaded")}
						
 ###################################################################
 ## This function is taken from Daniel Adlers and Duncan Murdochs ##
 ##                         rgl package                           ## 	
 ###################################################################			
 		
  if (missing(trans) && !rgl::rgl.cur()) trans <- diag(4)
  degvec <- seq(0,2*pi,length=n)
  ecoord2 <- function(p) {
    c(rx*cos(p[1])*sin(p[2]),ry*sin(p[1])*sin(p[2]),rz*cos(p[2])) }
  v <- apply(expand.grid(degvec,degvec),1,ecoord2)
  if (qmesh) v <- rbind(v,rep(1,ncol(v))) ## homogeneous
  e <- expand.grid(1:(n-1),1:n)
  i1 <- apply(e,1,function(z)z[1]+n*(z[2]-1))
  i2 <- i1+1
  i3 <- (i1+n-1) %% n^2 + 1
  i4 <- (i2+n-1) %% n^2 + 1
  i <- rbind(i1,i2,i4,i3)
  if (!qmesh)
    rgl::quads3d(v[1,i],v[2,i],v[3,i],...)
  else return(rgl::rotate3d(rgl::qmesh3d(v,i,material=...),matrix=trans))

}


plotsubspace<-function(CA, CB=NULL, corr=FALSE, shadeCA=TRUE, shadeCB=TRUE, axes.lab=FALSE, ...){

  if(requireNamespace("rgl", quietly=TRUE)==FALSE){stop("rgl not loaded")}

  Avec<-eigen(CA)$vectors[,1:3]
  Aval<-eigen(CA)$values
  propA<-round(100*sum(Aval[1:3])/sum(Aval),0)
  Aval<-2*sqrt(Aval[1:3])
  
  rgl::clear3d()

  s1 <- ellipsoid3d(Aval[1], Aval[2],  Aval[3], qmesh = TRUE, trans = diag(4))

  if(dim(CA)[1]==3){
    s1<-rgl::rotate3d(s1,matrix=t(Avec))
  }

  if(shadeCA==TRUE){
    rgl::shade3d(s1, col = "red") 
  }else{
    rgl::wire3d(s1, col = "red")
  } 


  if(is.null(CB)==FALSE){
    if(dim(CA)[1]==3){
      Bvec<-eigen(CB)$vectors
      Bval<-2*sqrt(eigen(CB)$values)
    }else{
      Bproj<-t(Avec)%*%CB%*%Avec
      Bvec<-eigen(Bproj)$vectors
      Bval<-eigen(Bproj)$values
      propB<-round(100*sum(Bval)/sum(eigen(CB)$values),0)
      Bval<-2*sqrt(Bval)
    }	
    s2 <- ellipsoid3d(Bval[1], Bval[2],  Bval[3], qmesh = TRUE, trans = diag(4))
    s2<-rgl::rotate3d(s2, matrix=t(Bvec))

    if(shadeCB==TRUE){
      rgl::shade3d(s2, col = "blue")
    }else{
      rgl::wire3d(s2, col = "blue")
    }
   }

   if(dim(CA)[1]<4){
     xvec<-colnames(CA)[1]
     yvec<-colnames(CA)[2]
     zvec<-colnames(CA)[3]
   }else{
     xvec<-paste("[", paste(round(Avec[,1],2),collapse=","), "]", sep="")
     yvec<-paste("[", paste(round(Avec[,2],2),collapse=","), "]", sep="") 
     zvec<-paste("[", paste(round(Avec[,3],2),collapse=","), "]", sep="")
   }

   nameA<-as.character(substitute(CA))
   nameA<-nameA[length(nameA)]

   if(is.null(CB)){
     nameB<-NULL
   }else{
     nameB<-as.character(substitute(CB))
     nameB<-nameB[length(nameB)]
   }
   rgl::axes3d()
   if(dim(CA)[1]!=3){
      mvec<-paste(paste(nameA, "=", propA , "%", sep=""), paste(nameB, "=", propB, "%",sep=""), sep="    ")
      rgl::title3d(main=mvec)
   }
   if(axes.lab==TRUE){
     rgl::title3d(xlab=xvec, ylab=yvec, zlab=zvec)
   }
} 
