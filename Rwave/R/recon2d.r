#shift fft 1D data
#-----------------

fftshift <- function(x) {
    lng <- length(x)
    y <- numeric(lng)
    y[1:(lng/2)] <- x[(lng/2+1):(lng)]
    y[(lng/2+1):lng] <- x[1:(lng/2)]
    y
}


#b727 <- matrix(scan("b727s.dat"),ncol=512,byrow=TRUE)
#b727 <- matrix(scan("b727r.dat"),ncol=256,byrow=TRUE)

#ncol <- dim(b727)[2]/2
#b727r<- b727[,-1+2*(1:ncol)]   
#b727i<- b727[,2*(1:ncol)]       

#B727 <- t(Conj(complex(real=b727r, imag=b727i)))
#nrow <- dim(B727)[1]
#ncol <- dim(B727)[2]

#Generating the original image
#-----------------------------
#* apply fftshift twice to swap the center of the matrix


#o727 <- t(apply(apply(B727,2,fftshift),1,fftshift))
#o727 <- apply(o727,2,fft)
#o727 <- t(apply(apply(o727,2,fftshift),1,fftshift))



#*** only for display to compare with matlab output
#image(Mod(t(o727)))

#Applying Gabor Transform to ...
#-------------------------------
#gtime <- 100
#scale <- 50

#BB727 <- B727
#cgt727 <- array(0+0i,c(ncol,nrow,gtime))
#for(k in 1:ncol) cgt727[k,,] <-
#complex(real=cgt(Re(BB727[,k]),gtime,2/gtime,scale,plot=FALSE),imag=cgt(Im(BB727[,k]),gtime,2/gtime,scale,plot=FALSE))

#image(apply(apply(Mod(cgt727),c(1,3),mean),1,fftshift))


#* The following tow S-routine concludes our Radar experiments
#-------------------------------------------------------------
#b727 <- matrix(scan("b727s.dat"),ncol=512,byrow=TRUE)
#b727 <- matrix(scan("b727r.dat"),ncol=256,byrow=TRUE)

#Generating the original image
#-----------------------------

#e.g. B727 <- showRadar(b727)

showRadar<- function(x) {
   ncol <- dim(x)[2]/2
   xr<- x[,-1+2*(1:ncol)]   
   xi<- x[,2*(1:ncol)]       

   y <- t(Conj(xr + 1i*xi))

   oy <- t(apply(apply(y,2,fftshift),1,fftshift))
   oy <- apply(oy,2,fft)
   oy <- t(apply(apply(oy,2,fftshift),1,fftshift))
   image(Mod(t(oy)))
   oy
}

#Enhancing Radar Image by Gabor Transform
#----------------------------------------


#gtime <- 128
#scale <- 50

#e.g. B727 <- cgtRadar(b727,gtime,scale)
#e.g. 
#     ncol <- dim(b727)[2]/2
#     b727r<- b727[,-1+2*(1:ncol)]   
#     b727i<- b727[,2*(1:ncol)]       

#     b727 <- t(Conj(complex(real=b727r, imag=b727i)))
#     B727 <- cgtRadar(b727,gtime,scale,flag=FALSE)

cgtRadar <- function(x,gtime,scale,flag=TRUE) {

   y <- x

   if(flag) {
     ncol <- dim(x)[2]/2
     xr<- x[,-1+2*(1:ncol)]   
     xi<- x[,2*(1:ncol)]       

     y <- t(Conj(xr + 1i*xi))
   }

   nrow <- dim(y)[1]
   ncol <- dim(y)[2]
   cgty <- array(0+0i,c(ncol,nrow,gtime))
   for(k in 1:ncol) cgty[k,,] <-
     cgt(Re(y[,k]),gtime,2/gtime,scale,plot=FALSE) + 1i*cgt(Im(y[,k]),gtime,2/gtime,scale,plot=FALSE)

   oy <- apply(apply(Mod(cgty),c(1,3),mean),1,fftshift)
   for(k in 1:nrow) cgty[,k,] <- apply(cgty[,k,],1,fftshift)
   image(oy)
   list(output=oy,cgtout=cgty)
}






