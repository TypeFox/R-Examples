register.newfd <- function(yfd, Wfd,type=c('direct','monotone','periodic'))
{
coef  <- Wfd$coefs
shift <- coef[1]

if (all(coef == 0)) {
   if (type=='periodic') {
      if (shift == 0) {
         yregfd <- yfd
         return(yregfd)
      }
   } else {
      yregfd <- yfd
      return(yregfd)
   }
}

# Now evaluate on a fine grid:

Wrange = Wfd$basis$range
ybasis = yfd$basis
yrange = ybasis$range

if( type=='periodic' & any(Wrange != yrange) )
  stop('Registration functions and functions to be registered must have the same range')


neval    <- max(10*ybasis$nbasis + 1,101)
tfine  <- seq(yrange[1],yrange[2], length=neval)  

# Shift registration is easy

if( type=='periodic' ){ 
  yfine = eval.fd(tfine,yfd)
  yfine <- shifty(tfine,yfine,shift)
  ycoef  <- project.basis(yfine, tfine, ybasis, 1)
  yregfd <- fd(ycoef, ybasis)
  return(yregfd)
}

# On the other hand, if we have warping functions:

if( type=='direct') xfine = eval.fd(tfine,Wfd)

if( type=='monotone'){ 
  xfine = eval.monfd(tfine,Wfd)
  xfine = xfine%*%diag(1/xfine[neval,])*(Wrange[2]-Wrange[1])+Wrange[1]
}
xfine = xfine*( xfine>yrange[1] & xfine< yrange[2]) + yrange[2]*(xfine>=yrange[2]) + yrange[1]*(xfine<=yrange[1])
yfine = eval.fd(tfine,yfd)

xdim = dim(xfine)
ydim = dim(yfine)

# Check that we have the right dimensions, we can register multiple y dimensions
# to one warping function, but we must have as many warping functions as there
# are y replicates

if( xdim[2] != ydim[2] ) stop('There must be as many warping function replicates as y replicates')
if( length(ydim) == 3 & length(xdim)==2 ) xfine = array(xfine,ydim)


# Now do the registration

ycoef = 0*yfd$coef
cdim = dim(ycoef)

if(length(cdim)==1) {
  yfine = eval.fd(xfine,yfd)
  ycoef = project.basis(yfine,tfine,ybasis)
}

if(length(cdim)==2){
  for(i in 1:cdim[2]){
    yfine = eval.fd(xfine[,i],yfd[i])
    ycoef[,i] = project.basis(yfine,tfine,ybasis)
  }
}
    
if(length(cdim)==3){
  for(j in 1:cdim[3]){
    for(i in 1:cdim[2]){
       yfine = eval.fd(xfine[,i,j],yfd[i,j])  
       ycoef[,i,j] = project.basis(yfine,tfine,ybasis)
    }
  }
}

yregfd <- fd(ycoef, ybasis)
return(yregfd)
}