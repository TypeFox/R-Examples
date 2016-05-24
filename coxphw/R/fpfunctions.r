        powM2 <- function(z) z^(-2)
        powM1 <- function(z) z^(-1)
        powM0.5 <- function(z) z^(-0.5)
        pow2 <- function(z) z^2
        pow3 <- function(z) z^3


        ## define repeated powers
        RI <- function(z) z * log(z)
        RpowM2 <- function(z) z^(-2) * log(z)
        RpowM1 <- function(z) z^(-1) * log(z)
        RpowM0.5 <- function(z) z^(-0.5) * log(z)
        Rlog <- function(z) log(z) * log(z)
        Rsqrt <- function(z) sqrt(z) * log(z)
        Rpow2 <- function(z) z^2 * log(z)
        Rpow3 <- function(z) z^3 * log(z)

        ## pretransformation function
        PT <- function(z) {
          obj <- fp.scale(z)
         (z + obj$shift) / obj$scale
        }

fp.power<-function(z,a,b=NULL){
  if(!is.null(b)){
     if(b>a) {
       he<-b         # save b
       b<-a          # let new b = a
       a<-he         # let new a = saved b
       }
    if(!(b %in% c(-2,-1, -0.5, 0, 0.5, 1,2,3))) stop("wrong power, must be one of   c(-2,-1, -0.5, 0, 0.5, 1,2,3)\n")
  }
  if(!(a %in% c(-2,-1, -0.5, 0, 0.5, 1,2,3))) stop("wrong power, must be one of   c(-2,-1, -0.5, 0, 0.5, 1,2,3)\n")
  l<-length(z)
  if(a==0) v1<-log(z)
  else v1<-z^a
  v<-v1
  if(!is.null(b)) {
     if(b!=a) {
      if(b==0) v2<-log(z)
      else v2<-z^b
      }
      else v2<-log(z) * v1
      v<-c(v1,v2)
   }

    v <- matrix(v,l,length(v)/l)
    if(!is.null(b)) dimnames(v)[[2]]<-c(paste("fp(",a,")", sep=""), paste("fp(",if(b==a)paste("R"),b,")", sep=""))
    else dimnames(v)[[2]]<-list(paste("fp(",a,")", sep=""))
    dimnames(v)[[1]]<-z
     v
}
