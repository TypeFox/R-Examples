dwiMD <- function(object, ...) cat("No Mean Diffusivity calculation defined for this class:",class(object),"\n")

setGeneric("dwiMD", function(object,  ...) standardGeneric("dwiMD"))

setMethod( "dwiMD", "dtiData", function( object, eps=.05) {
   grad <- object@gradient[,-object@s0ind]
   ngrad <- dim(grad[2])
##
##  first get optimal gradient triple 
##
   z <- abs(t(grad)%*%grad)
   zsorted <- sort(as.dist(z))
   ntry <- ngrad*(ngrad-1)/2
   krit <- 3
   l <- 1
   while(krit < eps & l < ntry ){
      ij <- (1:ngrad)[apply(zsorted[l]==z,1,any)]
      zz <- c(1,1)%*%z[ij,]
      k <- (1:ngrad)[zz==min(zz)][1]
      nkrit <- (sum(z[c(ij,k),c(ij,k)])-3)/2
      if(nkrit < krit){
         krit <- nkrit
         ijk <- c(ij,k)
      }
      l <- l+1
   }
##
##  gradient indices in ijk
##
   si <- object@si[,,,(1:object@ngrad)[-object@s0ind][ijk]]
   s0 <- object@si[,,,object@s0ind]
   bv <- object@bvalue[(1:object@ngrad)[-object@s0ind][ijk]]
   if(length(object@s0ind)>1) s0 <- apply(s0,1:3,mean)
   ADC <- sweep(-log1p(sweep(si,1:3,s0,"/")),4,bv,"/")
   ADC[ADC<0] <- 0
   MD <- apply(ADC,1:3,mean)
   MD
}
)

setMethod( "dwiMD", "dtiTensor", function(object) {
   extract(object,what="md")$md
}
)

 


