`syrjala` <-
function(coords=NULL,var1=NULL,var2=NULL,nperm=999){
   datanames <- c(deparse(substitute(var1)), 
                     deparse(substitute(var2)))
   x <- coords$x
   y <- coords$y
   if (is.null(x) || is.null(y)){
      stop("`coords` contains no $x or $y component.")
   }
   if (any(is.na(x)) || any(is.na(y))){
      stop("`coords` contains missings.")
   }
   nd <- length(x)
   if (length(y) != nd){
      stop("length of `coords$x` and `coords$y` differ.")
   }
   if (length(var1) != nd){
      stop("number of `coords` and `variable 1` differ.")
   }
   if (length(var2) != nd){
      stop("number of `coords` and `variable 2` differ.")
   }
   na <- which(is.na(var1) | is.na(var2))
   if (length(na) > 0 ){
      x <- x[-na]
      y <- y[-na]
      var1 <- var1[-na]
      var2 <- var2[-na]
   }
   if(any(var1 < 0)||any(var2 < 0)){
      stop('negative data found!')
   }
   if (all(var1 == var2)){
      stop('all `variable 1` and `variable 2` are equal!')
   }
   possible <- duplicated(x) & duplicated(y)
   if (any(possible)){
      result <- possible
      xvals <- unique(x[possible])
      for (xvalue in xvals) {
         subs <- (x == xvalue)
         result[subs] <- duplicated(y[subs])
      }
      if(any(result)){
         result <- which(result)
         if(length(result)>1){
            stop(paste('points',result,'are location duplicates!'))
         } else {
            stop(paste('point',result,'is a location duplicate!'))
         }
      }
   }
#    subroutine syrjala(x,y,var1,var2,nd,nperm,cvm,ks)
   ans <- .Fortran('syrjala',x=as.double(x),y=as.double(y),
            var1=as.double(var1),var2=as.double(var2),
            nd=as.integer(nd),nperm=as.integer(nperm),
            cvm=double(nperm+1),ks=double(nperm+1),            
            PACKAGE='ecespa')
   ans <- list(cvm.obs=ans$cvm[1],cvm.sim=ans$cvm,ks.obs=ans$ks[1],ks.sim=ans$ks,datanames=datanames,nperm=nperm)
   class(ans) <- 'syrjala.test'
   return(ans)
}

