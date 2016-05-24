##################################################################
##########     Print output                  #####################

mix.print <- function(model,digits=3, ...){

  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Class of family",class(model),"not recognized.",sep=" "))
  if (is.list(model$mu)) stop("The mix.print function is only appropriate for the univariate analysis.\n")
  
   sdmu <- NULL
   sdsigma <- NULL
   sdshape <- NULL
   
cat("\n Number of observations: ",model$n, "\n\n")

if((class(model) != "Skew.normal") && (class(model) != "Normal"))	cat("\n Hyperparameter(nu): ",model$nu, "\n\n")

if((class(model) == "Normal") || (class(model) == "t")){
   
   if(length(model$im.sdev) != 0 ) {
      size <- length(model$im.sdev)
      posmu <- seq(1, size, by=3)
      possig <- seq(2, size, by=3)
      sdmu <- model$im.sdev[posmu]
      sdsig <- model$im.sdev[possig]
     estimates <-  rbind(model$mu,sdmu,model$sigma2,sdsig)
   }

   if(length(model$im.sdev) == 0 )   estimates <-  rbind(model$mu,model$sigma2)

   nam <- NULL
   for(i in 1:length(model$mu)) nam <- c(nam,paste("group",i))
   dimnames(estimates)[[2]] <- nam
   if(length(model$im.sdev) != 0 ) dimnames(estimates)[[1]] <- c("mu","sd","sigma2","sd")
   if(length(model$im.sdev) == 0 ) dimnames(estimates)[[1]] <- c("mu","sigma2")
}
if((class(model) != "Normal") && (class(model) != "t")){
   if(length(model$im.sdev) != 0 ) {
        size <- length(model$im.sdev)
        if(class(model) != "Skew.normal"){
            size <- size-1
            if(class(model) == "Skew.cn") size <- size-1
        }
        posmu <- seq(1, size, by=4)
        possig <- seq(2, size, by=4)
        possh <- seq(3, size, by=4)
        sdmu <- model$im.sdev[posmu]
        sdsig <- model$im.sdev[possig]
        sdsh <- model$im.sdev[possh]
        estimates <-  rbind(model$mu,sdmu,model$sigma2,sdsig,model$shape,sdsh)
   }

   if(length(model$im.sdev) == 0 )   estimates <-  rbind(model$mu,model$sigma2,model$shape)

   nam <- NULL
   for(i in 1:length(model$mu)) nam <- c(nam,paste("group",i))
   dimnames(estimates)[[2]] <- nam
   if(length(model$im.sdev) != 0 ) dimnames(estimates)[[1]] <- c("mu","sd","sigma2","sd","shape","sd")
   if(length(model$im.sdev) == 0 ) dimnames(estimates)[[1]] <- c("mu","sigma2","shape")
}
print(round(estimates,digits),...)
if(length(model$aic) != 0){
   	cat("\n AIC: ",model$aic,"\n")
   	cat("\n BIC: ",model$bic,"\n")
   	cat("\n EDC: ",model$edc,"\n")
   	cat("\n ICL: ",model$icl,"\n")
} 
   	cat("\n EM iterations: ",model$iter, "\n\n")
    invisible(model)
}
