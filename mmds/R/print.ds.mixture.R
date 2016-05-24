print.ds.mixture<-function(x,...){
   # default plotting method for ds.mixture objects
   # stops all of the object being printed
   cat("\nds.mixture object \n")

   ttype<-"line"
   if(x$pt) ttype<-"point"

   cat(x$mix.terms,ttype,"transect model\n")
   cat(length(x$data$distance),"observations\n")
   cat("Truncation at",x$width,"\n")
   cat("AIC =", x$aic, "\n")

   invisible()
}
