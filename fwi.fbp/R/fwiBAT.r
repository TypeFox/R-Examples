fwiBAT<-function(input,init=c(ffmc_yda=85,dmc_yda=6,dc_yda=15,lat=55),out="all",lat.adjust=TRUE){
  .Deprecated(new="fwi", 
              package="cffdrs", 
              msg="The 'fwi.fbp' package and contained functions are being deprecated and replaced by the 'cffdrs' package, please update your code to use the 'cffdrs' package.",
              old="fwiBAT")
  out.fwi0<-out.fwi<-NULL
  n0<-round(nrow(input)/365)
  n<-ifelse(365*n0>=nrow(input),n0,n0+1)
  for (i in 1:n){
   input0<-input[(365*(i-1)+1):(365*i),]
   input0<-input0[!is.na(input0[,1]),]
   if (i==1){
        out.fwi0<-.fwiBAT0(input0,init=init,out=out,lat.adjust=lat.adjust)
     } else {
        out.fwi0<-.fwiBAT0(input0,init=out.fwi[nrow(out.fwi),c("ffmc","dmc","dc","lat")],out=out,lat.adjust=lat.adjust)
     }
   out.fwi<-rbind(out.fwi,out.fwi0)
  }
  out.fwi
}