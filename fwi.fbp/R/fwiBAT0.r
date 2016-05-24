.fwiBAT0<-function(input,init=c(ffmc_yda=85,dmc_yda=6,dc_yda=15,lat=55),out="all",lat.adjust=TRUE){
  if(!is.na(charmatch("input",search()))) {detach(input)} 
  out.fwi<-NULL
  for (j in 1:nrow(input)){
        if (j==1){out.fwi0<-fwi(input[j,],init=init,out=out,lat.adjust=lat.adjust)} else {
           out.fwi0<-fwi(input=input[j,],yda.fwi=out.fwi0,init=init,out=out,lat.adjust=lat.adjust)}
           out.fwi<-rbind(out.fwi,out.fwi0)
        }
  out.fwi
}
