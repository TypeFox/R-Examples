Sub.filename<-function(filename){
   if(grepl(".xls",filename)){
      name<-sub(".xls","",filename)
   }
   
   if(grepl(".xlsx",filename)){
      name<-sub(".xlsx","",filename)
   }
   
    if(grepl(".csv",filename)){
      name<-sub(".csv","",filename)
   }
   
   return(name)
}