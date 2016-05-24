
odds<-function(x, model=c("log-additive", "dominant", "recessive", "overdominant", "codominant"), sorted=c("no","p-value","or")){
   tables<-attr(x,"tables")

   if(is.null(tables)) stop("An object fitted with WGassociation is needed. scanWGassociation doesn't estimate ORs.")

   models<-c("log-additive", "additive" , "dominant", "recessive", "overdominant", "codominant")
   mo<-pmatch(tolower(model[1]),models)
   if (is.na(mo)) stop("Incorrect model chosen")

   add<-lapply(tables, function(o){
                       if(is.null(dim(o))){
                          c(NA,NA,NA,NA)
                       } else {
                          tag<- c("0,1,2", "0,1,2", "Dominant", "Recessive", "Overdominant","Codominant")
                          lag<- c(0,0,2,2,2,2)
                          r<-match(tag[mo],rownames(o))+lag[mo]
                          row<-o[r,5:8]
                          if(mo>2) row[4]<-o[r-1,8]
                          if (mo==6) {
                              lab<-names(row)
                              row<-c(row[1:3],o[r+1,5:7],row[4])
                              dim(row)<-c(7,1)
                              rownames(row)<-c(paste(lab[1:3],rep("1",3),sep="."),paste(lab[1:3],rep("2",3),sep="."),lab[4])
                          }
                          row
                       }
                  }
        )
   add<-as.data.frame(add)
   add<-t(add)

   pvals<-attr(x,"pvalues")

   codo<- match("codominant", colnames(pvals))
   
   if(!is.na(codo)){
      pv<-match("p-value", colnames(add))
      add[is.na(add[,pv]),pv]<-pvals[is.na(add[,pv]),codo]
   }

   so<-pmatch(tolower(sorted[1]),c("p-value","or"))
   if(!is.na(so)){
      if(so==1) col<-match("p-value", colnames(add))
      if(so==2 & mo< 6) col<-1
      if(so==2 & mo==6) col<-4
      add<-add[order(add[,col]),]
   }
   colnames(add)[length(colnames(add))]<-paste("p-value",models[mo],sep=".")
   as.data.frame(add)
}
