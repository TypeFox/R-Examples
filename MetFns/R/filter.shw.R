filter.shw<-function(data, shw)
{
  if(!is.data.frame(data) || !is.character(shw) || nchar(shw)!=3 || !grepl("^[A-Z]+$",shw) )
      stop("invalid input parameter(s) specification: check data/shw ")
      
  if(!is.na(match("zero",colnames(data)))){          
     data.shw<-data[data$Shw==shw,]
 }
  else{
       k<-match("SPO",names(data))
       data.shw<-as.data.frame(setNames(replicate(20,numeric(0)), names(data)[1:(k+2)]))
       for(i in seq((k+1), ncol(data)-1,by=2)) {
              data.n<-data[which(data[,i]==shw),c(1:k,i:(i+1))]
              names(data.n)[(k+1):(k+2)]<-c("Shw","N")
              data.shw<-rbind(data.shw,data.n)
      }
 }
 data.shw
}
