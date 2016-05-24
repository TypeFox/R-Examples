#reverse coding
revcode<-function(data,vars){
x<-length(vars)
 for (i in 1:x){
 mx<-max(data[,vars[i]])
 mn<-min(data[,vars[i]])
 data[,vars[i]]<-(mx)-data[,vars[i]]+(mn)
}
return(data)
}
