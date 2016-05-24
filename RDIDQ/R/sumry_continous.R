sumry_continous <-
function(cont_var_test){
  t=cont_var_test
  non_missing<-colSums(!is.na(t))
  missing_values <-colSums(is.na(t))
  
  #Define Matrix
  n=length(cont_var_test)
  z=names(t)
  
  x1=matrix(0,nrow =n , ncol = 15) #final Matrix 
  # Variable name
  for (i in 1:n){
    x1[i,1]=z[i]
  }
  colnames(x1)<-c("Variable name","non-missing","missingvalues","min","mean","max","5%" ,"10%", "25%",  "50%",  "75%",  "90%",	"95%",	"99%",	"100%")
  x=matrix(0,nrow =n , ncol = 9)  #Quantile matrix
  rownames(x)<-z
  colnames(x)<-c("5%" ,"10%", "25%",  "50%",  "75%",	"90%",	"95%",	"99%",	"100%")
  
  
  #Final DATASHEET of DIDQ
  
  final=data.frame(x1)
  names(final)<-c("Variable name","non-missing","missingvalues","min","mean","max","5%" ,"10%", "25%",  "50%",	"75%",	"90%",	"95%",	"99%",	"100%")
  names(x)
  x1[,2]=unlist(non_missing)
  x1[,3]=unlist(missing_values)
  class(final)
  final[,2]=unlist(non_missing)
  final[,3]=unlist(missing_values)
  
  #Calculate Minimum value
  for(i in 1:ncol(t)){
    x1[i,4]=min(t[,i],na.rm=TRUE)
  }
  
  #Calculate Mean value
  for(i in 1:ncol(t)){
    x1[i,5]=mean(t[,i],na.rm=TRUE)
  }
  
  #Calculate Maximum value
  for(i in 1:ncol(t)){
    x1[i,6]=max(t[,i],na.rm=TRUE)
  }
  
  #Calculate Percentile
  
  
  for (i in 1:ncol(t)){
    x[i,]=quantile(t[,i], c(.05, .10, .25, .50, .75, .90, .95, .99,1),na.rm=TRUE)
  }
  
  for (i in 1:nrow(x)){
    l=6
    for (j in 1:ncol(x)){
      k=l+j
      x1[i,k]=x[i,j]
    }
  }
  final_cont_data=as.data.frame(x1)
  return(final_cont_data)
}
