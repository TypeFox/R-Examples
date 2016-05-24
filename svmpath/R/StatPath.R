StatPath<-function(f,y,Elbow){
    yhat<-y*f
  error<-sum(yhat<0)
  margin<-sum((1-yhat)[yhat<1])
  selbow<-length(Elbow)
list(error=error,margin=margin,selbow=selbow)
  }
