print.svmpath<-function(x,digits=6,maxsteps=length(x$Step),...){
  for( i in seq(maxsteps)){
    step<-x$Step[i]
    stats<-list(selbow=x$Size.Elbow[step],error=x$Error[step],margin=x$SumEps[step])
PrintPath(TRUE,x$Step[i],x$Obs.step[i],x$Moveto[i],x$Movefrom[i],x$lambda[step],digits,stats)
  }
  invisible()
}
  
