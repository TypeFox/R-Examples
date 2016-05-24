PrintPath<-function(trace,step,obs,moveto,movefrom,lambda,digits,stats){
  if(trace){
    moveto<-rep(moveto,length=length(obs))
    movefrom<-rep(movefrom,length=length(obs))
  for(i in seq(along=obs)){
    cat(step,":\tObs ",obs[i],"\t",movefrom[i],"->",moveto[i],"  lambda = ",format(lambda, digits=digits,nsmall=digits),"  Sum Eps = ",format(round(stats$margin,2))," Elbow = ",stats$selbow," Error = ",stats$error,"\n",sep="")
  }
}
  invisible()
}
