#' Plots profiles of the groups.
#' 
#' @param X list containing the data matrices of all groups
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X identifying the measured data
#' @param xlab label of the x-axis of the plot
#' @param ylab label of the y-axis of the plot
#' @return Plots profiles of the groups.
#' @example R/example_plot.txt
#' @export
hrm.plot = function(X, group , factor1, subject, data, xlab="dimension", ylab="means" ){

  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1),is.character(data),is.character(xlab),is.character(ylab))
  #stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  
  group = as.character(group)
  factor1 = as.character(factor1)
  subject = as.character(subject)
  xlab=as.character(xlab)
  ylab=as.character(ylab)
  X = split(X, X[,group], drop=TRUE)
  a = length(X)
  d = nlevels(X[[1]][,factor1])
  n = rep(0,a) 
  
  means = data.frame(dimension=1:d)
  groupnames = c()
  
  
  for(i in 1:a){
    groupnames[i]=as.character(X[[i]][1,group])
    X[[i]] = X[[i]][ order(X[[i]][,subject], X[[i]][,factor1]), ]
    X[[i]]=X[[i]][,data]
    X[[i]] = matrix(X[[i]],ncol=d,byrow=TRUE)
    n[i] = dim(X[[i]])[1]
    means[,(i+1)]=colMeans(X[[i]])
  }
  
  colnames(means)=c("dimension",groupnames)
  
  means=melt(means, id.vars="dimension")
  colnames(means)=c("dimension", "group", "value")

  ggplot() +
     geom_line(data=means, aes(x=means$dimension, y=means$value,group=means$group,colour=means$group)) +
     geom_point(data=means, aes(x=means$dimension, y=means$value,group=means$group,colour=means$group),size=1.5) +
    xlab(xlab) + ylab(ylab)

}

# hrm.plot end ------------------------------------------------------------
