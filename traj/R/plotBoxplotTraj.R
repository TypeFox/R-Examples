plotBoxplotTraj <-
function(x, clust.num = NULL,  ...)
{
  clust = x$clusters
  data = as.matrix(x$data)
  time = as.matrix(x$time)
  
  data = data[,-1]
  time = time[,-1]

  unique.clusters = sort(unique(clust[,2]))
  
  num.clust = length(unique.clusters)
  
  
  # Create multiplot canvas
  if(is.null(clust.num))
  {
    num.plot.x = round(sqrt(num.clust))
    if(num.plot.x^2 < num.clust)
      num.plot.y = num.plot.x + 1
    else
      num.plot.y = num.plot.x
    
    par(mfrow=c(num.plot.y, num.plot.x), oma = c(0,0,4,0) )
    
    
  }
  else unique.clusters = clust.num
  
  min.y = min(data, na.rm=T)
  max.y = max(data, na.rm=T)
  
  min.x = min(time, na.rm=T)
  max.x = max(time, na.rm=T)
  
  xlab = "time"
  ylab = "y"
  args = list(...)
  params = names(args)
  
  
  # Loops through all clusters
  for(i_clust in unique.clusters)
  {
    
    clust.data.pos = which(clust$cluster == i_clust)
    
    clust.data = data[clust.data.pos,]
    
    clust.time = time[clust.data.pos,]


    unique.clust.time = as.numeric(names(table(clust.time)))


    #if(!("xlim" %in% params))
    #  args[["xlim"]] = c(min.x, max.x)
    if(!("ylim" %in% params))
      args[["ylim"]] = c(min.y, max.y)
    if(!("main" %in% params))
      args[["main"]] = paste("Cluster ", i_clust, sep = "")
    if(!("xlab" %in% params))
      args[["xlab"]] = xlab
    if(!("ylab" %in% params))
      args[["ylab"]] = ylab
    
    clust.data = as.vector(clust.data)
    clust.time = as.vector(clust.time)
    
    clust.bind = cbind(clust.data, clust.time)
    
    args[[""]] = formula(clust.data ~ clust.time)
    args[["data"]] = clust.bind
    
    do.call(boxplot, args)    

  }
  
  if(is.null(clust.num))
  {
    text = "Boxplots for Every Cluster"
    mtext(text, side = 3, line = 0.5, outer = TRUE, cex = 1.2) 
  }
  
  if(is.null(clust.num)){
    par(mfrow = c(1,1))
  }
  
}
