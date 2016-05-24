plotCombTraj <-
function(x, stat.type = "mean", colored = FALSE, ...)
{
  
  clust = x$clusters
  data = as.matrix(x$data)
  time = as.matrix(x$time)
  
  
  data = data[,-1]
  time = time[,-1]
  
  unique.clusters = sort(unique(clust[,2]))
  
  num.clust = length(unique.clusters)
  
  min.y = min(data, na.rm=T)
  max.y = max(data, na.rm=T)
  
  min.x = min(time, na.rm=T)
  max.x = max(time, na.rm=T)
  
  
  max.y = max.y + num.clust * 0.04 * (max.y - min.y) 
  
  xlab = "time"
  ylab = "y"
  args = list(...)
  params = names(args)
  
  
  for(i_clust in unique.clusters)
  {
    
    clust.data.pos = which(clust$cluster == i_clust)
    
    clust.data = data[clust.data.pos,]
    
    clust.time = time[clust.data.pos,]
    
    unique.clust.time = as.numeric(names(table(clust.time)))
    
    comp.stat = data.frame(matrix(ncol = 3, nrow = 1))
    
    for(i_time in unique.clust.time)
    {
      time.pos = which(clust.time == i_time)
      
      data.time = clust.data[time.pos]
      
      if(stat.type == "mean"){
        stat.data = mean(data.time)
        centiles = 0
        
        comp.stat = rbind(comp.stat, c(stat.data, centiles[1], centiles[2]))
      }
      else if(stat.type == "median"){
        stat.data = median(data.time)
        centiles = 0
        
        comp.stat = rbind(comp.stat, c(stat.data, centiles[1], centiles[2]))
      }
    }
    
    comp.stat = comp.stat[-1,]
    
    stat.name = paste(toupper(substr(stat.type, 1, 1)), substr(stat.type, 2, nchar(stat.type)), sep = "")
    
    # Setup args for plotting
    args[["x"]] = unique.clust.time
    args[["y"]] = comp.stat[,1]
    if(!("xlim" %in% params))
      args[["xlim"]] = c(min.x, max.x)
    if(!("ylim" %in% params))
      args[["ylim"]] = c(min.y, max.y + num.clust * 0.05 * (max.y - min.y) + 2 * 0.05 * (max.y - min.y))
    if(!("main" %in% params))
      args[["main"]] = paste(stat.name," Trajectory of all Clusters", sep = "")
    if(!("xlab" %in% params))
      args[["xlab"]] = xlab
    if(!("ylab" %in% params))
      args[["ylab"]] = ylab
    if(colored)
      args[["col"]] = i_clust
    
      args[["pch"]] = i_clust
    # Generates the cluster plot with the first sampled data.
    if(i_clust == 1)
      do.call(plot, args)
    else if(colored)
      points(unique.clust.time, comp.stat[,1], col = i_clust, pch = i_clust)
    else
      points(unique.clust.time, comp.stat[,1], pch = i_clust)
    
    if(colored)
      lines(unique.clust.time, comp.stat[,1], col = i_clust)
    else
      lines(unique.clust.time, comp.stat[,1], lty = i_clust)
  }
  
  if(colored){
    legend(min.x + 0.05 * (max.x - min.x), max.y + num.clust * 0.07 * (max.y - min.y),
           unique.clusters,
           lty=rep(1, num.clust),
           lwd= rep(2.5, num.clust),
           col=unique.clusters,
           pch = unique.clusters)
  }
  else{
    legend(min.x + 0.05 * (max.x - min.x), max.y + num.clust * 0.07 * (max.y - min.y),
           unique.clusters,
           lty=unique.clusters,
           lwd= rep(2.5, num.clust),
           col=rep(1, num.clust),
           pch = unique.clusters)
  }
}
