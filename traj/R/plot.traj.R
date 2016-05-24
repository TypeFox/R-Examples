plot.traj <-
function(x, num.samples = 10, clust.num = NULL, color.vect = NULL, ...)
{

  clust = x$clusters
  data = as.matrix(x$data)
  time = as.matrix(x$time)
  
  
  if(!is.null(color.vect))
  {
    if(length(color.vect) != num.samples)
      stop("Need to specify as many colors as samples to plot.")
  }
  
  
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
  else 
  {
    if(clust.num < 1) stop("clust.num must be larger than 0.")
    unique.clusters = clust.num
  }
  
  if(is.null(color.vect))
    color.vect = sample(colors()[-10],num.samples)
  
  

  min.y = min(data, na.rm=T)
  max.y = max(data, na.rm=T)

  min.x = min(time,na.rm=T )
  max.x = max(time, na.rm=T)
  
  xlab = "time"
  ylab = "y"
  args = list(...)
  params = names(args)
  
  if(!("xlab" %in% params))
    args[["xlab"]] = xlab
  if(!("ylab" %in% params))
    args[["ylab"]] = ylab
                  
  # Loops through all clusters
  for(i_clust in unique.clusters)
  {
    clust.data.pos = which(clust$cluster == i_clust)
    
    
    if(length(clust.data.pos) < num.samples){
      sampled.data.pos = clust.data.pos
    }
    else
      sampled.data.pos = sample(clust.data.pos, num.samples)

    # Setup args for plotting
    args[["x"]] = time[sampled.data.pos[1],]
    args[["y"]] = data[sampled.data.pos[1],]
    if(!("xlim" %in% params))
      args[["xlim"]] = c(min.x, max.x)
    if(!("ylim" %in% params))  
      args[["ylim"]] = c(min.y, max.y)
    if(!("main" %in% params))
      args[["main"]] = paste("Cluster ", i_clust, sep = "")
    args[["col"]] = color.vect[1]

    
    # Generates the cluster plot with the first sampled data.
    do.call(plot, args)


    lines(time[sampled.data.pos[1],], data[sampled.data.pos[1],], col = color.vect[1])
    
    color.pos = 2
    # Loops throuth the rest of the samples in cluster and plots them on the same plot.
    for(i_id in sampled.data.pos[-1])
    {
      points(time[i_id,], data[i_id,], col = color.vect[color.pos])
      lines(time[i_id,], data[i_id,], col = color.vect[color.pos])
      color.pos = color.pos + 1
    }
  }
  # Add main title to the multiplot image. 
  if(is.null(clust.num))
  {
    if(num.samples > 1)
      text = paste("Cluster plots of data vs. time of", num.samples, " samples")
    else
      text = paste("Cluster plots of data vs. time of", num.samples, " sample")
    
    mtext(text, side = 3, line = 0.5, outer = TRUE, cex = 1.2) 
  } 
  
  if(is.null(clust.num)){
    par(mfrow = c(1,1))
  }
  
}
