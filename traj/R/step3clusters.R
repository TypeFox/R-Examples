step3clusters <-
function(trajFactors, nclusters = NULL, nstart = 50, criteria = "ccc", forced.factors = NULL)
{

  if(is.null(forced.factors)){
    data = trajFactors$factors
    
  }
  else{
    data = trajFactors$measurments[,c("ID", forced.factors)]
    colnames(data)[1] = "output"
  }
  
  # Error checking
  if(class(data) != "data.frame")
    stop("data muts be a data.frame")
  
  if(nclusters > nrow(data) && !is.null(nclusters))
    stop("Requesting more clusters in 'nclusters' than available rows in data.")
  
  #Sizing data
  dim.of.data = dim(data)
  sample.size = dim.of.data[1]
  
  
  # Deal with IDs
    IDvector = data[1] 
    data = data[-1]
  
  max.num.obs  = dim(data)[2]
  
  cluster.est = NULL

  # Calculate the number of clusters to use
  if(is.null(nclusters)){
    cluster.est = NbClust(data, method = "kmeans", index = criteria)
    
    all.criteria = as.data.frame(cluster.est$All.index)
    all.criteria.x = as.integer(rownames(all.criteria))

    par(mfrow = c(1,2))

    plot(all.criteria.x, as.matrix(all.criteria),
         main = paste(criteria, " criteria " , "versus Clusters"),
         xlab = "Clusters", 
         ylab = "Criteria")
    
    num.clust = cluster.est$Best.nc[1]
          
    wss <- (nrow(data)-1)*sum(apply(data,2,var))
    
    for (i in 2:15) wss[i] <- sum(kmeans(data, 
                                         centers=i)$withinss)
    plot(1:15, wss, type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares", main = "Scree Plot for Number of Clusters")
    
  }
  else
  {
    if(nclusters < 1 ) stop("forcec.clust must be larger than 0")
    
    num.clust = round(nclusters)
  }
  
  # use k-means to split the data into the designated number of clusters
  cluster.data = kmeans(data, centers = num.clust, nstart=nstart)
  
  # Bind the cluster position to ID vector
  output = cbind(IDvector, cluster.data$cluster)
  output = as.data.frame(output)
  names(output) = c("ID", "cluster")
  
  table.output = rbind(table(output$cluster) , table(output$cluster) / sum(table(output$cluster)) * 100)
  rownames(table.output) = c("(#)", "(%)")
  
  data = cbind(IDvector, data)
  names(data)[1] = "output"
  
  # Create object "traj" to export
  structure(list(clusters = output, clust.distr = table(output$cluster), 
                 clust.estim = as.data.frame(cluster.est$All.criteria),  
                 factors  = data, 
                 e.values = trajFactors$e.values, princ.fact = trajFactors$princ.fact,
                 measurments = trajFactors$measurments, data = trajFactors$data, 
                 time = trajFactors$time ),class='traj')
  
}
