# clustering
# experiment: dataset and metadata list
# method: "hc" (hierarchical clustering) or "kmeans"
# distance: "euclidean", "manhattan", spearman" or "pearson" methods
# clustMethod: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
# hc.type: "samples" or "variables"
# num.clusters: number of clusters in kmeans
clustering = function(dataset, method = "hc", distance = "euclidean",  
                      type = "samples", num.clusters = 5, clustMethod = "complete"){
	if (method == "hc") {
		result = hierarchical_clustering(dataset, distance, hc.type = type, 
                                     clustMethod = clustMethod)
	} 
  else if (method == "kmeans") {
		result = kmeans_clustering(dataset, num.clusters, type)
	}
	result
}


# HIERARCHICAL CLUSTERING
hierarchical_clustering <- function(dataset, distance='euclidean', clustMethod='complete', 
                                    hc.type = "samples")
{
	if (hc.type == "samples"){
		datamat = t(dataset$data)
	}
  else datamat = dataset$data
  
	if (distance == 'pearson' || distance == 'spearman'){
		dist.matrix = dist(1-cor(t(datamat), method = distance))
	} else {
		dist.matrix = dist(datamat, method = distance);
	}

	hc.tree = hclust(dist.matrix, method = clustMethod);
	hc.tree
}

# K-MEANS
kmeans_clustering <- function(dataset, num.clusters, type = "samples")
{	
  if (type == "samples"){
    datamat = t(dataset$data)
  }
  else datamat = dataset$data
  
	result.kmeans = kmeans(datamat, centers = num.clusters, nstart = 100, iter.max = 50)
	result.kmeans
}

kmeans_result_df = function(kmeans.result){
	kmeans.result.df = NULL
	for (i in 1:max(kmeans.result$cluster)){
		kmeans.result.df = rbind(kmeans.result.df, data.frame(cluster = i, samples = paste(names(kmeans.result$cluster[kmeans.result$cluster==i]),collapse=' ',sep=" ")))
	}
	kmeans.result.df = data.frame(kmeans.result.df)
	kmeans.result.df$samples = as.character(kmeans.result.df$samples)
	kmeans.result.df
}



######################### CLUSTERING PLOTS #########################


#kmeans plot
kmeans_plot = function(dataset, kmeans.result){
  num.clusters = max(kmeans.result$cluster)
  par(mfrow=c(num.clusters, 1), mar = c(2,2,2,2))
  for (i in 1:num.clusters){
    matplot(as.matrix(dataset$data[,kmeans.result$cluster == i]), type="l", col="grey", 
            main = paste("Cluster ", i, ", n = ", kmeans.result$size[i], sep =""), axes = F)
    lines(apply(as.matrix(dataset$data[,kmeans.result$cluster == i,drop=F]), 1, median), type="l", col="blue",lwd=1)
    axis(2)
    axis(1, 1:nrow(dataset$data), rownames(dataset$data))
  }
}


#dendrogram
"dendrogram_plot" = function(dataset, hc.result, column.metadata = 1, labels = NULL, ...){
  if (!is.null(labels)){
    labels.hc = labels
  } else {
	if (is.null(column.metadata)){
		labels.hc = colnames(dataset$data)
	} else {
		labels.hc = dataset$metadata[,column.metadata]
	}
  }
  hc.result$labels = labels.hc
  ggdendro::ggdendrogram(hc.result, ...)
}

"dendrogram_plot_col" = function(dataset, hc.result, classes.col, title = "", lab.cex = 1.0, leg.pos = "topright", ...) 
{
  classes = dataset$metadata[,classes.col]
  cluster = as.dendrogram(hc.result)
  cluster = dendrapply(cluster, color_leaf, dataset, classes, lab.cex)
  plot(cluster, main = title, horiz = FALSE, ...)
  leg.txt = levels(classes)
  leg.col = 1:length(levels(classes))
  leg.txt = c("Key", leg.txt)
  leg.col = c("black", leg.col)
  if (leg.pos != "none")
    legend(leg.pos, leg.txt, text.col = leg.col, bty = "n")
}

"color_leaf" = function (n, dataset, classes, lab.cex = 1.0) {
    if (is.leaf(n)) {
      a <- attributes(n)
      i <- match(a$label, get_sample_names(dataset))
      attr(n, "nodePar") <- c(a$nodePar, list(lab.col = as.integer(classes[i]), 
                                              lab.cex = lab.cex, pch = NA))
    }
    n
}
