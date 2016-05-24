#This function turns Cluster[] from Java into a list of more suitable S3 Objects.
r_clusters_from_java_clusters <- function(clus) {
  
  if(rJava::is.jnull(clus)) {
    warning("An error occured in the clustering function. Therefore, NULL is returned.")
    return(NULL)
  }
  
  subspace_matrix <- rJava::.jcall("ClusteringApplier",returnSig="[[Z",method="extract_subspace",clus,simplify=T)
  objects_matrix <- rJava::.jcall("ClusteringApplier",returnSig="[[I",method="extract_objects",clus,simplify=T)
  
  if(nrow(subspace_matrix)==0){
    warning("No subspace Clusters were generated. NULL is being returned. This is probably due to the parameters given to the clustering algorithm. Try a set of parameters that is more likely to produce many clusters")
    return(NULL)
  }
  
  res <- lapply(1:nrow(subspace_matrix),function(index){
    objects <- as.vector(objects_matrix[index,])
    #Add 1 because Java uses 0 as first index but R uses 1
    objects <- objects+1
    #Filter out those indices that were added in "extract_objects" to make the objects matrix rectangular
    objects <- objects[objects>0]
    return(subspace_cluster(subspace=as.vector(subspace_matrix[index,]),
                     objects=objects))
  })
  class(res) <- append(class(res),"subspace_clustering")
  return(res)
}
