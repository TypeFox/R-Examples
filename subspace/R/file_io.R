
cluster_from_string <- function(str,dim) {
  values <- str %>% strsplit(split=" ") %>% unlist() %>% as.numeric()
  subspace <- values[1:dim]==1
  objects <- values[(dim+1):length(values)]
  return(subspace_cluster(subspace=subspace,objects=objects))
}
cluster_to_string <- function(cluster) {
  subspacevector <- sapply(cluster$subspace,function(bool){if(bool){1}else{0}})
  res <- c(subspacevector,cluster$objects) %>% as.character() %>% paste(collapse=" ")
  return(res)
}

#'Reads a Clustering Object from a File
#'
#'Reads a File and creates an object of class \emph{subspace_clustering} from
#'that.
#'
#'@param file_path The path to the File from which the clustering should be read
#'@param index_starts_at The index used in the file to refer to the first object
#'  in The data matrix
#'@param dim if provided, overrides any value for dim that is found in the first
#'  line of the file
#'  
#'@note Files must have the following Format: The first line contains the
#'  substring "DIM=*dim*;" where *dim* is the number of dimensions of the data
#'  set.
#'  
#'  Each subsequent line corresponds to a cluster and contains only numbers
#'  separated by spaces. The first *dim* of these values have to be either '0'
#'  or '1' and indicate in which subspace a cluster exists. All other values in
#'  the line have to be the row numbers of the objects that the cluster
#'  contains. Row numbers in the file are assumed to be 0-indexed and are
#'  changed to 1-indexed as they are loaded into R. This can be changed with the
#'  parameter \emph{index_starts_at}. E.g. a clustering for a three-dimensional
#'  dataset with one cluster that is in the first and third dimension and 
#'  contains the first, second and 1337-th object has to be represented as:
#'  
#'  DIM=3;\cr 1 0 1 0 1 1336
#'@seealso \code{\link{clustering_to_file}}
#'@export
clustering_from_file <- function(file_path,index_starts_at=0,dim=NULL) {
  #read the file into memory. We should not have memory size problems with clustering files, so 
  #just read it line-wise as of big vector of strings
  handle <- file(file_path)
  file_content_array <- readLines(handle)
  close(handle)
  
  if(is.null(dim)){
    #Now we get the number of dimensions "dim" of the clustering. This is necessary because each line just contains
    #space-separated numbers, the first "dim" of which represent the dimensions in which a cluster is a cluster.
    dim_first_index <- stringr::str_locate_all(file_content_array[1],pattern="DIM=")[[1]][,"end"]+1
    dim <- file_content_array[1] %>% 
      substr(start=dim_first_index,stop=stringr::str_length(file_content_array[1])) %>%
      strsplit(split=";") %>% 
      unlist() %>%
      utils::head(n=1) %>%
      as.numeric()
  }
  #Now we have determined dim and turn all the other lines into a subspace_clustering
  clustering <- lapply(file_content_array[-1],cluster_from_string,dim=dim)
  #we also need to increase the index of objects to match the starting-index specified in index_starts_at.
  #In the default case this means that we increase all of the values by one because the indexes in 
  #R start with 1 and the indexes in the file start at 0.
  clustering <- lapply(clustering,function(cluster){
    cluster$objects <- cluster$objects+(1-index_starts_at)
    return(cluster)})
  class(clustering) <- append(class(clustering),"subspace_clustering")
  return(clustering)
}
#Turns a vector of assignments into a clustering. The vector should be as long
#as the data set has rows. Each element of the vector should be a positive
#integer value. Negative values are assumed to be noise. If two elements of the
#vector have the same values then that means that the corresponding objects are
#in the same cluster. E.g. the vector c(1,-1,5,6,1,2) means that the first and
#the fifth object in the data set form a cluster, the second is a noise object
#and the others are the only objects of their clusters. Note that because of the
#vector format, each object can be in exactly one cluster and no information
#about subspaces is encoded. The resulting clustering will therefore claim that
#every cluster is in the full dataspace. If the subspace information is
#important, the dimension of the data space should therefore also be passed to 
#this function.
clustering_from_vector <- function(vec,dim=5) {
  number_map <- integer(0)
  better_assignment_vector <- sapply(1:length(vec), function(i) {
    if(vec[i]<0 | is.na(vec[i])) {
      return(-1)
    } else if (vec[i] %in% number_map) {
      return(which(number_map==vec[i]))
    } else {
      #Note: <<- is global assign so that this function can have a side effect despite being used in an sapply
      number_map <<- append(number_map,vec[i])
      return(which(number_map==vec[i]))
    }
  })
  clustering <- lapply(1:max(better_assignment_vector),function(val){
    return(subspace_cluster(objects=which(better_assignment_vector==val),subspace=rep(T,dim)))
  })
  class(clustering) <- append(class(clustering),"subspace_clustering")
  return(clustering)
}

#'Write a Subspace Clustering to Disk
#'
#'@param clustering a subspace clustering object as generated by one of the
#'  functions from the \emph{subspace} package
#'@param file_path the path to the file into which the clustering should be
#'  written
#'@param index_should_start_at the value that is used to refer to the first
#'  value in the dataset.
#'@note By default, R uses the value \emph{1} when referring to the first object
#'  in a data frame or array, while most other languages use \emph{0}. To make
#'  working with this convention easy, clusterings written to disk are converted
#'  to this 0-indexing System. The standard parameter for the corresponding
#'  function \code{\link{clustering_from_file}} is set in such a way that files
#'  read will automatically be converted to 1-indexes, which means that you
#'  should never need to change this parameter if you work exclusively with the
#'  \emph{subspace} package.
#'@seealso \code{\link{clustering_from_file}}
#'@export
clustering_to_file <- function(clustering,file_path,index_should_start_at=0) {
  if(is.null(clustering)) {
    return()
  }
  dim <- length(clustering[[1]]$subspace)
  #make sure the indexes are correct. We assume that the clustering passed here is already 1-indexed, i.e. the first
  #row in the data table has index 1. Then, the indexes that are written to a
  #file are changed so that they start with the number specified in the index_should_start_at parameter.
  clustering <- lapply(clustering,function(cluster){
    cluster$objects <- cluster$objects -(1-index_should_start_at)
    return(cluster)})
  
  strings_to_write <- sapply(clustering,cluster_to_string)

  handle <- file(file_path,"w")
  cat(paste("DIM=",as.character(dim),";\n",sep=""),file=handle)
  cat(strings_to_write,file=handle,sep="\n")
  close(handle)
  
}


