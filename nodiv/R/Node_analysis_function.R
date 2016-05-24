

# #INTERNAL FUNCTION  

 
nodiv_anal <- function(node, nodiv_data, repeats, method)
{
  # return a vector of NAs if the node is the basal node (as that has no parent node)
  if (node == basal_node(nodiv_data))  
    # forloop res_list <- append(res_list, NA) else {
    return(data.frame(ses = rep(NA, Nsites(nodiv_data)),  pval = rep(NA, Nsites(nodiv_data)), nodeemp = rep(NA, Nsites(nodiv_data)), nodemeans = rep(NA, Nsites(nodiv_data)), nodesds = rep(NA, Nsites(nodiv_data))))
  
  parent_data <- subsample.distrib_data(nodiv_data, species = Node_species(nodiv_data, Parent(node, nodiv_data)))
  
  parNode_sites <- match(parent_data$coords$sites, nodiv_data$coords$sites)
  # a vector indicating which of all sites are considered for this node
  
  # a boolean vector indicating which of species descending from the parent node that descend from the focal node
  Node_sp <- Node_species(nodiv_data, node)  
  
  res_object <- Nodesig(parent_data, Node_sp, repeats, method, show = F)

  res_object <- lapply(res_object, function(x) #ensure that dimensions are the same for merging into matrices
  {
    ret <- rep(NA, Nsites(nodiv_data))
    ret[parNode_sites] <- x
    ret
  })
  
  return(as.data.frame(res_object))
}


## EXPORTED FUNCTION

Node_analysis <- function(nodiv_data, repeats = 100, method = c("rdtable", "quasiswap"), cores = 1, log_parallel_progress = FALSE)
{
  if(!inherits(nodiv_data, "nodiv_data")) stop("This function only works on objects of class nodiv_data")
  method <- match.arg(method)
 
  paral <- FALSE
  if(cores > 1)
    if(requireNamespace("parallel"))
      paral <- TRUE
  if(paral)
  {
   cl <- parallel::makeCluster(cores)  
   pb <- txtProgressBar(min = Nspecies(nodiv_data) + 1, max = Ntip(nodiv_data$phylo) + Nnode(nodiv_data$phylo), style = 3)
  
   results <- parallel::parLapply(cl, nodenumbers(nodiv_data), function(node) 
   {
     #setTxtProgressBar(pb, node)
     ret <- nodiv_anal(node, nodiv_data, repeats, method)
     if(log_parallel_progress)
       save(ret, file = paste(node,".rda", sep = ""))
     ret[,1:2]
   })
  } else { 
    pb <- txtProgressBar(min = Nspecies(nodiv_data) + 1, max = Ntip(nodiv_data$phylo) + Nnode(nodiv_data$phylo), style = 3)
    results <- lapply(nodenumbers(nodiv_data), function(node)
    {
      setTxtProgressBar(pb, node)
      ret <- nodiv_anal(node, nodiv_data, repeats, method)
      ret[,1:2] #to keep object size down, I am currently only using these two variables for the big analysis
    })   
  } 
  nodiv_result <- nodiv_res(results, nodiv_data, repeats, method)
  return(nodiv_result)
}
  


