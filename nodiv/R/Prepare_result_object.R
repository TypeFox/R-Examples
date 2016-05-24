######################################################
# Functions to interpret the representation analysis #
######################################################


## FUNCTIONS TO CREATE THE DATA OBJECT - all internal (not exported)



parent_representation = function(node_number, rep_matrix, nodiv_data)
{
  desc = Descendants(node_number, nodiv_data)
  if(desc[1] < basal_node(nodiv_data) | desc[2] < basal_node(nodiv_data))
    return(rep(NA, nrow(rep_matrix)))
  desc1row = desc[1] - Nspecies(nodiv_data)
  desc2row = desc[2] - Nspecies(nodiv_data)
  
  return(rowMeans(cbind(rep_matrix[,desc1row], 1-rep_matrix[, desc2row])))
}

nodenames <- function(nodiv_data)
{
  ret <- nodenumbers(nodiv_data)
  if(!is.null(nodiv_data$phylo$node.label))
    ret <- paste(ret, nodiv_data$phylo$node.label)
  ret
}

nodiv_res <- function(results, nodiv_data, repeats, method)
{
  ret <- nodiv_data
  class(ret) <- c("nodiv_result", class(nodiv_data))
  ret$method <- method
  ret$repeats <- repeats
  
  SR <- sapply(results, "[[", 1)  
  ret$SOS <- sapply(nodenumbers(nodiv_data), function(node) parent_representation(node, SR, nodiv_data))
  colnames(ret$SOS) <- nodenames(nodiv_data)
  rownames(ret$SOS) <- nodiv_data$coords$sites
  
  rval <- sapply(results, "[[", 2)
  par_rval <- sapply(nodenumbers(nodiv_data), function(node) parent_representation(node, rval, nodiv_data))
  pval <- apply(par_rval, 2, pval_sig)
  
  #making sure that none of the values are more extreme than merited by the number of repeats
  pval[pval > 1-2/repeats] <- 1-2/repeats
  pval[pval < 2/repeats] <- 2/repeats
  
  ret$GND <- apply(pval, 2, function(x) 1-inv_logit(mean(logit(x), na.rm = T)))
  
  return(ret)
}





# Summarize.sites <- function(dispersion)
#   #Summarizes the information in the representation matrix in a by-sites manner
# {
#   # dispersion 		: a representation matrix resulting from running measure_dispersion()
#   # dat.LL			: the geographical coordinates of the sites
#   
#   cell <- rownames(dispersion)  
#   nodeeff <- apply(dispersion, 1, function(x) mean(abs(x),na.rm = T))
#   mapdata <- data.frame(cell, nodeeff)
#   return(mapdata)
# }
# 
# parent_representation = function(node_number, rep_matrix, nodiv_data)
#   # takes the representation matrix, and summarizes at the parent node (because sister species in the representation matrix are mirror images)
# {
#   desc = Descendants(node_number, nodiv_data)
#   if(desc[1] < basal_node(nodiv_data) | desc[2] < basal_node(nodiv_data))
#     return(rep(NA, nrow(rep_matrix)))
#   desc1row = which(colnames(rep_matrix) == as.character(desc[1]))
#   desc2row = which(colnames(rep_matrix) == as.character(desc[2]))
#   return(rowMeans(cbind(rep_matrix[,desc1row], -rep_matrix[, desc2row])))
# }



# parent_pval_representation_matrix <- function(rep_matrix, nodiv_data)
# {
#   retmat <- sapply(nodenumbers(nodiv_data), parent_pval_representation, rep_matrix = rep_matrix, nodiv_data)
#   colnames(retmat) <- colnames(rep_matrix)
#   rownames(retmat) <- rownames(rep_matrix)
#   return(retmat)
# }
# 
# parent_representation_matrix <- function(rep_matrix, nodiv_data)
# {
#   retmat <- sapply(nodenumbers( nodiv_data), parent_representation, rep_matrix = rep_matrix, nodiv_data = nodiv_data)
#   colnames(retmat) <- colnames(rep_matrix)
#   rownames(retmat) <- rownames(rep_matrix)
#   return(retmat)
# }

# 
# create.datalist <- function(dispersion, sitestatistics, coords)
# {
#   # a function to create the final result object after the representation analysis has completed
#   datalist <- list()
#   dispersion[abs(dispersion) == Inf] <- NA
#   datalist$siteresults <- merge(sitestatistics, Summarize.sites(dispersion), by = "cell")
#   datalist$siteresults <- merge(datalist$siteresults, coords, by = "cell")
#   datalist$siteresults <- datalist$siteresults[match(as.character(coords$cell), datalist$siteresults$cell),]
#   datalist$noderesults <- Summarize.nodes(dispersion)
#   datalist$rep_matrix <- dispersion[match(as.character(datalist$siteresults$cell), rownames(dispersion)),]
#   datalist$parent_rep_matrix <- parent_representation_matrix(datalist$rep_matrix)
#   datalist$parent_rep_matrix[abs(datalist$parent_rep_matrix) == Inf] <- NA
#   pod <- apply(datalist$parent_rep_matrix,2, function(x) mean(abs(x),na.rm = T))  #I guess this really belongs in the summarize.nodes function
#   pod[pod == Inf] <- NA
#   datalist$noderesults$parent_overdisp <- pod
#   return(datalist)
# }


