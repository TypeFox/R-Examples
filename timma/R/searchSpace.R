#' Generate search space
#' 
#' A function to generate the search space for sffs
#' 
#' @param drug_number an integer to specify the number of drugs
#' @param k_set a vector to specify the selected target set
#' @param profile_data drug-target interaction data
#' @param y_actual the drug sensitivity data
#' 
#' @return a list of the following components:
#' \item{IM_d}{search space of identical sets}
#' \item{IM_superset}{search space of supersets}
#' \item{IM_subset}{search space of subsets}
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' num<-length(tyner_sensitivity[,1])
#' k_set<-rep(0, dim(tyner_interaction_binary)[2])
#' k_set[1]<-1
#' space<-searchSpace(num, k_set, tyner_interaction_binary, tyner_sensitivity[,1])

searchSpace <- function(drug_number, k_set, profile_data, y_actual) {
    # parameter 1: drug_num, the number of drugs 
    # parameter 2: k_set, the current kinase set 
    # parameter 3: profile_data, the complete drug-target-profile data
    # parameter 4: y_actual, the actual efficacy
    
    # extend one col by 0
    col0 <- rep(0, drug_number)
    prof0 <- cbind(profile_data[, which(k_set == 1)], col0)
    
    
    # extend one col by 1
    col1 <- rep(1, drug_number)
    prof1 <- cbind(profile_data[, which(k_set == 1)], col1)
    prof <- unique(rbind(prof1, prof0))
    dec_prof <- apply(prof, 1, function(x) strtoi(paste(x, collapse = ""), base = 2))
    
    
    dec <- apply(prof0, 1, function(x) strtoi(paste(x, collapse = ""), base = 2))
    # for identical
    col_num <- length(dec_prof)
    identical_idx <- sapply(dec, function(x) which(dec_prof == x))
    
    
    
    IM_d <- array(NA, dim = c(drug_number, col_num, 2))
    IM_subset <- array(Inf, dim = c(drug_number, col_num, 2))
    IM_superset <- array(-Inf, dim = c(drug_number, col_num, 2))
    
    for (i in 1:drug_number) {
        # get the decimal
        IM_d[i, identical_idx[i], 1] <- y_actual[i]
        
        
        # get the binary set: superset and subset
        bin_set <- binarySet(prof0[i, ])
        # ismember function R version: match
        subset_index <- dec_prof %in% bin_set$subset
        IM_subset[i, subset_index, 1] <- y_actual[i]
        
        superset_index <- dec_prof %in% bin_set$superset
        IM_superset[i, superset_index, 1] <- y_actual[i]
        
    }
    
    
    dec <- apply(prof1, 1, function(x) strtoi(paste(x, collapse = ""), base = 2))
    # for identical
    col_num <- length(dec_prof)
    identical_idx <- sapply(dec, function(x) which(dec_prof == x))
    
    for (i in 1:drug_number) {
        # get the decimal
        IM_d[i, identical_idx[i], 2] <- y_actual[i]
        
        
        # get the binary set: superset and subset
        bin_set <- binarySet(prof1[i, ])
        # ismember function R version: match
        subset_index <- dec_prof %in% bin_set$subset
        IM_subset[i, subset_index, 2] <- y_actual[i]
        
        superset_index <- dec_prof %in% bin_set$superset
        IM_superset[i, superset_index, 2] <- y_actual[i]
        
    }
    
    return(list(d = IM_d, i = IM_superset, o = IM_subset))
} 
