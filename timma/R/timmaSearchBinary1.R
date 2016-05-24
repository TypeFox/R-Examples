#'Prediction in the search space with two.sided TIMMA model 
#'
#'A function to return the prediction error in the search space for sffs
#'
#'@param profile_k current selected drug-target interaction data
#'@param space the search space returned by \code{\link{searchSpace}} function
#'@param sens drug sensitivity data
#'@param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#'@return the prediction error
#'@author Liye He \email{liye.he@@helsinki.fi}
#' @examples 
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' num<-length(tyner_sensitivity[,1])
#' k_set<-rep(0, dim(tyner_interaction_binary)[2])
#' k_set[1]<-1
#' space<-searchSpace(num, k_set, tyner_interaction_binary, tyner_sensitivity[,1])
#' profile_k<-tyner_interaction_binary[, which(k_set==1)]
#' error<-timmaSearchBinary1(profile_k, space, tyner_sensitivity[,1])

timmaSearchBinary1 <- function(profile_k, space, sens, loo = TRUE) {
    # space is 3d array
    dim_info <- dim(space$d)
    rows <- dim_info[1]
    cols <- dim_info[2]
    
    
    # IM_d<-array(NA, dim=dim_info[1:3])
    IM_d <- array(NA, dim = dim_info[1:2])
    IM_superset <- array(-Inf, dim = dim_info[1:2])
    # IM_i<-arrayminfcpp(dim_info[1], dim_info[2], dim_info[3])
    IM_subset <- array(Inf, dim = dim_info[1:2])
    # IM_o<-arrayinfcpp(dim_info[1], dim_info[2], dim_info[3]) S_r<-rep(0,rows)
    identical_idx <- rep(0, rows)
    for (i in 1:rows) {
        index <- profile_k[i] + 1
        IM_d[i, ] <- space$d[i, , index]
        IM_superset[i, ] <- space$i[i, , index]
        IM_subset[i, ] <- space$o[i, , index]
        # S_r[i]<-which((!is.na(IM_d[i,]))==TRUE)
        identical_idx[i] <- which((!is.na(IM_d[i, ])) == TRUE)
    }
    
    # M_d<-apply(IM_d, MARGIN=c(1, 2), sum, na.rm=TRUE)/apply(IM_d,MARGIN=c(1,2), function(x){sum(!is.na(x))})
    M_d <- sumcpp1(IM_d, rows, cols)
    # M_d<-apply(IM_d,2, sum, na.rm=T)/apply(IM_d, 2, function(x) {sum(!is.na(x))})
    # min_subset<-apply(IM_subset, 2, min) min_index<-apply(IM_subset, 2, which.min)
    # max_superset<-apply(IM_superset, 2, max) max_index<-apply(IM_superset, 2, which.max)
    maxval <- maxcpp1(IM_superset, rows, cols)
    minval <- mincpp1(IM_subset, rows, cols)
    min_subset <- minval$min
    min_index <- minval$min_idx
    max_superset <- maxval$max
    max_index <- maxval$max_idx
    
    
    
    # find cell which needs maximization averaging
    cell <- is.nan(M_d) & is.finite(max_superset)  # is.nan or is.na????????
    cell <- which(cell == TRUE)
    if (length(cell) != 0) {
        for (i in cell) {
            # row<-((i-1) %% rows) + 1 col<-floor((i-1) / rows)+1 print(ii) the drug sets that are the subset of the
            # cell
            drug_sub_cell <- !is.infinite(IM_superset[, i])
            # the drug index which achieves max sensitivity
            index <- max_index[i]
            # the dec of the drug with max sensitivity
            dec_maxsens <- identical_idx[index]
            
            # find the supersets of S(index,:) in S that has smaller sensitivity
            supersets_small <- IM_subset[, dec_maxsens] < max_superset[i]
            
            # find common item with drug_sub_cell and supersets_small
            common_cell <- which(drug_sub_cell & supersets_small)
            
            if (length(common_cell) != 0) {
                # max averaging
                k <- 1
                for (j in common_cell) {
                  max_superset[i] <- (max_superset[i] * k + sens[j])/(k + 1)
                  k <- k + 1
                }
            }
        }
    }
    
    cell2 <- is.nan(M_d) & is.finite(min_subset)
    cell2 <- which(cell2 == TRUE)
    if (length(cell2) != 0) {
        for (i in cell2) {
            # row<-((i-1) %% rows) + 1 col<-floor((i-1) / rows)+1 the drug sets that are the superset of the cell
            drug_sub_cell <- !is.infinite(IM_subset[, i])
            # the drug index which achieves min sensitivity
            index <- min_index[i]
            
            # the dec of the drug with min sensitivity
            dec_minsens <- identical_idx[index]
            
            # find the subsets of S(index,:) in S that has higher sensitivity
            subsets_small <- IM_superset[, dec_minsens] > min_subset[i]
            # find common item with drug_sub_cell and supersets_small
            if (length(subsets_small) == 0) {
                common_cell2 <- vector("numeric")
            } else {
                common_cell2 <- which(drug_sub_cell & subsets_small)
            }
            if (length(common_cell2) != 0) {
                # min averaging
                k <- 1
                for (j in common_cell2) {
                  min_subset[i] <- (min_subset[i] * k + sens[j])/(k + 1)
                  k <- k + 1
                }
            }
        }
    }
    
    M <- M_d
    M[cell] <- (max_superset[cell] + 1)/2
    M[cell2] <- (min_subset[cell2] + 0)/2
    # cels that not only have lower boundery and also have upper boundary
    average_index <- intersect(cell, cell2)
    M[average_index] <- (max_superset[average_index] + min_subset[average_index])/2
    
    # predicted error
    error_predict <- rep(NA, rows)
    # predicted efficacy
    pred <- rep(NA, rows)
    if (loo == FALSE) {
        pred <- M[identical_idx]
        error_predict <- abs(pred - sens)
        
    } else {
        for (i in 1:rows) {
            # remove drug i, namely remove the i-th row
            
            # get the dim info dim_IMd<-dim(IM_d) dim_IMd[3]<-dim_IMd[3]-1
            dim_IMd <- c(rows - 1, cols)
            
            IM_d_loo <- array(IM_d[-i, ], dim = dim_IMd)
            IM_subset_loo <- array(IM_subset[-i, ], dim = dim_IMd)
            # IM_subset_loo<-newarray(IM_subset[,,-i], dim_IMd[1], dim_IMd[2], dim_IMd[3])
            IM_superset_loo <- array(IM_superset[-i, ], dim = dim_IMd)
            # IM_superset_loo<-newarray(IM_superset[,,-i], dim_IMd[1], dim_IMd[2], dim_IMd[3])
            sens_loo <- sens[-i]
            drug_idx_loo <- identical_idx[-i]
            
            # M_d_loo<-apply(IM_d_loo, 2, sum, na.rm=TRUE)/apply(IM_d_loo,2, function(x){sum(!is.na(x))})
            M_d_loo <- sumcpp1(IM_d_loo, rows - 1, cols)
            M_loo <- M_d_loo
            # min_subset_loo<-apply(IM_subset_loo, c(1,2), min) min_index_loo<-apply(IM_subset_loo, c(1,2), which.min)
            # max_superset_loo<-apply(IM_superset_loo, c(1,2), max) max_index_loo<-apply(IM_superset_loo, c(1,2),
            # which.max)
            
            maxval <- maxcpp1(IM_superset_loo, rows - 1, cols)
            minval <- mincpp1(IM_subset_loo, rows - 1, cols)
            min_subset_loo <- minval$min
            min_index_loo <- minval$min_idx
            max_superset_loo <- maxval$max
            max_index_loo <- maxval$max_idx
            
            
            cell <- is.nan(M_d_loo) & is.finite(max_superset_loo)
            cell <- which(cell == TRUE)
            
            cell2 <- is.nan(M_d_loo) & is.finite(min_subset_loo)
            cell2 <- which(cell2 == TRUE)
            
            # does the removed drug need max averaging j_max<-which(cell==drug_index[i])
            j_max <- which(cell == identical_idx[i])
            # does the removed drug need min averaging j_min<-which(cell2==drug_index[i])
            j_min <- which(cell2 == identical_idx[i])
            
            if (length(j_max) != 0 && length(j_min) == 0) {
                # index for the cell
                cell_index <- cell[j_max]
                
                # row<-((cell_index-1) %% rows) + 1 col<-floor((cell_index-1) / rows)+1
                
                drug_sub_cell <- !is.infinite(IM_superset_loo[, cell_index])
                # the drug index which achieves max sensitivity
                index <- max_index_loo[cell_index]
                
                # the index of the dec of the drug with max sensitivity
                dec_maxsens <- drug_idx_loo[index]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets_small <- IM_subset_loo[, dec_maxsens] < max_superset_loo[cell_index]
                
                # find common item with drug_sub_cell and supersets_small
                common_cell <- which(drug_sub_cell & supersets_small)
                # cat(common_cell,'\n') max averaging
                if (length(common_cell) != 0) {
                  k <- 1
                  for (j in common_cell) {
                    max_superset_loo[cell_index] <- (max_superset_loo[cell_index] * k + sens_loo[j])/(k + 
                      1)
                    k <- k + 1
                  }
                }
                
                pred[i] <- (max_superset_loo[identical_idx[i]] + 1)/2
                error_predict[i] <- abs(pred[i] - sens[i])
                
            } else if (length(j_max) == 0 && length(j_min) != 0) {
                cell2_index <- cell2[j_min]
                
                drug_sub_cell <- !is.infinite(IM_subset_loo[, cell2_index])
                
                index <- min_index_loo[cell2_index]
                dec_minsens <- drug_idx_loo[index]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets_small <- IM_superset_loo[, dec_minsens] > min_subset_loo[cell2_index]
                # find common item with drug_sub_cell and supersets_small
                common_cell <- which(drug_sub_cell & supersets_small)
                if (length(common_cell) != 0) {
                  # max averaging
                  k <- 1
                  for (j in common_cell) {
                    min_subset_loo[cell2_index] <- (min_subset_loo[cell2_index] * k + sens_loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
                
                pred[i] <- (min_subset_loo[identical_idx[i]] + 0)/2
                error_predict[i] <- abs(pred[i] - sens[i])
            } else if (length(j_max) != 0 && length(j_min) != 0) {
                
                
                cell_index <- cell[j_max]
                
                # row<-((cell_index-1) %% rows) + 1 col<-floor((cell_index-1) / rows)+1
                
                drug_sub_cell <- !is.infinite(IM_superset_loo[, cell_index])
                # the drug index which achieves max sensitivity
                index <- max_index_loo[cell_index]
                
                # the dec of the drug with max sensitivity
                dec_maxsens <- drug_idx_loo[index]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets_small <- IM_subset_loo[, dec_maxsens] < max_superset_loo[cell_index]
                
                # find common item with drug_sub_cell and supersets_small
                common_cell <- which(drug_sub_cell & supersets_small)
                # cat(common_cell,'\n') max averaging
                if (length(common_cell) != 0) {
                  k <- 1
                  for (j in common_cell) {
                    max_superset_loo[cell_index] <- (max_superset_loo[cell_index] * k + sens_loo[j])/(k + 
                      1)
                    k <- k + 1
                  }
                }
                
                cell2_index <- cell2[j_min]
                
                drug_sub_cell <- !is.infinite(IM_subset_loo[, cell2_index])
                
                index <- min_index_loo[cell2_index]
                dec_minsens <- drug_idx_loo[index]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets_small <- IM_superset_loo[, dec_minsens] > min_subset_loo[cell2_index]
                # find common item with drug_sub_cell and supersets_small
                common_cell <- which(drug_sub_cell & supersets_small)
                if (length(common_cell) != 0) {
                  # max averaging
                  k <- 1
                  for (j in common_cell) {
                    min_subset_loo[cell2_index] <- (min_subset_loo[cell2_index] * k + sens_loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
                
                
                pred[i] <- (max_superset_loo[identical_idx[i]] + min_subset_loo[identical_idx[i]])/2
                error_predict[i] <- abs(pred[i] - sens[i])
                
            } else {
                # length(j_max)==0 && length(j_min)==0
                pred[i] <- M_loo[identical_idx[i]]
                error_predict[i] <- abs(pred[i] - sens[i])
            }
        }
        
    }
    return(error_predict)
} 
