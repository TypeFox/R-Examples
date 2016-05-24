#' Find the same columns from two matrices
#' 
#' A function to find the same columns from two matrices
#' 
#' @param profile the drug-target interaction data matrix
#' @param selected_list the selected drug-target matrix
#' @param kinase_name the names of the targets
#' @return a vector of combined selected target names
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' x<-data.frame(tyner_interaction_binary)
#' kinase_names<-dimnames(tyner_interaction_binary)
#' float<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[,1])
#' k_select <- float$k_sel
#' select_kinase_names <- findSameSet(x, k_select, kinase_names)
#' }
findSameSet <- function(profile, selected_list, kinase_name) {
    # parameter selected_list: the index for selected kinase
    selected_kinases <- c()
    for (i in 1:length(selected_list)) {
        same_set <- findSameCol(profile, profile[, selected_list[i]])
        selected_kinases <- c(selected_kinases, paste(kinase_name[which(same_set == 1)], collapse = ";"))
    }
    return(selected_kinases)
} 
