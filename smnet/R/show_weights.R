show_weights <- function(SSNobject, adjacency, netID = 1){
  rid_data        <- SSNobject@data
  ord             <- order(rid_data$rid)
  rid_data        <- rid_data[ord,]
  rid_data        <- rid_data[rid_data$netID == netID,]
  variable_names  <- names(rid_data)
  weight_type     <- vector("character", length = length(variable_names))
  # now cycle through columns in turn and work out which is a weight
  for(i in 1:length(variable_names)){
    current_var <- as.vector(as.matrix(rid_data[, i]))
    if(class(current_var) %in% c("integer", "numeric")){
      weight_type[i] <- check_weight(current_var, adj = adjacency, silent = T) 
    } else {
      weight_type[i] <- "unrecognised"
    }
  }
  weight.poss     <- variable_names[which(!weight_type == "unrecognised")]
  weight.poss     <- paste(paste(weight.poss, collapse = ", "), ".", sep = "")
  weight.poss     <- paste("\nThe following recognised weights were found: \n", 
                           "-------------------------------------------- \n", 
                           weight.poss, sep = "")
  cat(weight.poss)
}