Prepare_Network_and_Covariates <- function(formula,
                                    possible_structural_terms,
                                    possible_covariate_terms,
                                    possible_network_terms,
                                    covariate_data = NULL,
                                    normalization_type = c("log","division"),
                                    is_correlation_network = FALSE,
                                    is_directed = FALSE
                                    ){

  node_covariates_list <- Parse_Formula_Object(formula = formula,
     possible_structural_terms =possible_structural_terms,
     possible_covariate_terms =possible_covariate_terms,
     possible_network_terms =possible_network_terms,
     raw_network = NULL,
     theta = NA,
     terms_to_parse = "covariate")
  network_covariates_list <- Parse_Formula_Object(formula = formula,
     possible_structural_terms =possible_structural_terms,
     possible_covariate_terms =possible_covariate_terms,
     possible_network_terms =possible_network_terms,
     raw_network = NULL,
     theta = NA,
     terms_to_parse = "network")

  raw_network <- node_covariates_list$network

  # Calculate number of nodes in network
  num_nodes <- nrow(raw_network)
  node_names <- colnames(raw_network)

  # Make sure that the diagnoal terms in the network are zeros
  diag(raw_network) <- 0

  # Check if network is square
  if (nrow(raw_network) != ncol(raw_network)) {
    stop("Network must be a square matrix or data frame!")
  }

  # Just to make sure the code does not break!
  if(is.null(normalization_type)){
    normalization_type <- "log"
  }

  # If network is a data.frame, make it a matrix
  if (is.data.frame(raw_network)) {
    raw_network <- as.matrix(raw_network)
  }

  #determine whether covariates were provided and the total number of covariates
  if(length(node_covariates_list) < 2){
    node_covariates_provided <- FALSE
  }else{
    node_covariates_provided <- TRUE
  }
  if(length(network_covariates_list) < 2){
    network_covariates_provided <- FALSE
  }else{
    network_covariates_provided <- TRUE
  }

  # transformed covariates will have one additional slice of all ones (we are omitting this for now)
  num_covariates <- 0

  if(!node_covariates_provided){
    cat("No node level covariates provided.\n")
  }else{
    node_covariates_list <- node_covariates_list[-length(node_covariates_list)]
    for(i in 1:length(node_covariates_list)){
      type <- node_covariates_list[[i]]$term
      if(type == "sender" | type == "receiver"| type == "absdiff"| type == "nodecov" | type == "intercept" | type == "nodematch"){
        num_covariates <- num_covariates + 1
      }else if(type == "nodemix"){
        # need to get the number of levels
        covar <- node_covariates_list[[i]]$covariate
        index <- which(colnames(covariate_data) == covar)
        node_covariates_list[[i]]$levels <- unique(covariate_data[,index])
        num_levels <- length(node_covariates_list[[i]]$levels)
        node_covariates_list[[i]]$num_levels <- num_levels
        if(node_covariates_list[[i]]$base == "NULL"){
          # if the user specified base = NULL then we use all levels for matching
          num_covariates <- num_covariates + num_levels*num_levels
          node_covariates_list[[i]]$base_index <- 0
        }else{
          num_covariates <- num_covariates + num_levels*num_levels - 1
          base_index <- which(node_covariates_list[[i]]$levels == node_covariates_list[[i]]$base)
          node_covariates_list[[i]]$base_index <- base_index
        }
      }else{
        stop(paste("You specified a node level covariate term:",node_covariates_list[[i]]$term, "Node level covariate effects must be one of: 'sender', 'receiver', 'absdiff', 'nodecov', 'nodematch', or 'nodemix' , please respecify."))
      }
    }
    cat("You have specified", num_covariates - 1,"node level covariate effects.\n")
  }

  #determine the type and number of user provided network covariates (will not be altered)
  if(!network_covariates_provided){
    cat("No network covariates provided.\n")
    num_additional_covars <- 0
  }else{
    network_covariates_list <- network_covariates_list[-length(network_covariates_list)]
    num_covariates <- num_covariates + length(network_covariates_list)
    num_additional_covars <- length(network_covariates_list)
    cat("You have provided",num_additional_covars,"network covariates.\n")
  }

  generate_covariate_effect_matrix <- function(num_nodes,
                                               node_names,
                                               covariates,
                                               covariate_column,
                                               effect_type,
                                               level = NA,
                                               level2 = NA){
    return_matrix <- matrix(0,num_nodes,num_nodes)
    if(effect_type == "intercept"){
      return_matrix <- matrix(1,num_nodes,num_nodes)
      diag(return_matrix) <- 0
    }
    if(effect_type == "sender"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
            return_matrix[j,k] <- covariates[row,covariate_column]
          }
        }
      }
    }
    if(effect_type == "receiver"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
            return_matrix[j,k] <- covariates[row,covariate_column]
          }
        }
      }
    }
    if(effect_type == "absdiff"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
            row2 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
            return_matrix[j,k] <- abs(covariates[row1,covariate_column] - covariates[row2,covariate_column])
          }
        }
      }
    }
    if(effect_type == "nodecov"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
            row2 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
            return_matrix[j,k] <- covariates[row1,covariate_column] + covariates[row2,covariate_column]
          }
        }
      }
    }
    if(effect_type == "nodematch"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
            row2 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
            # handle both numeric and categorical values
            if(is.numeric(covariates[row1,covariate_column]) & is.numeric(covariates[row2,covariate_column])){
              check <- abs(covariates[row1,covariate_column] - covariates[row2,covariate_column])
            }else{
              check <- 1
              if(covariates[row1,covariate_column] == covariates[row2,covariate_column]){
                check <- 0
              }
            }
            if(check == 0){
              return_matrix[j,k] <- 1
            }
          }
        }
      }

    }
    if(effect_type == "nodemix"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            col1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
            row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
            colval <- covariates[col1,covariate_column]
            rowval <- covariates[row1,covariate_column]
            if(level == colval & level2 == rowval){
              return_matrix[j,k] <- 1
            }
          }
        }
      }
    }
    return(return_matrix)
  }

  #########################################################
  # Construct transformed_covariates: An array of covariates which parameterize
  # the latent space.

  # Generate array which covaries will be transformed into
  transformed_covariates <- array(0,dim=c(num_nodes,num_nodes,num_covariates))

  #omit this since we are adding in an intercept term
  #transformed_covariates[,,1] <- 1
  # Set a slice counter to keep track of where we should add the covariates
  # in the resulting covariate array.
  slice_counter <- 1

  # generate vector to store covariate names
  if(!node_covariates_provided & !network_covariates_provided){
    # do nothing
  }else{
    slice_names <- rep("",length = num_covariates)

    #remove since we are specifying an intercept
    #slice_names[1] <- "network"
  }



  #1. generate sender and reciever effects
  if(node_covariates_provided){

    # Determine if type_of_effect and covariates_to_use have the same length or
    # if no covariates_to_use is provided, that the number of columns in
    # covariate_data is the same length.
    node_covariates_list[[i]]$num_levels
    # Loop through covariates
    for(i in 1:length(node_covariates_list)){
      if(node_covariates_list[[i]]$term == "intercept"){
        col_index <- 0
      }else{
        col_index <- which(tolower(colnames(covariate_data)) ==
                             tolower(node_covariates_list[[i]]$covariate))
      }
      if(length(col_index) == 0){
        stop(paste("There is no matching column name in covariate_data for:",tolower(node_covariates_list[[i]]$covariate)))
      }

      if(node_covariates_list[[i]]$term == "nodemix"){
        levels_to_include <- node_covariates_list[[i]]$num_levels
        for (j in 1:levels_to_include){
          for (k in 1:levels_to_include){
            if(node_covariates_list[[i]]$levels[j] == node_covariates_list[[i]]$base & node_covariates_list[[i]]$levels[k] == node_covariates_list[[i]]$base){
              #do nothing since we do not include a term for the base
            }else{
              add <- generate_covariate_effect_matrix(
                num_nodes = num_nodes,
                node_names = node_names,
                covariates = covariate_data,
                covariate_column = col_index,
                effect_type = "nodemix",
                level = node_covariates_list[[i]]$levels[j],
                level2 = node_covariates_list[[i]]$levels[k])
              #print(add)
              transformed_covariates[,,slice_counter] <- add
              slice_names[slice_counter] <- paste(
                node_covariates_list[[i]]$covariate,
                node_covariates_list[[i]]$term,
                node_covariates_list[[i]]$levels[j],
                node_covariates_list[[i]]$levels[k],
                sep = "_")
              slice_counter <- slice_counter + 1
            }
          }
        }
      }else{
          add <- generate_covariate_effect_matrix(
            num_nodes = num_nodes,
            node_names = node_names,
            covariates = covariate_data,
            covariate_column = col_index,
            effect_type = node_covariates_list[[i]]$term)
        transformed_covariates[,,slice_counter] <- add
        if(node_covariates_list[[i]]$term == "intercept"){
          slice_names[slice_counter] <- "intercept"
        }else{
          slice_names[slice_counter] <- paste(
            node_covariates_list[[i]]$covariate,
            node_covariates_list[[i]]$term,
            sep = "_")
        }
        slice_counter <- slice_counter + 1
      }
    } # End of loop over node level covariates
  } # End of condition that we have node level covariates



  #2. tack on any user supplied effects
  if(network_covariates_provided){
    for(i in 1:num_additional_covars){
      transformed_covariates[,,slice_counter] <- network_covariates_list[[i]]$network_matrix_object
      slice_names[slice_counter] <- paste(
        network_covariates_list[[i]]$network,
        network_covariates_list[[i]]$term,
        sep = "_")
      slice_counter <- slice_counter + 1
    }
  }




  #3. return list object
  #
  #
  if(!is_directed){
    if(isSymmetric(raw_network) == FALSE){
      warning("You provided an asymmetric network when you specified that the network was undirect. Setting the upper triangle equal to the lower triangle. If the lower is larger, then it will be coppied, otherwise the upper will be coppied. To avoid problems, please provide a symmetric matrix.")

      lower_to_upper <- function(m) {
        m[upper.tri(m)] <- t(m)[upper.tri(m)]
        m
      }

      upper_to_lower <- function(m) {
        m[lower.tri(m)] <- t(m)[lower.tri(m)]
        m
      }

      if(sum(upper.tri(raw_network)) > sum(lower.tri(raw_network))){
        raw_network <- upper_to_lower(raw_network)
      }else{
        raw_network <- lower_to_upper(raw_network)
      }
      cat("Symmetrized Raw Network...\n")
      print(round(raw_network,1))
    }
  }

  if(!node_covariates_provided & !network_covariates_provided){
    # If no covariates were provided, then make sure the network lives on the
    # [0,1] interval and standardize it by one of the provided methods if it
    # does not. Then return the network and no covariates.

    # deal with the case where we have a correlation matrix provided
    if(is_correlation_network){
      network <- raw_network
      diag(network) <- 1
    }else{
      if(min(raw_network) < 0){
        raw_network <- raw_network - min(raw_network)
      }
      if(max(raw_network) > 1){
        if(normalization_type[1] == "log"){
          raw_network <- raw_network + 1
          network <- log(raw_network)
          network <- network/max(network)
        }

        if(normalization_type[1] == "division"){
          network <- raw_network/max(raw_network)
        }

        diag(network) <- 0
      }else{
        network <- raw_network
      }
    }
    cat("Transformed Network...\n")
    print(round(network,3))
    return(list(network = network))
  }else{
    #4. standardize covariates
    if(num_covariates > 1){
      for(i in 2:num_covariates){
        transformed_covariates[,,i] <- (transformed_covariates[,,i]-mean(c(transformed_covariates[,,i])))/sd(c(transformed_covariates[,,i]))
      }
    }
    #assign the dimnames to the array object
    dimnames(transformed_covariates) <- list(node_names,
                                             node_names,
                                             slice_names)
    return(list(network = raw_network,
                transformed_covariates = transformed_covariates,
                gpar.names = slice_names))
  }

} # End of function definition.
