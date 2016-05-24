Parse_Formula_Object <- function(formula,
                                 possible_structural_terms,
                                 possible_covariate_terms,
                                 possible_network_terms,
                                 raw_network = NULL,
                                 theta = NULL,
                                 terms_to_parse = "structural") {

  num_stats <- length(possible_structural_terms)
  # parse the formula
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable

  # get the actual response data
  if(is.null(raw_network)){
#     temp_net <- mget(lhs, envir = environment(),ifnotfound = list(not_found = NA))
#     if(class(temp_net) == "list"){
#       temp_net <- mget(lhs, envir = globalenv(),ifnotfound = list(not_found = NA))
#     }
    net <- dynGet(as.character(lhs),
                  ifnotfound = get(as.character(lhs)))
  }else{
    net <- raw_network
  }

  # check that the response is a matrix
  if (class(net) != "matrix") {
    stop("Response must be a matrix object.")
  }

  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", "", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, "\\+")[[1]]  # parse separate formula elements
  parsed_rhs <- vector (length = length(rhs), mode = "list")
  rhs_term_names <- rep("", length(rhs))
  alpha <- rep(1, length(rhs))
  threshold <- rep(0, length(rhs))
  for (i in 1:length(rhs)){
    parsed_rhs[[i]] <- parse_formula_term(rhs[i],
                                         possible_structural_terms,
                                         possible_covariate_terms,
                                         possible_network_terms)
    rhs_term_names[i] <- parsed_rhs[[i]]$term
    alpha[i] <- as.numeric(parsed_rhs[[i]]$weight)
    threshold[i] <- as.numeric(parsed_rhs[[i]]$threshold)
  }
  # if we are parsing the structural terms out of the formula
  if(terms_to_parse == "structural"){
    # remove all node level covariate terms
    remove <- which(rhs_term_names %in% possible_covariate_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    remove <- which(rhs_term_names %in% possible_network_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    # check that the names of the statistics match those that are possible
    possible <- 1:length(rhs_term_names)
    actual <- which(rhs_term_names %in% possible_structural_terms)
    if (length(possible) != length(actual)) {
      stop(paste("the specified structural term",
                 possible_structural_terms[setdiff(possible, actual)],
                 "is not an available statistic.", sep = " "))
    }
    # if theta is NULL, assume all ones
    if (is.null(theta) == TRUE) {
      theta <- rep(0, length(rhs_term_names))
    }

    # check that theta and the number of statistics are equal
    if (length(rhs_term_names) != length(theta)) {
      stop("'theta' must be the same length as the number of statistics")
    }
    stat.indx <- which(possible_structural_terms %in% rhs_term_names)
    statistics <- rep(0, num_stats)
    statistics[stat.indx] <- 1
    alphas <- rep(1, num_stats)
    thetas <- rep(0, num_stats)
    thresholds <- rep(0, num_stats)
    for (i in 1:length(rhs_term_names)) {
      alphas[which(rhs_term_names[i] == possible_structural_terms)] <- alpha[i]
      thetas[which(rhs_term_names[i] == possible_structural_terms)] <- theta[i]
      thresholds[which(rhs_term_names[i] == possible_structural_terms)] <- threshold[i]
    }
    return(list(net = net,
                statistics = statistics,
                alphas = alphas,
                thetas = thetas,
                thresholds = thresholds))
  }
  if(terms_to_parse == "covariate"){
    # if we are parsing covariate terms out of the formula
    remove <- which(rhs_term_names %in% possible_structural_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    remove <- which(rhs_term_names %in% possible_network_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    parsed_rhs <- append(parsed_rhs,list(network = net))
    return(parsed_rhs)
  }
  if(terms_to_parse == "network"){
    # if we are parsing network terms out of the formula
    remove <- which(rhs_term_names %in% possible_structural_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    remove <- which(rhs_term_names %in% possible_covariate_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    parsed_rhs <- append(parsed_rhs,list(network = net))
    return(parsed_rhs)
  }
}
