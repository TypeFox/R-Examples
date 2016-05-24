# create a gergm from a formula object (getgergm)
Create_GERGM_Object_From_Formula <- function(object,
                                             theta.coef,
                                             possible_structural_terms,
                                             possible_covariate_terms,
                                             possible_network_terms,
                                             raw_network,
                                             together = 1,
                                             transform.data = NULL,
                                             lambda.coef = NULL,
                                             transformation_type,
                                             is_correlation_network = FALSE,
                                             is_directed = TRUE
                                             ){

  res1 <- Parse_Formula_Object(object,
                               possible_structural_terms,
                               possible_covariate_terms,
                               possible_network_terms,
                               raw_network = raw_network,
                               theta = theta.coef,
                               terms_to_parse = "structural")
  thetas <- res1$thetas
  network <- res1$net
  alphas <- res1$alphas
  statistics <- res1$statistics
  thresholds <- res1$thresholds

  # for now we are not going to allow any covariates
  if(is_correlation_network){
    if(!is.null(lambda.coef)){
      stop("Covariate effects are currently not supported for correlation networks. Please respecify without covariates.")
    }
  }

  # create the network based on the transform family
  # if there are no lambda.coefficients, we assume there is no transformation
  # if there is a transformation specified, transform the observed network
  if (!is.null(lambda.coef) == TRUE){
    if(transformation_type == "logcauchy" | transformation_type == "lognormal"){
      if(min(network) <= 0){
        stop(paste("You have selected either a log-Cauchy or log-normal transformation but you have provided a network with values that are less than or equal to zero. Please ensure that the minimum value of the network you provide is greater than zero, or select a cauchy or normal transformation. The minimum value of the network provided is:",min(network)))
      }
      network <- log(network)
    }
    beta <- lambda.coef[1:(length(lambda.coef) - 1)]
    sig <- 0.01 + exp(lambda.coef[length(lambda.coef)])
    BZ <- 0
    if(is.na(dim(transform.data)[3]) == TRUE){
      BZ = BZ + beta * transform.data
    }
    if(!is.na(dim(transform.data)[3]) == TRUE){
      for (j in 1:(dim(transform.data)[3])) {
        BZ <- BZ + beta[j] * transform.data[, , j]
      }
    }
    if(transformation_type == "logcauchy" | transformation_type == "cauchy"){
      bounded.network <- pst(network, BZ, sig, 1)
    }
    if(transformation_type == "lognormal" | transformation_type == "gaussian"){
      bounded.network <- pst(network, BZ, sig, Inf)
    }

  }
  if (is.null(lambda.coef) == TRUE) {
    bounded.network <- network
    lambda.coef <- as.data.frame(0)
  }
  if (is.null(lambda.coef) != TRUE){
    lambda.coef <- as.data.frame(rbind(lambda.coef,NA))
    rownames(lambda.coef) <- c("est", "se")
  }

  # if we are providing a correlation network, transform it
  if(is_correlation_network){
    diag(network) <- 1
    print(round(network,2))
    bounded.network <- transform.correlations(network)
  }


  thetas <- t(as.matrix(thetas))
  thetas <- rbind(thetas, NA)
  colnames(thetas) <- possible_structural_terms
  rownames(thetas) <- c("est", "se")
  thetas <- as.data.frame(thetas)

  object <- Create_GERGM_Object(network = network,
                                bounded.network = bounded.network,
                                formula = object,
                                thetas = thetas,
                                lambda = lambda.coef,
                                alpha = alphas,
                                together = together,
                                possible.stats = possible_structural_terms,
                                thresholds = thresholds)
  object@stats_to_use <- statistics
  return(object)
}
