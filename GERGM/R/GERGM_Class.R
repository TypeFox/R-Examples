# an S4 class for gergm objects
setClass(Class = "gergm",
         representation = representation(
           network = "matrix",
           bounded.network = "matrix",
           formula = "formula",
           stats = "matrix",
           theta.coef = "data.frame",
           lambda.coef = "data.frame",
           weights = "numeric",
           num_nodes = "numeric",
           MCMC_output = "list",
           observed_network = "matrix",
           observed_bounded_network = "matrix",
           data_transformation = "array",
           stats_to_use = "numeric",
           previous_theta.coef = "data.frame",
           previous_lambda.coef = "data.frame",
           reduced_weights = "numeric",
           theta.par = "numeric",
           lambda.par = "numeric",
           theta_estimation_converged = "logical",
           simulated_vs_observed_p_values = "numeric",
           acceptable_fit = "logical",
           lambda_estimation_converged  = "logical",
           observed_simulated_t_test = "data.frame",
           console_output = "character",
           print_output = "logical",
           is_correlation_network = "logical",
           directed_network = "logical",
           simulation_only = "logical",
           thresholds = "numeric",
           BZ = "matrix",
           BZstdev = "numeric",
           transformation_type = "character",
           downweight_statistics_together = "logical",
           hyperparameter_optimization = "logical",
           target_accept_rate = "numeric",
           proposal_variance = "numeric",
           estimation_method = "character",
           number_of_simulations = "numeric",
           thin = "numeric",
           burnin = "numeric",
           MPLE_gain_factor = "numeric",
           simulated_statistics_for_GOF = "data.frame"
         ),
         validity = function(object) {
           if (!"matrix" %in% class(object@network) & is.null(object@network)
               == FALSE) {
             stop("'network' must be a 'matrix' object or 'NULL'.")
           }
           if (!"matrix" %in% class(object@bounded.network)) {
             stop("'bounded.network' must be a 'matrix' object.")
           }
           if (!is.data.frame(object@theta.coef)) {
             stop("'theta.coef' must be a 'data.frame' object.")
           }
           if (!"formula" %in% class(object@formula)) {
             stop("'formula' is not a 'formula' object.")
           }
           if (!is.data.frame(object@lambda.coef) & is.null(object@lambda.coef)
               == FALSE) {
             stop("'lambda.coef' must be a 'data.frame' object or 'NULL'.")
           }
           if (!"matrix" %in% class(object@stats)) {
             stop("'stats' must be a 'matrix' object.")
           }
           if (!is.numeric(object@weights) & is.null(object@weights)
               == FALSE) {
             stop("'weights' must be a 'data.frame' object or 'NULL'.")
           }
           if (!is.numeric(object@num_nodes) & is.null(object@num_nodes)
               == FALSE) {
             stop("'num_nodes' must be a 'data.frame' object or 'NULL'.")
           }
           if (!"list" %in% class(object@MCMC_output)) {
             stop("'MCMC_output' must be a 'list' object.")
           }
           if (!"matrix" %in% class(object@observed_network)) {
             stop("'observed_network' must be a 'matrix' object.")
           }
           if (!"matrix" %in% class(object@observed_bounded_network)) {
             stop("'observed_bounded_network' must be a 'matrix' object.")
           }
           if (!"array" %in% class(object@data_transformation)) {
             stop("'data_transformation' must be a 'array' object.")
           }
           if (!is.numeric(object@stats_to_use) & is.null(object@stats_to_use)
               == FALSE) {
             stop("'stats_to_use' must be a 'numeric' object or 'NULL'.")
           }
           if (!is.data.frame(object@previous_lambda.coef) & is.null(object@previous_lambda.coef)
               == FALSE) {
             stop("'previous_lambda.coef' must be a 'data.frame' object or 'NULL'.")
           }
           if (!is.data.frame(object@previous_theta.coef) & is.null(object@previous_theta.coef)
               == FALSE) {
             stop("'previous_theta.coef' must be a 'data.frame' object or 'NULL'.")
           }
           if (!is.numeric(object@reduced_weights) & is.null(object@reduced_weights)
               == FALSE) {
             stop("'reduced_weights' must be a 'numeric' object or 'NULL'.")
           }
           if (!is.numeric(object@theta.par) & is.null(object@theta.par)
               == FALSE) {
             stop("'theta.par' must be a 'numeric' object or 'NULL'.")
           }
           if (!is.numeric(object@lambda.par) & is.null(object@lambda.par)
               == FALSE) {
             stop("'lambda.par' must be a 'numeric' object or 'NULL'.")
           }
           if (!is.numeric(object@simulated_vs_observed_p_values)
               & is.null(object@simulated_vs_observed_p_values) == FALSE) {
             stop("'simulated_vs_observed_p_values' must be a 'numeric' object or 'NULL'.")
           }
           if (!"logical" %in% class(object@theta_estimation_converged)) {
             stop("'theta_estimation_converged' must be a 'logical' value.")
           }
           if (!"logical" %in% class(object@lambda_estimation_converged)) {
             stop("'lambda_estimation_converged' must be a 'logical' value.")
           }
           if (!"logical" %in% class(object@acceptable_fit)) {
             stop("'acceptable_fit' must be a 'logical' value.")
           }
           if (!is.data.frame(object@observed_simulated_t_test)
               & is.null(object@observed_simulated_t_test)
               == FALSE) {
             stop("'observed_simulated_t_test' must be a 'data.frame' object or 'NULL'.")
           }
           if (!is.character(object@console_output) & is.null(object@console_output)
               == FALSE) {
             stop("'console_output' must be a 'character' object or 'NULL'.")
           }
           if (!is.character(object@console_output) & is.null(object@console_output)
               == FALSE) {
             stop("'console_output' must be a 'character' object or 'NULL'.")
           }
           return(TRUE)
         }
)

# define coef for pretty output of gergm object
setMethod(f = "coef", signature = "gergm", definition = function(object, ...) {
  return(list(Theta  = object@theta.coef, Lambda = object@lambda.coef))
}
)

# define 'show' to get a pretty output of the gergm
setMethod(f = "show", signature = "gergm", definition = function(object){
  message("Theta:")
  print(object@theta.coef)
  message("Lambda:")
  print(object@lambda.coef)
  message("Weights:")
  print(object@weights)
  message("Network Statistics:")
  print(object@stats)
}
)
