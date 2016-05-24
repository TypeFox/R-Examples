#' @title A Function to estimate a GERGM.
#' @description The main function provided by the package.
#'
#' @param formula A formula object that specifies the relationship between
#' statistics and the observed network. Currently, the user may specify a model
#'  using any combination of the following statistics: `out2stars(alpha = 1)`,
#'  `in2stars(alpha = 1)`, `ctriads(alpha = 1)`, `mutual(alpha = 1)`,
#'  `ttriads(alpha = 1)`, `absdiff(covariate = "MyCov")`,
#'  `edgecov(covariate = "MyCov")`, `sender(covariate = "MyCov")`,
#'  `reciever(covariate = "MyCov")`, `nodematch(covariate)`,
#'   `nodemix(covariate, base = "MyBase")`, `netcov(network)`. To use
#'  exponential downweighting for any of the network level terms, simply
#'  specify a value for alpha less than 1. The `(alpha = 1)` term may be omitted
#'   from the structural terms if no exponential downweighting is required. In
#'   this case, the terms may be provided as: `out2star`, `in2star`, `ctriads`,
#'   `recip`, `ttriads`. If the network is undirected the user may only specify
#'   the following terms: `twostars(alpha = 1)`,  `ttriads(alpha = 1)`,
#'   `absdiff(covariate = "MyCov")`, `edgecov(covariate = "MyCov")`,
#'   `sender(covariate = "MyCov")`, `nodematch(covariate)`,
#'   `nodemix(covariate, base = "MyBase")`, `netcov(network)`. An intercept
#'   term is included by default, but can be omitted by setting
#'   omit_intercept_term = TRUE. If the user specifies
#'   `nodemix(covariate, base = NULL)`, then all levels of the covariate
#'   will be matched on.
#' @param covariate_data A data frame containing node level covariates the user
#' wished to transform into sender or reciever effects. It must have row names
#' that match every entry in colnames(raw_network), should have descriptive
#' column names.  If left NULL, then no sender or reciever effects will be
#' added.
#' @param normalization_type If only a raw_network is provided and
#' omit_intercept_term = TRUE then, the function
#' will automatically check to determine if all edges fall in the [0,1] interval.
#' If edges are determined to fall outside of this interval, then a trasformation
#' onto the interval may be specified. If "division" is selected, then the data
#' will have a value added to them such that the minimum value is atleast zero
#' (if necessary) and then all edge values will be divided by the maximum to
#' ensure that the maximum value is in [0,1]. If "log" is selected, then the data
#' will have a value added to them such that the minimum value is atleast zero
#' (if necessary), then 1 will be added to all edge values before they are logged
#' and then divided by the largest value, again ensuring that the resulting
#' network is on [0,1]. Defaults to "log" and need not be set to NULL if
#' providing covariates as it will be ignored.
#' @param network_is_directed Logical specifying whether or not the observed
#' network is directed. Default is TRUE.
#' @param use_MPLE_only Logical specifying whether or not only the maximum pseudo
#' likelihood estimates should be obtained. In this case, no simulations will be
#' performed. Default is FALSE.
#' @param transformation_type Specifies how covariates are transformed onto the
#' raw network. When working with heavly tailed data that are not strictly
#' positive, select "Cauchy" to transform the data using a Cauchy distribution.
#' If data are strictly positive and heavy tailed (such as financial data) it is
#' suggested the user select "LogCauchy" to perform a Log-Cauchy transformation
#' of the data. For a tranformation of the data using a Gaussian distribution,
#' select "Gaussian" and for strictly positive raw networks, select "LogNormal".
#' The Default value is "Cauchy".
#' @param estimation_method Simulation method for MCMC estimation. Default is
#' "Gibbs" which will generally be faster with well behaved networks but will not
#' allow for exponential downweighting.
#' @param maximum_number_of_lambda_updates Maximum number of iterations of outer
#' MCMC loop which alternately estimates transform parameters and ERGM
#' parameters. In the case that data_transformation = NULL, this argument is
#' ignored. Default is 10.
#' @param maximum_number_of_theta_updates Maximum number of iterations within the
#' MCMC inner loop which estimates the ERGM parameters. Default is 100.
#' @param number_of_networks_to_simulate Number of simulations generated for
#' estimation via MCMC. Default is 500.
#' @param thin The proportion of samples that are kept from each simulation. For
#' example, thin = 1/200 will keep every 200th network in the overall simulated
#' sample. Default is 1.
#' @param proposal_variance The variance specified for the Metropolis Hastings
#' simulation method. This parameter is inversely proportional to the average
#' acceptance rate of the M-H sampler and should be adjusted so that the average
#' acceptance rate is approximately 0.25. Default is 0.1.
#' @param downweight_statistics_together Logical specifying whether or not the
#' weights should be applied inside or outside the sum. Default is TRUE and user
#' should not select FALSE under normal circumstances.
#' @param MCMC_burnin Number of samples from the MCMC simulation procedure that
#' will be discarded before drawing the samples used for estimation.
#' Default is 100.
#' @param seed Seed used for reproducibility. Default is 123.
#' @param convergence_tolerance Threshold designated for stopping criterion. If
#' the difference of parameter estimates from one iteration to the next all have
#' a p -value (under a paired t-test) greater than this value, the parameter
#' estimates are declared to have converged. Default is 0.01.
#' @param MPLE_gain_factor Multiplicative constant between 0 and 1 that controls
#' how far away the initial theta estimates will be from the standard MPLEs via
#' a one step Fisher update. In the case of strongly dependent data, it is
#' suggested to use a value of 0.10. Default is 0.
#' @param acceptable_fit_p_value_threshold A p-value threshold for how closely
#' statistics of observed network conform to statistics of networks simulated
#' from GERGM parameterized by converged final parameter estimates. Default value
#' is 0.05.
#' @param force_x_theta_updates Defaults to 1 where theta estimation is not
#' allowed to converge until thetas have updated for x iterations . Useful when
#' model is not degenerate but simulated statistics do not match observed network
#' well when algorithm stops after first y updates.
#' @param force_x_lambda_updates Defaults to 1 where lambda estimation is not
#' allowed to converge until lambdas have updated for x iterations . Useful when
#' model is not degenerate but simulated statistics do not match observed network
#' well when algorithm stops after first y updates.
#' @param output_directory The directory where you would like output generated
#' by the GERGM estimation proceedure to be saved (if output_name is specified).
#' This includes, GOF, trace, and parameter estimate plots, as well as a summary
#' of the estimation proceedure and an .Rdata file containing the GERGM object
#' returned by this function. May be left as NULL if the user would prefer all
#' plots be printed to the graphics device.
#' @param output_name The common name stem you would like to assign to all
#' objects output by the gergm function. Default value of NULL will not save any
#' output directly to .pdf files, it will be printed to the console instead. Must
#' be a character string or NULL. For example, if "Test" is supplied as the
#' output_name, then 4 files will be output: "Test_GOF.pdf", "Test_Parameter_Estim
#' ates.pdf", "Test_GERGM_Object.Rdata", "Test_Estimation_Log.txt", and
#' "Test_Trace_Plot.pdf"
#' @param generate_plots Defaults to TRUE, if FALSE, then no diagnostic or
#' parameter plots are generated.
#' @param verbose Defaults to TRUE (providing lots of output while model is
#' running). Can be set to FALSE if the user wishes to see less output.
#' @param omit_intercept_term Defualts to FALSE, can be set to TRUE if the
#' user wishes to omit the model intercept term.
#' @param hyperparameter_optimization Logical indicating whether automatic
#' hyperparameter optimization should be used. Defaults to FALSE. If TRUE, then
#' the algorithm will automatically seek to find an optimal burnin and number of
#' networks to simulate, and if using Metropolis Hasings, will attempt to select
#' a proposal variance that leads to a acceptance rate within +-0.05 of
#' target_accept_rate. Furthermore, if degeneracy is detected, the algorithm
#' will attempt to adress the issue automatically. WARNING: This feature is
#' experimental, and may greatly increase runtime. Please monitor console
#' output!
#' @param target_accept_rate The target Metropolis Hastings acceptance rate.
#' Defaults to 0.25
#' @param ... Optional arguments, currently unsupported.
#' @return A gergm object containing parameter estimates.
#' @examples
#' \dontrun{
#' set.seed(12345)
#' net <- matrix(rnorm(100,0,20),10,10)
#' colnames(net) <- rownames(net) <- letters[1:10]
#' formula <- net ~  mutual + ttriads
#'
#' test <- gergm(formula,
#'               normalization_type = "division",
#'               network_is_directed = TRUE,
#'               use_MPLE_only = FALSE,
#'               estimation_method = "Metropolis",
#'               number_of_networks_to_simulate = 40000,
#'               thin = 1/10,
#'               proposal_variance = 0.5,
#'               downweight_statistics_together = TRUE,
#'               MCMC_burnin = 10000,
#'               seed = 456,
#'               convergence_tolerance = 0.01,
#'               MPLE_gain_factor = 0,
#'               force_x_theta_updates = 4)
#' }
#' @export
gergm <- function(formula,
                  covariate_data = NULL,
                  normalization_type = c("log","division"),
                  network_is_directed = c(TRUE, FALSE),
                  use_MPLE_only = c(FALSE, TRUE),
                  transformation_type = c("Cauchy","LogCauchy","Gaussian","LogNormal"),
                  estimation_method = c("Gibbs", "Metropolis"),
                  maximum_number_of_lambda_updates = 10,
                  maximum_number_of_theta_updates = 10,
                  number_of_networks_to_simulate = 500,
                  thin = 1,
                  proposal_variance = 0.1,
                  downweight_statistics_together = TRUE,
                  MCMC_burnin = 100,
                  seed = 123,
                  convergence_tolerance = 0.01,
                  MPLE_gain_factor = 0,
                  acceptable_fit_p_value_threshold = 0.05,
                  force_x_theta_updates = 1,
                  force_x_lambda_updates = 1,
                  output_directory = NULL,
                  output_name = NULL,
                  generate_plots = TRUE,
                  verbose = TRUE,
                  omit_intercept_term = FALSE,
                  hyperparameter_optimization = FALSE,
                  target_accept_rate = 0.25,
                  ...
                  ){

  # pass in experimental features through elipsis
  using_correlation_network = FALSE
  object <- as.list(substitute(list(...)))[-1L]
  if (length(object) > 0) {
    if (!is.null(object$using_correlation_network)) {
      if (object$using_correlation_network) {
        using_correlation_network <- TRUE
        cat("Using experimental correlation network feature...\n")
      }
    }
  }

  # hard coded possible stats
  possible_structural_terms <- c("out2stars",
                                 "in2stars",
                                 "ctriads",
                                 "mutual",
                                 "ttriads",
                                 "edges")
  possible_structural_terms_undirected <- c("edges",
                                            "twostars",
                                            "ttriads")
  possible_covariate_terms <- c("absdiff",
                                "nodecov",
                                "nodematch",
                                "sender",
                                "receiver",
                                "intercept",
                                "nodemix")
  possible_network_terms <- "netcov"
  possible_transformations <- c("cauchy",
                                "logcauchy",
                                "gaussian",
                                "lognormal")

  # check terms for undirected network
  if (!network_is_directed) {
    formula <- parse_undirected_structural_terms(
      formula,
      possible_structural_terms,
      possible_structural_terms_undirected)
  }

  # automatically add an intercept term unless omit_intercept_term is TRUE
  if (!omit_intercept_term) {
    formula <- add_intercept_term(formula)
    #check for an edges statistic
    form <- as.formula(formula)
    parsed <- deparse(form)
    if (length(parsed) > 1) {
      parsed <- paste0(parsed, collapse = " ")
    }
    if (grepl("edges",parsed)) {
      stop("You may not specify an edges statistic if omit_intercept_term == FALSE as this will introduce two identical intercept terms and instability in the model. An intercept term is automatically added in the lambda transformation step unless omit_intercept_term == TRUE, and we have found this method of adding an intercept to be less prone to degeneracy.")
    }
  }

  # set logical values for whether we are using MPLE only, whether the network
  # is directed, and which estimation method we are using as well as the
  # transformation type
  use_MPLE_only <- use_MPLE_only[1] #default is FALSE
  network_is_directed <- network_is_directed[1] #default is TRUE
  estimation_method <- estimation_method[1] #default is Gibbs
  transformation_type <- transformation_type[1] #default is "Cauchy"
  transformation_type <- tolower(transformation_type)
  normalization_type <- normalization_type[1]

  # if we are using a correlation network, then the network must be undirected.
  if (using_correlation_network) {
    network_is_directed <- FALSE
  }

  if (is.null(output_directory) & !is.null(output_name)) {
    stop("You have specified an output file name but no output directory. Please
         specify both or neither.")
  }

  if (length(which(possible_transformations %in% transformation_type  == T)) != 1) {
    stop("You have specified a transformation that is not recognized. Please
         specify one of: Cauchy, LogCauchy, Gaussian, or LogNormal")
  }

  #make sure proposal variance is greater than zero
  if (proposal_variance <= 0.001) {
    proposal_variance <- 0.001
    cat("You supplied a proposal variance that was less than or equal to zero. It has been reset to 0.001, considder respecifying...\n")
  }

  formula <- as.formula(formula)

  #0. Prepare the data
  Transformed_Data <- Prepare_Network_and_Covariates(
     formula,
     possible_structural_terms,
     possible_covariate_terms,
     possible_network_terms,
     covariate_data = covariate_data,
     normalization_type = normalization_type,
     is_correlation_network = using_correlation_network,
     is_directed = network_is_directed)

  data_transformation <- NULL
  if (!is.null(Transformed_Data$transformed_covariates)) {
    data_transformation <- Transformed_Data$transformed_covariates
  }
  gpar.names <- c(Transformed_Data$gpar.names, "dispersion")

  #1. Create GERGM object from network

  GERGM_Object <- Create_GERGM_Object_From_Formula(
     formula,
     theta.coef = NULL,
     possible_structural_terms,
     possible_covariate_terms,
     possible_network_terms,
     raw_network = Transformed_Data$network,
     together = 1,
     transform.data = data_transformation,
     lambda.coef = NULL,
     transformation_type = transformation_type,
     is_correlation_network = using_correlation_network,
     is_directed = network_is_directed
     )

  GERGM_Object@theta_estimation_converged <- FALSE
  GERGM_Object@lambda_estimation_converged <- FALSE
  GERGM_Object@observed_network  <- GERGM_Object@network
  GERGM_Object@observed_bounded_network <- GERGM_Object@bounded.network
  GERGM_Object@simulation_only <- FALSE
  GERGM_Object@transformation_type <- transformation_type
  GERGM_Object@downweight_statistics_together <- downweight_statistics_together
  GERGM_Object@directed_network <- network_is_directed


  if (!is.null(data_transformation)) {
    GERGM_Object@data_transformation <- data_transformation
  }

  if (is.null(output_name)) {
    GERGM_Object@print_output <- FALSE
  }else{
    GERGM_Object@print_output <- TRUE
  }

  # if we are using a correlation network then set field to TRUE.
  GERGM_Object@is_correlation_network <- using_correlation_network


  # set adaptive metropolis parameters
  GERGM_Object@hyperparameter_optimization <- hyperparameter_optimization
  GERGM_Object@target_accept_rate <- target_accept_rate
  GERGM_Object@proposal_variance <- proposal_variance
  GERGM_Object@estimation_method <- estimation_method
  GERGM_Object@number_of_simulations <- number_of_networks_to_simulate
  GERGM_Object@thin <- thin
  GERGM_Object@burnin <- MCMC_burnin
  GERGM_Object@MPLE_gain_factor <- MPLE_gain_factor

  #2. Estimate GERGM
  GERGM_Object <- Estimate_GERGM(formula,
                                 MPLE.only = use_MPLE_only,
                                 max.num.iterations = maximum_number_of_lambda_updates,
                                 mc.num.iterations = maximum_number_of_theta_updates,
                                 seed = seed,
                                 tolerance = convergence_tolerance,
                                 possible.stats = possible_structural_terms,
                                 GERGM_Object = GERGM_Object,
                                 force_x_theta_updates = force_x_theta_updates,
                                 verbose = verbose,
                                 force_x_lambda_updates = force_x_lambda_updates)

  #3. Perform degeneracy diagnostics and create GOF plots
  if (!GERGM_Object@theta_estimation_converged) {
    warning("Estimation procedure did not detect convergence in Theta estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
    GERGM_Object <- store_console_output(GERGM_Object,"Estimation procedure did not detect convergence in Theta estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
  }
  if (!GERGM_Object@lambda_estimation_converged) {
    warning("Estimation procedure did not detect convergence in Lambda estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
    GERGM_Object <- store_console_output(GERGM_Object,"Estimation procedure did not detect convergence in Lambda estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
  }

  #now simulate from last update of theta parameters
  GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                 seed1 = seed,
                                 possible.stats = possible_structural_terms)

  colnames(GERGM_Object@lambda.coef) = gpar.names
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))

  # change back column names if we are dealing with an undirected network
  if (!network_is_directed) {
    change <- which(colnames(GERGM_Object@theta.coef) == "in2stars")
    if (length(change) > 0) {
      colnames(GERGM_Object@theta.coef)[change] <- "twostars"
    }
  }

  init.statistics <- NULL
  if (GERGM_Object@is_correlation_network) {
    init.statistics <- h2(GERGM_Object@network,
                          triples = triples,
                          statistics = rep(1, length(possible_structural_terms)),
                          alphas = GERGM_Object@weights,
                          together = downweight_statistics_together)
  }else{
    init.statistics <- h2(GERGM_Object@bounded.network,
                          triples = triples,
                          statistics = rep(1, length(possible_structural_terms)),
                          alphas = GERGM_Object@weights,
                          together = downweight_statistics_together)
  }

  # fix issue with the wrong stats being saved
  GERGM_Object@stats[2,] <- init.statistics
  hsn.tot <- GERGM_Object@MCMC_output$Statistics

  # save these statistics so we can make GOF plots in the future, otherwise
  # they would be the transformed statistics which would produce poor GOF plots.
  GERGM_Object@simulated_statistics_for_GOF <- hsn.tot
  #thin statsitics
  hsn.tot <- Thin_Statistic_Samples(hsn.tot)

  #calculate t.test p-values for calculating the difference in the means of
  # the newly simulated data with the original network
  statistic_test_p_values <- rep(NA,length(possible_structural_terms))
  for (i in 1:length(possible_structural_terms)) {
    statistic_test_p_values[i] <- round(t.test(hsn.tot[, i],
                                      mu = init.statistics[i])$p.value,3)
  }

  stats.data <- data.frame(Observed = init.statistics,
                           Simulated = colMeans(hsn.tot))
  rownames(stats.data) <- possible_structural_terms
  cat("Statistics of observed network and networks simulated from final theta parameter estimates:\n")
  GERGM_Object <- store_console_output(GERGM_Object,"Statistics of observed network and networks simulated from final theta parameter estimates:\n")

  GERGM_Object <- store_console_output(GERGM_Object, toString(stats.data))

  statistic_test_p_values <- data.frame(p_values = statistic_test_p_values)
  rownames(statistic_test_p_values) <- possible_structural_terms
  cat("\nt-test p-values for statistics of observed network and networks simulated from final theta parameter estimates:\n \n")
  GERGM_Object <- store_console_output(GERGM_Object,"\nt-test p-values for statistics of observed network and networks simulated from final theta parameter estimates:\n \n")
  print(statistic_test_p_values)
  GERGM_Object <- store_console_output(GERGM_Object, toString(statistic_test_p_values))


  colnames(statistic_test_p_values) <- "p_values"
  GERGM_Object@observed_simulated_t_test <- statistic_test_p_values

  #test to see if we have an acceptable fit
  acceptable_fit <- statistic_test_p_values[which(GERGM_Object@stats_to_use == 1), 1]

  if (min(acceptable_fit) > acceptable_fit_p_value_threshold) {
    GERGM_Object@acceptable_fit <- TRUE
    message("Parameter estimates simulate networks that are statistically indistinguishable from observed network on the statistics specified by the user. ")
    GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates simulate networks that are statistically indistinguishable from observed network on the statistics specified by the user. ")
  }else{
    GERGM_Object@acceptable_fit <- FALSE
    message("Parameter estimates simulate networks that are statistically distinguishable from observed network. Check GOF plots to determine if the model provides a reasonable fit . This is a very stringent test for goodness of fit, so results may still be acceptable even if this criterion is not met.")
    GERGM_Object <- store_console_output(GERGM_Object, "Parameter estimates simulate networks that are statistically distinguishable from observed network. Check GOF plots to determine if the model provides a reasonable fit . This is a very stringent test for goodness of fit, so results may still be acceptable even if this criterion is not met.")
  }

  #4. output everything to the appropriate files and return GERGM object.
  if (generate_plots) {
    # only generate output if output_name is not NULL
    if (!is.null(output_name)) {
      if (is.null(output_directory)) {
        output_directory <- getwd()
      }
      current_directory <- getwd()
      setwd(output_directory)

      pdf(file = paste(output_name,"_GOF.pdf",sep = ""), height = 4, width = 8)
      GOF(GERGM_Object)
      dev.off()

      pdf(file = paste(output_name,"_Parameter_Estimates.pdf",sep = ""), height = 4, width = 5)
      Estimate_Plot(GERGM_Object)
      dev.off()

      pdf(file = paste(output_name,"_Trace_Plot.pdf",sep = ""), height = 4, width = 6)
      Trace_Plot(GERGM_Object)
      dev.off()

      save(GERGM_Object, file = paste(output_name,"_GERGM_Object.Rdata",sep = ""))

      write.table(GERGM_Object@console_output,file = paste(output_name,"_Estimation_Log.txt",sep = ""),row.names = F,col.names = F,fileEncoding = "utf8", quote = F)

      setwd(current_directory)
    } else{
      # if we are not saving everything to a directory then just print stuff to
      # the graphics device
      GOF(GERGM_Object)
      Sys.sleep(2)
      Estimate_Plot(GERGM_Object)
      Sys.sleep(2)
      Trace_Plot(GERGM_Object)
    }
  }

  # transform networks back to observed scale
  cat("Transforming networks simulated via MCMC as part of the fit diagnostics back on to the scale of observed network. You can access these networks through the '@MCMC_output$Networks' field returned by this function...\n")
  GERGM_Object <- Convert_Simulated_Networks_To_Observed_Scale(GERGM_Object)

  return(GERGM_Object)
}
