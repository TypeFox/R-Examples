#' For each set of parameters, the simulation is ran once to obtain the value of estimator and true value to make sure everything in ezsim is properly defined. The test results will be shown in the console. The test will be ran automatically when you create an ezsim object.
#' @name test.ezsim
#' @aliases test.ezsim
#' @title Perform a Test for an ezsim Object
#' @method test ezsim
#' @param x An ezsim Object
#' @param return_name Whehter to return the name of estimator
#' @param print_result Whehter to print the return
#' @param \dots unused
#' @return Optional: names of estimator.
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \code{\link{ezsim}}
#' @export 
#' @examples         
#' \dontrun{ 
#' ezsim_basic<-ezsim(
#'     m             = 100,
#'     run           = FALSE,
#'     run_test      = FALSE,
#'     display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
#'     parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
#'     dgp           = function() rnorm(n,mu,sigma),
#'     estimator     = function(x) c(mean_hat = mean(x), 
#'                                  sd_mean_hat=sd(x)/sqrt(length(x)-1)),
#'     true_value    = function() c(mu, sigma / sqrt(n-1))
#' )
#' 
#' test(ezsim_basic,print_result=TRUE)
#' }
test.ezsim <-
function(x,return_name=TRUE,print_result=FALSE,...){
	parameter_list <- generate(x$parameter_def)
	
	## Decide whether we need to stop the cluster after the test
	create_cluster_flag <- FALSE
	
	if (x$use_core > 1 & is.null(x$cluster) ){
		x$cluster <- makeCluster(x$use_core)
		if (!is.null(x$cluster_packages)){
			for (i in 1:length(x$cluster_packages))
				eval(substitute( 
					clusterEvalQ(x$cluster , require(w, character.only=TRUE)   )  ,
					list( w=x$cluster_packages[i] )
				))
		}
		clusterSetRNGStream(x$cluster,x$use_seed)
		create_cluster_flag<-TRUE
	}

	if (x$use_core == 1) 
		set.seed(x$use_seed)

	i=NULL
	tryCatch({
		## test for estimator
		cat("Testing for estimator...")
		compute_estimates <- function(i,ezsim_object) {
			ezsim_object$estimator(Jmisc::evalFunctionOnList(ezsim_object$dgp,i))
		}
			
		test_estimates<-
		if (!is.null(x$cluster)){
			parLapply(x$cluster,parameter_list,fun=compute_estimates,ezsim_object=x)
		} else{
			lapply(parameter_list,FUN=compute_estimates,ezsim_object=x)
		}
		cat("Passed\n")
		
		## test for estimator parser
		cat("Testing for estimator parser...")
		test_estimates_parsed<-
			foreach( i = test_estimates ) %do%{
				x$estimator_parser(i)
			}		
		if (any(!sapply(test_estimates_parsed,is.vector)))
			stop("estimator parser do not return a vector")

		length_estimates<-sapply(test_estimates_parsed,length)	
		if (!all(length_estimates==length_estimates[[1]]))
			stop("length of estimates are not the same")
		cat("Passed\n")
			
		# Check true value	
		cat("Testing for true value...")
		if (is.function(x$true_value)){
			compute_true_value <- function(i,ezsim_object) {
				Jmisc::evalFunctionOnList(ezsim_object$true_value,i)
			}

			test_true_value<-
			if (!is.null(x$cluster)){
				parLapply(x$cluster,parameter_list,fun=compute_true_value,ezsim_object=x)
			} else{
				lapply(parameter_list,FUN=compute_true_value,ezsim_object=x)
			}
			
			## check length
			length_true_value<-sapply(test_true_value,length)
			if (!all(length_true_value==length_true_value[[1]]))
				stop("length of true value are not the same")
			
			if (length_true_value[[1]]!=length_estimates[[1]])
				stop("length of true value and estimates are not the same")
		}
		cat("Passed\n")

	}, finally = {
		if (create_cluster_flag){
			tryCatch({
				stopCluster(x$cluster)
			}, finally = {
				x$cluster<-NULL
			})
		}
	})
}
