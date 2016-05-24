#' @title Active learning reliability method combining Kriging and Monte Carlo
#' Simulation
#'
#' @description Estimate a failure probability with the AKMCS method.
#'
#' @author Clement WALTER \email{clement.walter@cea.fr}
#'
#' @details AKMCS strategy is based on a original Monte-Carlo population which
#' is classified
#' with a kriging-based metamodel. This means that no sampling is done during
#' refinements
#' steps. Indeed, it tries to classify this Monte-Carlo population with a
#' confidence greater
#' than a given value, for instance \sQuote{distance} to the failure should be
#' greater than
#' \code{crit_min} standard deviation.
#'
#' Thus, while this criterion is not verified, the point minimizing it is added to
#' the learning database and then evaluated.
#'
#' Finally, once all points are classified or when the maximum number of calls
#' has been reached, crude Monte-Carlo is performed. A final test controlling
#' the size of this population regarding the targeted coefficient of variation
#' is done; if it is too small then a new population of sufficient size
#' (considering ordre of magnitude of found probability) is generated, and
#' algorithm run again.
#'
#' @note Problem is supposed to be defined in the standard space. If not,
#' use \code{\link{UtoX}} to do so. Furthermore, each time a set of vector
#' is defined as a matrix, \sQuote{nrow} = \code{dimension} and
#' \sQuote{ncol} = number of vector to be consistent with \code{as.matrix}
#' transformation of a vector.
#'
#' Algorithm calls lsf(X) (where X is a matrix as defined previously) and
#' expects a vector in return. This allows the user to optimise the computation
#' of a batch of points, either by vectorial computation, or by the use of
#' external codes (optimised C or C++ codes for example) and/or parallel
#' computation; see examples in \link{MonteCarlo}.
#'
#'
#' @references
#' \itemize{
#' \item
#' B. Echard, N. Gayton, M. Lemaire:\cr
#' \emph{AK-MCS : an Active learning reliability method combining Kriging and
#' Monte Carlo Simulation}\cr
#' Structural Safety, Elsevier, 2011.\cr
#'
#' \item
#' B. Echard, N. Gayton, M. Lemaire and N. Relun:\cr
#' \emph{A combined Importance Sampling and Kriging reliability method for
#' small failure probabilities with time-demanding numerical models}\cr
#' Reliability Engineering \& System Safety,2012\cr
#'
#' \item
#' B. Echard, N. Gayton and A. Bignonnet:\cr
#' \emph{A reliability analysis method for fatigue design}\cr
#' International Journal of Fatigue, 2014\cr
#' }
#'
#' @seealso
#' \code{\link{SubsetSimulation}}
#' \code{\link{MonteCarlo}}
#' \code{\link{MetaIS}}
#' \code{\link[DiceKriging]{km}} (in package \pkg{DiceKriging})
#'
#' @examples
#' \dontrun{
#' res = AKMCS(dimension=2,lsf=kiureghian,plot=TRUE)
#'
#' #Compare with crude Monte-Carlo reference value
#' N = 500000
#' dimension = 2
#' U = matrix(rnorm(dimension*N),dimension,N)
#' G = kiureghian(U)
#' P = mean(G<0)
#' cov = sqrt((1-P)/(N*P))
#' }
#'
#' #See impact of kernel choice with serial function from Waarts:
#' waarts = function(u) {
#'   u = as.matrix(u)
#'   b1 = 3+(u[1,]-u[2,])^2/10 - sign(u[1,] + u[2,])*(u[1,]+u[2,])/sqrt(2)
#'   b2 = sign(u[2,]-u[1,])*(u[1,]-u[2,])+7/sqrt(2)
#'   val = apply(cbind(b1, b2), 1, min)
#' }
#'
#' \dontrun{
#' res = list()
#' res$matern5_2 = AKMCS(2, waarts, plot=TRUE)
#' res$matern3_2 = AKMCS(2, waarts, kernel="matern3_2", plot=TRUE)
#' res$gaussian  = AKMCS(2, waarts, kernel="gauss", plot=TRUE)
#' res$exp       = AKMCS(2, waarts, kernel="exp", plot=TRUE)
#' 
#' #Compare with crude Monte-Carlo reference value
#' N = 500000
#' dimension = 2
#' U = matrix(rnorm(dimension*N),dimension,N)
#' G = waarts(U)
#' P = mean(G<0)
#' cov = sqrt((1-P)/(N*P))
#' }
#'
#' @import ggplot2
#' @import DiceKriging
#' @importFrom grDevices dev.off gray pdf rainbow
#' @importFrom graphics hist lines points
#' @importFrom stats chisq.test kmeans optimize pchisq pnorm ppois qbeta qgamma qnorm quantile qunif qweibull rnorm runif sd var
#' @importFrom utils capture.output flush.console
#' @export

AKMCS = function(dimension,
                 #' @param dimension dimension of the input space.
                 lsf,
                 #' @param lsf the function defining the failure/safety domain.
                 N    = 500000,
                 #' @param N Monte-Carlo population size.
                 N1   = 10*dimension,
                 #' @param N1 size of the first DOE.
                 Nmax = 200,
                 #' @param Nmax maximum number of calls to the LSF.
                 learn_db  = NULL,
                 #' @param learn_db coordinates of already known points.
                 lsf_value = NULL,
                 #' @param lsf_value value of the LSF on these points.
                 failure   = 0,
                 #' @param failure failure threshold.
                 precision = 0.05,
                 #' @param precision maximum desired cov on the Monte-Carlo estimate.
                 bayesian = TRUE,
                 #' @param bayesian estimate the conditional expectation E_X [ P[meta(X)<failure] ].
                 meta_model = NULL,
                 #' @param meta_model provide here a kriging metamodel from km if wanted.
                 kernel = "matern5_2",
                 #' @param kernel specify the kernel to use for km.
                 learn_each_train = TRUE,
                 #' @param learn_each_train specify if kernel parameters are re-estimated at each train.
                 crit_min = 2,
                 #' @param crit_min minimum value of the criteria to be used for refinement.
                 lower.tail = TRUE,
                 #' @param lower.tail as for pxxxx functions, TRUE for estimating P(lsf(X) < failure), FALSE
                 #' for P(lsf(X) > failure)
                 limit_fun_MH = NULL,
                 #' @param limit_fun_MH define an area of exclusion with a limit function.
                 failure_MH = 0,
                 #' @param failure_MH the theshold for the limit_fun_MH function.
                 sampling_strategy = "MH",
                 #' @param sampling_strategy either MH for Metropolis-Hastings of AR for accept-reject.
                 first_DOE = "Gaussian",
                 #' @param first_DOE either Gaussian or Uniform, to specify the population on which
                 #' clustering is done.
                 seeds = NULL,
                 #' @param seeds if some points are already known to be in the appropriate subdomain.
                 seeds_eval = limit_fun_MH(seeds),
                 #' @param seeds_eval value of the metamodel on these points.
                 burnin = 30,
                 #' @param burnin burnin parameter for MH.
                 plot = FALSE,
                 #' @param plot set to TRUE for a full plot, ie refresh at each iteration.
                 limited_plot = FALSE,
                 #' @param limited_plot set to TRUE for a final plot with final DOE, metamodel and LSF.
                 add = FALSE,
                 #' @param add if plots are to be added to a current device.
                 output_dir = NULL,
                 #' @param output_dir if plots are to be saved in jpeg in a given directory.
                 verbose = 0) {
                 #' @param verbose either 0 for almost no output, 1 for medium size output and 2 for all
                 #' outputs.


cat("==========================================================================================\n")
cat("                              Beginning of AK-MCS algorithm\n")
cat("==========================================================================================\n\n")

# Fix NOTE issue with R CMD check
x <- y <- z <- ..level.. <- crit <- NULL

## STEP 0 : INITIALISATION

Ncall  = 0
cov    = Inf
Nfailure = 0

if(lower.tail==FALSE){
  lsf_dec = lsf
  lsf = function(x) -1*lsf_dec(x)
  failure <- -failure
}

# plotting part
if(plot==TRUE){

  if(!is.null(output_dir)) {
    fileDir = paste(output_dir,"_AKMCS.pdf",sep="")
    pdf(fileDir)
  }
  xplot <- yplot <- c(-80:80)/10
  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
  p <- ggplot2::ggplot(data = df_plot, aes(x,y)) +
    ggplot2::geom_contour(aes(z=z, color=..level..), breaks = failure) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlim(-8, 8) + ggplot2::ylim(-8, 8)

  if(!is.null(learn_db)){
    row.names(learn_db) <- c('x', 'y')
    p <- p + ggplot2::geom_point(data = data.frame(t(learn_db), z = lsf_value), ggplot2::aes(color=z))
  }
  if(!is.null(limit_fun_MH)) {
    df_plot_MH = data.frame(expand.grid(x=xplot, y=yplot), z = limit_fun_MH(t(expand.grid(x=xplot, y=yplot))))
    p <- p + ggplot2::geom_contour(data = df_plot_MH, aes(z=z, color=..level..), breaks = failure_MH)
  }
  print(p)
}

# while(cov>precision){

if(Nfailure==0){
	cat(" ============================================= \n")
	cat(" STEP 1 : GENERATION OF THE WORKING POPULATION \n")
	cat(" ============================================= \n\n")


	if(is.null(limit_fun_MH)){
		if(verbose>0){cat(" * Generate N =",N,"standard Gaussian samples\n\n")}
		U = matrix(rnorm(dimension*N),dimension,N)
	}
	else{
	  if(sampling_strategy=="MH"){
	    if(verbose>0){cat(" * Generate N =",N,"points from",dim(seeds)[2],"seeds with Metropolis-Hastings algorithm\n")}

	    K = function(x){
	      W = array(rnorm(x), dim = dim(x))
	      sigma = apply(x, 1, sd)*2
	      y = (x + sigma*W)/sqrt(1 + sigma^2)
	      return(y)
	    }

	    seed <- sample.int(n = length(seeds_eval), size = N, replace = TRUE)
	    U <- seeds[,seed]
	    U_eval <- seeds_eval[seed]

	    replicate(burnin, {
	      U_star <- K(U)
	      U_eval_star <- limit_fun_MH(U_star)
	      indU_star <- U_eval_star<failure_MH
	      U[,indU_star] <<- U_star[,indU_star]
	      U_eval[indU_star] <<- U_eval_star[indU_star]
	    })
	    rm(U_eval)
	  }
	  else{
	    if(verbose>0){cat(" * Generate the N =",N,"Monte-Carlo population with an accept-reject strategy on a standard Gaussian sampling\n")}
	    rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}
	    U = generateWithAR(dimension=dimension,
	                       N=N,
	                       limit_f=limit_fun_MH,
	                       failure=failure_MH,
	                       rand=rand)
	  }
	}

	cat(" ================== \n")
	cat(" STEP 2 : FIRST DoE \n")
	cat(" ================== \n\n")

	switch(first_DOE,
		Gaussian = {
		  if(verbose>0){cat(" * Get N1 =",N1,"points by clustering of the N =",N,"points\n")}
			DoE = t(kmeans(t(U), centers=N1,iter.max=20)$centers)
		},
		Uniform = {
		  if(verbose>0){cat(" * Get N1 =",N1,"points with a uniform sampling in a hyper sphere of radius max(radius(N points)) \n")}
			radius <- max(sqrt(rep(1,dim(U)[1])%*%U^2))
			DoE = t( kmeans( runifSphere(dimension,N,radius), centers=N1, iter.max=20)$centers )
		},
		No = {
		  if(verbose>0){cat(" * No first DoE requested \n")}
		},
		stop("Wrong first DOE sampling strategy\n")
	)

	if(first_DOE!="No"){
	  if(verbose>0){cat(" * Add points to the learning database\n")}
	  lsf_DoE <- lsf(DoE);Ncall = Ncall + N1
	  if(is.null(learn_db)){
	    learn_db = array(DoE, dim = dim(DoE), dimnames = list(rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]))
	    lsf_value = lsf_DoE
	  }
	  else{
	    learn_db = cbind(learn_db,DoE)
	    lsf_value = c(lsf_value,lsf_DoE)
	  }
	}

	if(plot==TRUE){
	  if(verbose>0){cat(" * 2D PLOT : First DoE \n")}
	    p <- p + geom_point(data = data.frame(t(learn_db), z = lsf_value), aes(color=z))
	    print(p)
	}
}

if(verbose>0){cat(" * Train the model :\n")}
	if(is.null(meta_model) || learn_each_train==TRUE || Nfailure>0) {
	  if(verbose>1){cat("    - Learn hyperparameters !!! \n")}
    meta = trainModel(design   = learn_db,
						  response = lsf_value,
						  kernel   = kernel,
						  type="Kriging")
    Nfailure = 0
	}
  else {
    if(verbose>1){cat("    - Use previous hyperparameters !!! \n")}
    meta = trainModel(meta_model,
                      updesign=DoE,
                      upresponse=lsf_DoE,
                      type="Kriging")
  }

  meta_model = meta$model
  meta_fun   = meta$fun

  if(verbose>0){cat(" * Evaluate criterion on the work population\n")}
  meta_pred = meta_fun(U)
  criterion <- abs(failure - meta_pred$mean)/meta_pred$sd
  minC = min(criterion)

  if(verbose>1){cat("    - minimum value of the criterion =",minC,"\n\n")}

  #plotting part
  if(plot == TRUE){
    if(verbose>0){cat(" * 2D PLOT : FIRST APPROXIMATED LSF USING KRIGING \n")}
    z_meta = meta_fun(t(df_plot[,1:2]))
    df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(failure - z_meta$mean)/z_meta$sd)
    print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = failure) +
      geom_contour(data = df_plot_meta, aes(z=crit, color=..level.., alpha = 0.5), linetype = 4, breaks = crit_min))
  }

	cat(" ======================= \n")
	cat(" STEP 3 : UPDATE THE DoE \n")
	cat(" ======================= \n")

	k = 0;
	while(minC<crit_min & k<Nmax) {

		k = k+1;
		if(verbose>0){cat("\n * ITERATION ",k,"\n")}
		if(verbose>1){cat("   - min < 2 & k =",k,"< Nmax =",Nmax,"=> improve the model\n")}

		candidate = as.matrix(U[,which.min(criterion)])
		eval      = lsf(candidate);Ncall = Ncall + 1
		learn_db  = cbind(learn_db,candidate)
		lsf_value = c(lsf_value, eval)

		#plotting part
		if(plot==TRUE){
		  if(verbose>1){cat("   - 2D PLOT : UPDATE \n")}
		  p <- p + geom_point(data = data.frame(x=candidate[1], y=candidate[2], z = eval), aes(color = z))
			print(p_meta + geom_point(data = data.frame(x=candidate[1], y=candidate[2]), color = "red", size = 4))
		}

		if(verbose>1){cat("   - Train the model\n")}
		if(learn_each_train==TRUE) {
		  if(verbose>1){cat("   - Learn hyperparameters !!! \n")}
			meta = trainModel(design=learn_db,
							  response=lsf_value,
							  kernel=kernel,
							  type="Kriging")
      Nfailure = 0
		}
		else {
		  if(verbose>1){cat("     + Use previous hyperparameters !!! \n")}
			meta = trainModel(meta_model,
							  updesign=candidate,
							  upresponse=eval,
							  type="Kriging")
		}
		meta_model = meta$model
		meta_fun   = meta$fun

		if(verbose>0){cat(" * Evaluate criterion on the work population\n")}
		meta_pred = meta_fun(U)
		criterion <- abs(failure - meta_pred$mean)/meta_pred$sd
		minC = min(criterion)
		cat("     + minimum value of the criterion =",minC,"\n")

		#plotting part
		if(plot==TRUE){
		  if(verbose>0){cat("   - 2D PLOT : END ITERATION",k,"\n")}
		  z_meta = meta_fun(t(df_plot[,1:2]))
		  df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(failure - z_meta$mean)/z_meta$sd)
		  print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = failure) +
		          geom_contour(data = df_plot_meta, aes(z=crit, color=..level.., alpha = 0.5), linetype = 4, breaks = crit_min))
		}
	}

	cat(" ======================================================================================= \n")
	cat(" STEP 4 : EVALUATE FAILURE PROBABILITY WITH A MONTE-CARLO ESTIMATOR USING THE META-MODEL\n")
	cat(" ======================================================================================= \n")

	if(bayesian==TRUE){
	  P = mean(pnorm((failure-meta_pred$mean)/meta_pred$sd))
	}
	else{
	  P = mean(meta_pred$mean<failure)
	}
	cov = sqrt((1-P)/(N*P))

cat("   - p =",P,"\n")
cat("   - failure =", failure,"\n")
cat("   - 95% confidence interval on Monte Carlo estimate:",P*(1-2*cov),"< p <",P*(1+2*cov),"\n")
cat(" * cov =",cov,"\n")

	if( cov > precision) {
	  Nfailure = sum(lsf_value<failure)
		if(P>0){
			N = ceiling((1-P)/(precision^2*P))
			cat("   => cov too large ; this order of magnitude for the probability brings N =",N,"\n")
		}
		else {
			cat("   => cov too large, only",Nfailure,"failings points in the learn_db\n")
		}
	}
  else{
    cat(" * cov < precision =",precision,"; End of AKMCS algorithm\n")
    cat("   - Pf =",P,"\n")
    cat("   - cov =",cov,"\n")
    cat("   - Ncall =",Ncall,"\n")
  }
# }

#plotting part
if(limited_plot==TRUE){
cat("\n * 2D PLOT : LSF, FINAL DATABASE AND METAMODEL")
  z_meta = meta_fun(t(df_plot[,1:2]))
  df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(failure - z_meta$mean)/z_meta$sd)
  print(p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = failure) +
          geom_contour(data = df_plot_meta, aes(z=crit, color=..level.., alpha = 0.5), linetype = 4, breaks = crit_min) +
          geom_point(data = data.frame(t(learn_db), z = lsf_value), aes(color=z)))
}

if( (plot || limited_plot) & (add==FALSE) & !is.null(output_dir) ) { dev.off() }

points = U[,(meta_pred$mean<0)]

#' @return An object of class \code{list} containing the failure probability and some
#' more outputs as described below:
res = list(p=P,
           #' \item{p}{the estimated failure probability.}
           cov=cov,
           #' \item{cov}{the coefficient of variation of the Monte-Carlo probability estimate.}
           Ncall=Ncall,
           #' \item{Ncall}{the total number of calls to the \code{lsf}.}
           learn_db=learn_db,
           #' \item{learn_db}{the final learning database, ie. all points where \code{lsf} has
           #' been calculated.}
           lsf_value=(-1)^(!lower.tail)*lsf_value,
           #' \item{lsf_value}{the value of the \code{lsf} on the learning database.}
           meta_fun=function(x) {
             g = meta_fun(x)
             g$mean = (-1)^(!lower.tail)*g$mean
             return(g)
           },
           #' \item{meta_fun}{the metamodel approximation of the \code{lsf}. A call output is a
           #' list containing the value and the standard deviation.}
           meta_model=meta_model,
           #' \item{meta_model}{the final metamodel. An S4 object from \pkg{DiceKriging}. Note
           #' that the algorithm enforces the problem to be the estimation of P[lsf(X)<failure]
           #' and so using \sQuote{predict} with this object will return inverse values if
           #' \code{lower.tail==FALSE}; in this scope prefer using directly \code{meta_fun} which
           #' handles this possible issue.}
           points=points,
           #' \item{points}{points in the failure domain according to the metamodel.}
           meta_eval=meta_fun(points))
           #' \item{meta_eval}{evaluation of the metamodel on these points.}
if(plot+limited_plot) {
  res = c(res, list(z_meta=z_meta$mean))
          #' \item{z_meta}{if \code{plot}==TRUE, the evaluation of the metamodel on the plot grid.}
}
return(res)
}
