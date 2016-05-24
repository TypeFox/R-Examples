#' @title Metamodel based Impotance Sampling
#' 
#' @description Estimate failure probability by MetaIS method.
#' 
#' @author Clement WALTER \email{clement.walter@cea.fr}
#' 
#' @details
#' MetaIS is an Important Sampling based probability estimator. It makes use of
#' a kriging surogate to approximate the optimal density function, replacing the
#' indicatrice by its kriging pendant, the probability of being in the failure
#' domain. In this context, the normallizing constant of this quasi-optimal PDF
#' is called the \sQuote{augmented failure probability} and the modified
#' probability \sQuote{alpha}.
#' 
#' After a first uniform Design of Experiments, MetaIS uses an alpha
#' Leave-One-Out criterion combined with a margin sampling strategy to refine
#' a kriging-based metamodel. Samples are generated according to the weighted
#' margin probability with Metropolis-Hastings algorithm and some are selected
#' by clustering; the \code{N_seeds} are got from an accept-reject strategy on
#' a standard population.
#' 
#' Once criterion is reached or maximum number of call done, the augmented
#' failure probability is estimated with a crude Monte-Carlo. Then, a new
#' population is generated according to the quasi-optimal instrumenal PDF;
#' \code{burnin} and \code{thinning} are used here and alpha is evaluated.
#' While the coefficient of variation of alpha estimate is greater than a
#' given threshold and some computation spots still available (defined by
#' \code{Ncall_max}) the estimate is refined with extra calculus.
#' 
#' The final probability is the product of p_epsilon and alpha, and final
#' squared coefficient of variation is the sum of p_epsilon and alpha one's.
#' 
#' @return   An object of class \code{list} containing the failure probability
#' and some more outputs as described below:
#' \item{p}{The estimated failure probability.}
#' \item{cov}{The coefficient of variation of the Monte-Carlo probability
#' estimate.}
#' \item{Ncall}{The total number of calls to the \code{lsf}.}
#' \item{learn_db}{The final learning database, ie. all points where \code{lsf}
#' has been calculated.}
#' \item{lsf_value}{The value of the \code{lsf} on the learning database.}
#' \item{meta_fun}{The metamodel approximation of the \code{lsf}. A call output
#' is a list containing the value and the standard deviation.}
#' \item{meta_model}{The final metamodel. An S4 object from \pkg{DiceKriging}.
#' Note that the algorithm enforces the problem to be the estimation of
#' P[lsf(X)<failure] and so using \sQuote{predict} with this object will
#' return inverse values if \code{lower.tail==FALSE}; in this scope prefer
#' using directly \code{meta_fun} which handle this possible issue.}
#' \item{points}{Points in the failure domain according to the metamodel.}
#' \item{meta_eval}{Evaluation of the metamodel on these points.}
#' \item{z_meta}{If \code{plot}==TRUE, the evaluation of the metamodel on
#' the plot grid.}
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
#' @references
#' \itemize{
#'   \item
#'   V. Dubourg:\cr
#'   Meta-modeles adaptatifs pour l'analyse de fiabilite et l'optimisation sous
#'   containte fiabiliste\cr
#'   PhD Thesis, Universite Blaise Pascal - Clermont II,2011\cr
#'   
#'   \item
#'   V. Dubourg, B. Sudret, F. Deheeger:\cr
#'   Metamodel-based importance sampling for structural reliability analysis
#'   Original Research Article\cr
#'   Probabilistic Engineering Mechanics, Volume 33, July 2013, Pages 47-57\cr
#'   
#'   \item
#'   V. Dubourg, B. Sudret:\cr
#'   Metamodel-based importance sampling for reliability sensitivity analysis.\cr
#'   Accepted for publication in Structural Safety, special issue in the honor
#'   of Prof. Wilson Tang.(2013)\cr
#'   
#'   \item
#'   V. Dubourg, B. Sudret and J.-M. Bourinet:\cr
#'   Reliability-based design optimization using kriging surrogates and subset
#'   simulation.\cr
#'   Struct. Multidisc. Optim.(2011)\cr
#' }
#' 
#' @seealso
#' \code{\link{SubsetSimulation}}
#' \code{\link{MonteCarlo}}
#' \code{\link[DiceKriging]{km}} (in package \pkg{DiceKriging})
#' 
#' @examples 
#' kiureghian = function(x, b=5, kappa=0.5, e=0.1) {
#' x = as.matrix(x)
#' b - x[2,] - kappa*(x[1,]-e)^2
#' }
#' 
#' \dontrun{
#' res = MetaIS(dimension=2,lsf=kiureghian,plot=TRUE)
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
#' #See impact of kernel choice with Waarts function :
#' waarts = function(u) {
#'   u = as.matrix(u)
#'   b1 = 3+(u[1,]-u[2,])^2/10 - sign(u[1,] + u[2,])*(u[1,]+u[2,])/sqrt(2)
#'   b2 = sign(u[2,]-u[1,])*(u[1,]-u[2,])+7/sqrt(2)
#'   val = apply(cbind(b1, b2), 1, min)
#' }
#' 
#' \dontrun{
#' res = list()
#' res$matern5_2 = MetaIS(2,waarts,plot=TRUE)
#' res$matern3_2 = MetaIS(2,waarts,kernel="matern3_2",plot=TRUE)
#' res$gaussian = MetaIS(2,waarts,kernel="gauss",plot=TRUE)
#' res$exp = MetaIS(2,waarts,kernel="exp",plot=TRUE)
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
#' @import Matrix
#' @import mvtnorm
#' @export

MetaIS = function(dimension,
                  #' @param dimension of the input space
                  lsf,
                  #' @param lsf the failure defining the failure/safety domain
            			## Algorithm parameters
            			N = 500000,
            			#' @param N size of the Monte-Carlo population for P_epsilon estimate
            			N_alpha = 100,
            			#' @param N_alpha initial size of the Monte-Carlo population for alpha estimate
            			N_DOE = 10*dimension,
            			#' @param N_DOE size of the initial DOE got by clustering of the N1 samples
            			N1 = N_DOE*30,
            			#' @param N1 size of the initial uniform population sampled in a hypersphere of radius Ru
            			Ru = 8,
            			#' @param Ru radius of the hypersphere for the initial sampling
            			Nmin = 30, 
            			#' @param Nmin minimum number of call for the construction step
            			Nmax = 200,          
            			#' @param Nmax maximum number of call for the construction step
            			Ncall_max = 1000,     
            			#' @param Ncall_max maximum number of call for the whole algorithm
            			precision = 0.05,      
            			#' @param precision desired maximal value of cov
            			N_seeds = 2*dimension, 
            			#' @param N_seeds number of seeds for MH algoritm while generating into the margin (
            			#' according to MP*gauss)
            			Niter_seed = 10000,     
            			#' @param Niter_seed maximum number of iteration for the research of a seed for alphaLOO
            			#' refinement sampling
            			N_alphaLOO = 5000,      
            			#' @param N_alphaLOO number of points to sample at each refinement step
            			K_alphaLOO = 2*dimension, 
            			#' @param K_alphaLOO number of clusters at each refinement step
            			alpha_int = c(0.1,10),  
            			#' @param alpha_int range for alpha to stop construction step
            			k_margin = 1.96,        
            			#' @param k_margin margin width; default value means that points are classified with more
            			#' than 97,5\%
            			lower.tail = TRUE,      
            			#' @param lower.tail specify if one wants to estimate P[lsf(X)<failure] or P[lsf(X)>failure].
            			## Subset parameters
                  learn_db  = NULL,       
            			#' @param learn_db Coordinates of alredy known points
                  lsf_value = NULL,       
            			#' @param lsf_value Value of the LSF on these points
                  failure   = 0,          
            			#' @param failure Failure threshold
                  meta_model = NULL,      
            			#' @param meta_model Provide here a kriging metamodel from km if wanted     
                  kernel = "matern5_2",   
            			#' @param kernel Specify the kernel to use for km
                  learn_each_train = TRUE,
            			#' @param learn_each_train Specify if kernel parameters are re-estimated at each train
                  limit_fun_MH = NULL,    
            			#' @param limit_fun_MH Define an area of exclusion with a limit function
            			failure_MH = 0,         
            			#' @param failure_MH Threshold for the limit_MH function
                  sampling_strategy = "MH",
            			#' @param sampling_strategy Either MH for Metropolis-Hastings of AR for accept-reject
                  seeds = NULL,           
            			#' @param seeds If some points are already known to be in the appropriate subdomain
                  seeds_eval = limit_fun_MH(seeds), 
            			#' @param seeds_eval Value of the metamodel on these points
                  burnin = 20,            
            			#' @param burnin Burnin parameter for MH
                  ## plot parameter
                  plot = FALSE,           
            			#' @param plot Set to TRUE for a full plot, ie refresh at each iteration
                  limited_plot = FALSE,   
            			#' @param limited_plot Set to TRUE for a final plot with final DOE, metamodel and LSF
                  add = FALSE,            
            			#' @param add If plots are to be added to a current device
                  output_dir = NULL,      
            			#' @param output_dir If plots are to be saved in jpeg in a given directory
                  verbose = 0) {          
                  #' @param verbose Either 0 for almost no output, or 1 for medium size or 2 for all outputs


cat("==========================================================================================\n")
cat("                              Beginning of Meta-IS algorithm \n")
cat("==========================================================================================\n\n")

cat("===========================================================================\n")
cat(" STEP 1 : Adaptative construction of h the approximated optimal density \n")
cat("===========================================================================\n\n")

## Init
Ncall = 0
cov_epsilon = Inf
ITER = 0;

# Fix NOTE issue with R CMD check
x <- z <- ..level.. <- crit <- NULL

if(lower.tail==FALSE){
  lsf_dec = lsf
  lsf = function(x) -1*lsf_dec(x)
  failure <- -failure
}

# plotting part
if(plot==TRUE){
  
  if(!is.null(output_dir)) {
    fileDir = paste(output_dir,"_Meta_IS.pdf",sep="")
    pdf(fileDir)
  }
  xplot <- yplot <- seq(-Ru, Ru, l = 20*Ru)
  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
  p <- ggplot(data = df_plot, aes(x,y)) +
    geom_contour(aes(z=z, colour = ..level..), breaks = failure) +
    theme(legend.position = "none") +
    xlim(-8, 8) + ylim(-8, 8)
  
  if(!is.null(learn_db)){
    row.names(learn_db) <- rep(c('x', 'y'), length.out= dimension)
    p <- p + geom_point(data = data.frame(t(learn_db), z = lsf_value), aes(color=z))
  }
  if(!is.null(limit_fun_MH)) {
    df_plot_MH = data.frame(expand.grid(x=xplot, y=yplot), z = limit_fun_MH(t(expand.grid(x=xplot, y=yplot))))
    p <- p + geom_contour(data = df_plot_MH, aes(z=z, colour=..level..), breaks = failure_MH)
  }
  print(p)
}
	
# while(cov_epsilon > precision){

	cat("\n A- REFINEMENT OF PROBABILISTIC CLASSIFICATION FUNCTION PI \n")
	cat("    ======================================================== \n\n")

	ITER = ITER + 1
	cat(" ITERATION ",ITER,"\n")
	cat(" -------------\n\n")

	if(N_DOE>0){
	  if(verbose>0){cat(" * Generate N1 =",N1,"samples uniformly distributed in a hypersphere of radius Ru =",Ru,"\n")}
		U = runifSphere(dimension,N1,radius=Ru)

	  if(verbose>0){cat(" * Get N_DOE =",N_DOE,"points by clustering of the N1 =",N1,"points\n")}
		if(!is.null(limit_fun_MH)) {
			ind = limit_fun_MH(U$N1) #in a Subset algorithm, select points in Fi-1
			prop = sum(ind<failure_MH)/N1 #get accept-reject proportion
			U = runifSphere(dimension,ceiling(N1/prop),radius=Ru) #generate more points to take AR proportion into account
			ind = limit_fun_MH(U)
			U = U[,ind<failure_MH] #finally get ~N1 points uniformly distributed in the subdomain Fi-1
		}
		DoE = t(kmeans(U, centers=N_DOE,iter.max=20)$centers)
		rm(U)
		
	  if(verbose>0){cat(" * Assessment of performance function G on these points\n")}
		lsf_DoE = lsf(DoE);Ncall = Ncall + N_DOE
	
	  if(verbose>0){cat(" * Add points to the learning database\n")}
		if(is.null(learn_db)){
		  learn_db = array(DoE, dim = dim(DoE), dimnames = list(rep(c('x', 'y'), length.out = dimension)))
			lsf_value = lsf_DoE
		}
		else{
			learn_db = cbind(learn_db,DoE)
			lsf_value = c(lsf_value, lsf_DoE)
		}

		if(plot==TRUE){
		  if(verbose>0){cat(" * 2D PLOT : First DoE \n")}
		  p <- p + geom_point(data = data.frame(t(learn_db), z = lsf_value), aes(color=z))
		  print(p)
		}
	}

	if(verbose>0){cat(" * Train the model\n")}
	if(is.null(meta_model) || learn_each_train==TRUE) {
	  if(verbose>1){cat("   - Learn hyperparameters !!! \n")}
		meta = trainModel(design=learn_db,
		      response=lsf_value,
		      kernel=kernel,
		      type="Kriging")
	}
	else {
	  if(verbose>1){cat("   - Use previous hyperparameters !!! \n")}
		meta = trainModel(meta_model,
		      updesign=DoE,
		      upresponse=lsf_DoE,
		      type="Kriging")
	}
	
	if(verbose>0){cat("\n * UPDATE quantities based on kriging surrogate model : MP, wMP, pi\n")}
	meta_model = meta$model
	meta_fun = meta$fun
	MP = function(x,k=k_margin) {
		x = as.matrix(x)
		G_meta = meta_fun(x)
		res = pnorm((failure + k*G_meta$sd - G_meta$mean)/G_meta$sd) - pnorm((failure - k*G_meta$sd - G_meta$mean)/G_meta$sd)
		return(res)
	}
	wMP = function(x) {
		x = as.matrix(x)
		MP(x)*exp(-0.5*rep(1,dim(x)[1])%*%x^2)
	}
	pi = function(x) {
		x = as.matrix(x)
		G_meta = meta_fun(x)
		pnorm((failure-G_meta$mean)/G_meta$sd)
	}
	
	#plotting part
	if(plot == TRUE){
	  if(verbose>0){cat(" * 2D PLOT : FIRST APPROXIMATED LSF USING KRIGING \n")}
	  z_meta = meta_fun(t(df_plot[,1:2]))
	  df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(failure - z_meta$mean)/z_meta$sd)
	  print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level..), breaks = failure) +
	          geom_contour(data = df_plot_meta, aes(z=crit, color=..level..), linetype = 4, breaks = k_margin))
	}

	k = 0;
	Nmax = min(Nmax,Ncall_max - N_alpha);
	if(verbose>0){cat(" * Calculate alphaLOO \n")}
	LOO = leaveOneOut.km(meta_model, type="UK", trend.reestim=FALSE)
	piLOO = pnorm((failure-LOO$mean)/LOO$sd)
	notNull = piLOO>(10^(-16))
	if(verbose>1){cat("   -",sum(piLOO<10^(-16)),"samples not considered as pi<10^-16\n")}
	alphaLOO = mean(1*(lsf_value[notNull]<failure)/piLOO[notNull])
	if(verbose>1){cat("   - alphaLOO =",alphaLOO,"\n")}
	criterion = (alpha_int[1]<alphaLOO)*(alphaLOO<alpha_int[2])*(k>=Nmin) + (k>Nmax)

	while(!criterion) {
	  if(verbose>0){cat(" * Criterion not reached :\n")}
	  if(verbose>1){cat("   - alphaLOO :",alpha_int[1],"<",alphaLOO,"<",alpha_int[2],"\n")}
	  if(verbose>1){cat("   - k :",Nmin,"<",k,"<",Nmax,"\n")}
		ITER = ITER + 1
		cat("\n\n ITERATION ",ITER,"\n")
		cat(" -------------\n\n")
	  if(verbose>0){cat(" * Find seeds using accept-reject strategy on a standard population\n")}
		n     = 0;
		niter = 0;
		candidate = matrix(NA,dimension,N_seeds, dimnames = list(rep(c('x', 'y'), length.out = dimension)))
		while(n<N_seeds && niter<=Niter_seed) {
			niter = niter + 1
			missing_seeds = N_seeds - n
			tmp = matrix(rnorm(dimension*missing_seeds),dimension,missing_seeds)
			is_margin = wMP(tmp)
			is_margin = is_margin>(10^(-16))
			found_seeds = sum(is_margin)
			if(found_seeds>0) candidate[,n + 1:found_seeds] = tmp[,is_margin]
			n = sum(!is.na(candidate[1,])) #number of NAs stands for number of missing vectors
		}
		if(niter>Niter_seed && n<N_seeds){
			criterion = TRUE;
			if(verbose>1){cat("   - Only",n,"seeds found in",niter,"iterations, end of refinement step\n")}
		}
		else{
		  if(verbose>1){cat("   -",n,"seed(s) founded after",niter,"iterations\n")}
		  if(verbose>0){cat(" * Generate N_alphaLOO =",N_alphaLOO,"samples with the weighted margin probability\n")}
		  
			candidate = do.call(cbind, lapply(1:N_alphaLOO, function(iter){
			    W = array(rnorm(candidate), dim = dim(candidate))
			    sigma = apply(candidate, 1, sd)*2
			    y = candidate + sigma*W
			    ratio = wMP(y)/wMP(candidate)
			    sel = ratio>runif(ratio)
			    candidate[,sel] <<- y[,sel]
			    candidate
			}))

			#plotting part
			if(plot==TRUE){
			  if(verbose>0){cat(" * 2D PLOT \n")}
			  print(p_meta +
			          stat_bin2d(data = as.data.frame(t(candidate)), bins = 20*Ru) +
			          scale_fill_gradientn(colours = rainbow(4))
			        )
			}

		  if(verbose>0){cat(" * Get K_alphaLOO =",K_alphaLOO,"points by clustering\n")}
			candidate = tryCatch(as.matrix(t(kmeans(t(candidate), centers=K_alphaLOO, iter.max=30)$centers)),
						error = function(cond){
							message(cond);
							r = rankMatrix(candidate);
							res = as.matrix(t(kmeans(t(candidate), centers=r, iter.max=30)$centers))
							return(res)
						})

			#plotting part
			if(plot==TRUE){
			  if(verbose>0){cat(" * 2D PLOT \n")}
			  print(p_meta + geom_point(data = as.data.frame(t(candidate)), color = "red", size = 4) )
			}
			
		  if(verbose>0){cat(" * Calculate performance function G on candidates\n")}
			eval = lsf(candidate);Ncall = Ncall + dim(candidate)[2]

		  if(verbose>0){cat(" * Add points to he learning database\n")}
			learn_db = cbind(learn_db,candidate)
			lsf_value = c(lsf_value,eval)

		  if(verbose>0){cat(" * Train the model\n")}
			if(learn_each_train==TRUE) {
			  if(verbose>1){cat("   - Learn hyperparameters\n")}
				meta = trainModel(design=learn_db,
						  response=lsf_value,
						  kernel=kernel,type="Kriging")
			}
			else {
			  if(verbose>1){cat("   - Use previous hyperparameters\n")}
				meta = trainModel(meta_model,
						  updesign=candidate,
						  upresponse=eval,type="Kriging")
			}

		  if(verbose>0){cat("\n * UPDATE quantities based on kriging surrogate model : MP, wMP, pi\n")}
			meta_model = meta$model
			meta_fun = meta$fun

			k = k + K_alphaLOO
	
			if(verbose>0){cat(" * Calculate alphaLOO \n")}
			LOO = leaveOneOut.km(meta_model, type="UK", trend.reestim=FALSE)
			piLOO = pnorm((failure-LOO$mean)/LOO$sd)
			notNull = piLOO>(10^(-16))
			if(verbose>1){cat("   -",sum(piLOO<10^(-16)),"samples not considered as pi<10^-16\n")}
			alphaLOO = mean(1*(lsf_value[notNull]<failure)/piLOO[notNull])
			if(verbose>1){cat("   - alphaLOO =",alphaLOO,"\n")}
			criterion = (alpha_int[1]<alphaLOO)*(alphaLOO<alpha_int[2])*(k>=Nmin) + (k>Nmax)
			
			#plotting part
			if(plot==TRUE | (limited_plot && criterion) ){
			  if(verbose>0){cat(" * 2D PLOT \n")}
        p <- p + geom_point(data = data.frame(t(learn_db), z = lsf_value), aes(color=z))
        z_meta = meta_fun(t(df_plot[,1:2]))
        df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(failure - z_meta$mean)/z_meta$sd)
        print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = failure) +
                geom_contour(data = df_plot_meta, aes(z=crit, color=..level.., alpha = 0.5), linetype = 4, breaks = k_margin))
			}
		}
	}

	cat("\n B- ESTIMATE AUGMENTED FAILURE PROBABILITY USING MC ESTIMATOR  \n")
	cat("    =========================================================== \n\n")

	if(!is.null(limit_fun_MH)){
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
		  if(verbose>0){cat(" * Generate Monte-Carlo population with an accept-reject strategy on a standard gaussian sampling\n")}
		  rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}
		  U = generateWithAR(dimension=dimension,
		                     N=N,
		                     limit_f=limit_fun_MH,
		                     failure=failure_MH,
		                     rand=rand)

		  if(verbose>0){cat(" * Calculate Monte-Carlo estimate\n")}
			Ind = pi(U)
			P_epsilon <- MC_est <- mean(Ind)
			VA_var = var(Ind)
			MC_var = VA_var/N
			cov_epsilon = sqrt(MC_var)/MC_est
		}
	}
	else{
	  if(verbose>0){cat(" * Generate standard gaussian samples\n")}
		U = matrix(rnorm(dimension*N),dimension,N)
	  if(verbose>0){cat(" * Calculate Monte-Carlo estimate\n")}
		Ind = pi(U)
		P_epsilon <- MC_est <- mean(Ind)
		VA_var = var(Ind)
		MC_var = VA_var/N
		cov_epsilon = sqrt(MC_var)/MC_est
	}

	points=U[,Ind>0.5]
	cat(" P_epsilon =",P_epsilon,"\n")
	cat(" cov_epsilon =",cov_epsilon,"\n")

	if(cov_epsilon>precision) {
		N = ceiling(VA_var/(precision^2*MC_est^2))
		cat(" * cov_epsilon too large ; this order of magnitude for the probabilty brings N =",N,"\n")
	}
# }


#Adaptative importance sampling scheme
cat("\n===========================================================================\n")
cat("\n STEP 2 : Adaptative importance sampling scheme \n")
cat("\n===========================================================================\n\n")

if(verbose>0){cat(" * Define h PDF \n")}
	h = function(x) {
		x = as.matrix(x)
		res = pi(x)*dmvnorm(t(x))/P_epsilon
		return(res)
	}

if(verbose>0){cat(" * Generate samples according to h with Metropolis-Hastings and calculate MC estimator\n")}
	N_alpha_max = Ncall_max - Ncall
  if(verbose>1){cat("   - Calculate approximated optimal density on the learning database\n")}
	h_learn_db = h(learn_db)
	if(verbose>1){cat("   - Find samples whose value is > 0\n")}
	if(verbose>1){cat("     ? seeds =",sum(h_learn_db>0),"samples in learn_db | h > 0\n")}
	if(verbose>0){cat("   - Calculate quasi-optimal density fonction on previous Monte-Carlo population\n")}
	h_U = h(U)
	if(verbose>1){cat("     ? seeds =",sum(h_U>0),"samples in MC pop | h > 0\n")}
	h_U = head(order(h_U, decreasing = TRUE), N_alpha_max - sum(h_learn_db>0))
	U = cbind(as.matrix(learn_db[,h_learn_db>0]),as.matrix(U[,h_U]))

if(verbose>1){cat("   - Generate N_alpha_max =",N_alpha_max,"points from",dim(U)[2],"seeds with MH\n")}
	# U = do.call(cbind, lapply(1:(ceiling(N_alphaLOO/dim(U)[2])*burnin), function(iter){
	 replicate(burnin, { 
	  W = array(rnorm(U), dim = dim(U))
	  sigma = apply(U, 1, sd)
	  y = U + sigma*W
	  ratio = h(y)/h(U)
	  sel = ratio>runif(ratio)
	  U[,sel] <<- y[,sel]
	 })
	# })[c(1:ceiling(N_alphaLOO/dim(U)[2]))*burnin])
	
	inDB = which(duplicated(cbind(U, learn_db))) - dim(U)[2]
	U = unique(cbind(learn_db[,inDB], U), M = 2)
	G = lsf_value[inDB]

	if(plot==TRUE) {
	  if(verbose>1){cat("   - 2D PLOT \n")}
	  row.names(U) = rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]
	  print(p_meta + 
	          stat_bin2d(data = as.data.frame(t(U)), bins = 20*Ru) +
	          scale_fill_gradientn(colours = rainbow(4))
	  )
	}

	cov_alpha = Inf
	while((cov_alpha>precision)*(N_alpha<N_alpha_max)) {
	  if(verbose>1){cat("   - Monte-Carlo estimator of alpha\n")}
	  if(verbose>1){cat("     ? Select randomly N_alpha = ",N_alpha," in the working population\n")}
		indices = tryCatch(c(indices,sample(c(1:N_alpha_max)[-indices],(N_alpha-length(indices)),replace=F)),
				error=function(cond) {return(sample(c(1:N_alpha_max),N_alpha,replace=F))})

	  if(verbose>1){cat("     ? Evaluate lsf on these samples if necessary\n")}
		isNAinG = is.na(G[indices]);
	  if(verbose>1){cat("      +",sum(isNAinG),"samples to evaluate in N_alpha =",N_alpha," samples\n")}
		G[indices][isNAinG] = lsf(U[,indices[isNAinG]]); Ncall = Ncall + sum(isNAinG)
		learn_db = cbind(learn_db,U[,indices[isNAinG]])
		lsf_value = c(lsf_value,G[indices][isNAinG])

	  if(verbose>1){cat("     ? Evaluate meta-model derived indicatrice pi on these samples\n")}
		Ind_meta = pi(U[,indices])
		Ind_lsf = (G[indices]<failure)

	  if(verbose>0){cat("#Evaluate kriging indicatrice likeness\n")}
	  if(verbose>0){cat(" mean(Ind_meta) =",mean(Ind_meta),"\n")}

	  if(verbose>0){cat("#Evaluate alpha estimator\n")}
		Ind = Ind_lsf/Ind_meta
		alpha = mean(Ind)
		VA_var = var(Ind)
		MC_var = VA_var/N_alpha
		VA_cov = sd(Ind)/alpha
		cov_alpha <- MC_cov <- sqrt(1/N_alpha)*VA_cov
		cat(" alpha =",alpha,"\n")
		cat(" cov_alpha =",cov_alpha,"\n")

		if(cov_alpha>precision) {
			N_alpha = ceiling(VA_var/(precision*alpha)^2)
			cat("#cov_alpha too large ; this order of magnitude for alpha brings N_alpha =",N_alpha,"\n")
			if(N_alpha>N_alpha_max) {
				cat("#N_alpha =",N_alpha,"> N_alpha_max =",N_alpha_max," => N_alpha = N_alpha_max\n")
				N_alpha = N_alpha_max;
			}
		}
	}
	if(N_alpha==N_alpha_max){
		if(verbose>0){cat("#Evaluate lsf on these samples if necessary\n")}
		isNAinG = is.na(G);
		if(verbose>1){cat(sum(isNAinG),"samples to evaluate in N_alpha =",N_alpha," samples\n")}
		G[isNAinG] = lsf(U[,isNAinG]); Ncall = Ncall + sum(isNAinG)
		learn_db = cbind(learn_db,U[,isNAinG])
		lsf_value = c(G,G[isNAinG])

		if(verbose>0){cat("#Evaluate meta-model derived indicatrice pi on these samples\n")}
		Ind_meta = pi(U)
		Ind_lsf = (G<failure)

		if(verbose>0){cat("#Evaluate alpha estimator\n")}
		Ind = Ind_lsf/Ind_meta
		alpha = mean(Ind)
		VA_var = var(Ind)
		MC_var = VA_var/N_alpha
		VA_cov = sd(Ind)/alpha
		cov_alpha <- MC_cov <- sqrt(1/N_alpha)*VA_cov
		cat(" alpha =",alpha,"\n")
		cat(" cov_alpha =",cov_alpha,"\n")
	}

#Results
cat("==========================================================================================",
"                              End of Meta-IS algorithm",
"==========================================================================================",sep="\n")
P = P_epsilon*alpha
cov = sqrt(cov_epsilon^2 + cov_alpha^2 + cov_epsilon^2*cov_alpha^2)

cat("   - P_epsilon =",P_epsilon,"\n")
cat("   - 95% conf. interv. on P_epsilon:", P_epsilon*(1-2*cov_epsilon),"< p <", P_epsilon*(1+2*cov), "\n")
cat("   - alpha =",alpha,"\n")
cat("   - 95% conf. interv. on alpha:", alpha*(1-2*cov_alpha),"< alpha <", alpha*(1+2*cov_alpha), "\n")
cat("   - p =",P,"\n")
cat("   - 95% conf. interv. on p:", P*(1-2*cov),"< p <",P*(1+2*cov),"\n")

if(plot + limited_plot){
    res = list(proba=P,
               cov=cov,
               Ncall=Ncall,
               learn_db=learn_db,
               lsf_value=(-1)^(!lower.tail)*lsf_value,
               meta_fun=meta_fun,
               meta_model=meta_model,
               points=points);
} else {res = list(proba=P,
		cov=cov,
		Ncall=Ncall,
		learn_db=learn_db,
		lsf_value=(-1)^(!lower.tail)*lsf_value,
		meta_fun=meta_fun,
		meta_model=meta_model,
		points=points)}



return(res)

}
 
