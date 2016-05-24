#' @title Subset by Support vector Margin Algorithm for Reliability esTimation
#' 
#' @description \code{S2MART} introduces a metamodeling step at each subset simulation
#' threshold, making number of necessary samples lower and the probability estimation
#' better according to subset simulation by itself.
#' 
#' @author Clement WALTER \email{clement.walter@cea.fr}
#' 
#' @details
#' 
#'   S2MART algorithm is based on the idea that subset simulations conditional
#'   probabilities are estimated with a relatively poor precision as it
#'   requires calls to the expensive-to-evaluate limit state function and
#'   does not take benefit from its numerous calls to the limit state function
#'   in the Metropolis-Hastings algorithm. In this scope, the key concept is
#'   to reduce the subset simulation population to its minimum and use it only
#'   to estimate crudely the next quantile. Then the use of a metamodel-based
#'   algorithm lets refine the border and calculate an accurate estimation of
#'   the conditional probability by the mean of a crude Monte-Carlo.
#' 
#' In this scope, a compromise has to be found between the two sources of
#' calls to the limit state function as total number of calls = (\code{Nn} +
#' number of calls to refine the metamodel) x (number of subsets) :
#'  \itemize{
#'    \item{\code{Nn} calls to find the next threshold value : the bigger \code{Nn},
#'    the more accurate the \sQuote{decreasing speed} specified by the
#'    \code{alpha_quantile} value and so the smaller the number of subsets}
#'    \item{total number of calls to refine the metamodel at each threshold}
#'  }
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
#'   J.-M. Bourinet, F. Deheeger, M. Lemaire:\cr
#'   \emph{Assessing small failure probabilities by combined Subset Simulation and Support Vector Machines}\cr
#'   Structural Safety (2011)
#'   
#'   \item
#'   F. Deheeger:\cr
#'   \emph{Couplage m?cano-fiabiliste : 2SMART - m?thodologie d'apprentissage stochastique en fiabilit?}\cr
#'       PhD. Thesis, Universit? Blaise Pascal - Clermont II, 2008
#' 
#'     \item
#'       S.-K. Au, J. L. Beck:\cr
#'       \emph{Estimation of small failure probabilities in high dimensions by Subset Simulation} \cr
#'       Probabilistic Engineering Mechanics (2001)
#' 
#'     \item
#'       A. Der Kiureghian, T. Dakessian:\cr
#'       \emph{Multiple design points in first and second-order reliability}\cr
#'      Structural Safety, vol.20 (1998)
#'
#'     \item
#'       P.-H. Waarts:\cr
#'       \emph{Structural reliability using finite element methods: an appraisal of DARS:\cr Directional Adaptive Response Surface Sampling}\cr
#'      PhD. Thesis, Technical University of Delft, The Netherlands, 2000
#'   }
#' 
#' @seealso
#' \code{\link{SMART}}
#' \code{\link{SubsetSimulation}}
#' \code{\link{MonteCarlo}}
#' \code{\link[DiceKriging]{km}} (in package \pkg{DiceKriging})
#' \code{\link[e1071]{svm}} (in package \pkg{e1071})
#' 
#' @examples 
#' \dontrun{
#'   res = S2MART(dimension = 2,
#'                lsf = kiureghian,
#'                N1 = 1000, N2 = 5000, N3 = 10000,
#'                plot = TRUE)
#'   
#'   #Compare with crude Monte-Carlo reference value
#'   reference = MonteCarlo(2, kiureghian, N_max = 500000)
#' }
#' 
#' #See impact of metamodel-based subset simulation with Waarts function :
#' \dontrun{
#'   res = list()
#'   # SMART stands for the pure metamodel based algorithm targeting directly the
#'   # failure domain. This is not recommended by its authors which for this purpose
#'   # designed S2MART : Subset-SMART
#'   res$SMART = mistral:::SMART(dimension  = 2, lsf = waarts, plot=TRUE)
#'   res$S2MART = S2MART(dimension = 2,
#'                       lsf = waarts,
#'                       N1 = 1000, N2 = 5000, N3 = 10000,
#'                       plot=TRUE)
#'   res$SS = SubsetSimulation(dimension = 2, waarts, n_init_samples = 10000)
#'  res$MC = MonteCarlo(2, waarts, N_max = 500000)
#' }
#' 
#' @import ggplot2
#' @import e1071
#' @import Matrix
#' @import mvtnorm
#' @export

S2MART = function(dimension,
                  #' @param dimension the dimension of the input space
                  lsf,
                  #' @param lsf the function defining the failure domain. Failure is lsf(X) < \code{failure}
            			## Algorithm parameters
            			Nn = 100,
            			#' @param Nn number of samples to evaluate the quantiles in the subset step
            			alpha_quantile = 0.1,
            			#' @param alpha_quantile cutoff probability for the subsets
            			failure = 0,
            			#' @param failure the failure threshold
            			## Meta-model choice
            
            			...,
            			#' @param ... All others parameters of the metamodel based algorithm
            			## Plot information
            			plot=FALSE,
            			#' @param plot to produce a plot of the failure and safety domain. Note that this requires a lot of
            			#' calls to the \code{lsf} and is thus only for training purpose
            			output_dir=NULL,
            			#' @param output_dir to save the plot into the given directory. This will be pasted with "_S2MART.pdf"
                  verbose = 0) {
                  #' @param verbose either 0 for almost no output, 1 for medium size output and 2 for all outputs

cat("==========================================================================================\n")
cat("                 Beginning of Metamodelisation on Subset Simulation algorithm \n")
cat("==========================================================================================\n")

# Fix NOTE issue with R CMD check
x <- z <- ..level.. <- NULL

	i = 0;			# Subset number
	y0 = Inf;		# first level
	y = NA;			# will be a vector containing the levels
	meta_fun = list(NA);	# initialise meta_fun list, that will contain lsf approximation at each iteration
	z_meta = list(NA)	# initialise z_meta list, that will contain evaluation of the surrogate model on the grid
	P = 1;			# probability
	Ncall = 0;		# number of calls to the lsf
	delta2 = 0;		# theoretical cov
	z_lsf = NULL		# For plotting part

	#Define the list variables used in the core loop
	U = list(Nn=matrix(nrow=dimension,ncol=Nn))
	#G stands for the value of the limit state function on these points
	G = list(g=NA,#value on learn_db points, ie all the value already calculated at a given iteration
		Nn=NA*c(1:Nn))
	
	#beginning of the core loop
	while (y0>0) {

		i = i+1;
		cat("\n SUBSET NUMBER ",i,"\n")
		cat(" -------------\n\n")
		
		# plotting part
		if(plot==TRUE){
		  if(!is.null(output_dir)) {
		    fileDir = paste(output_dir,"_S2MART.pdf",sep="")
		    pdf(fileDir)
		  }
		  xplot <- yplot <- c(-80:80)/10
		  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
		  p <- ggplot(data = df_plot, aes(x,y)) +
		    geom_contour(aes(z=z, color=..level..), breaks = failure) +
		    theme(legend.position = "none") +
		    xlim(-8, 8) + ylim(-8, 8)
		  
		  print(p)
		}
		

		cat("===========================================================================\n")
		cat(" STEP 1 : Subset Simulation part \n")
		cat("===========================================================================\n\n")
	
		if(i==1){
		  if(verbose>0){cat(" * Generate Nn =",Nn," standard gaussian samples\n")}
			U$Nn = matrix(rnorm(dimension*Nn, mean=0, sd=1),dimension,Nn)
		}
		else {
		  if(verbose>0){cat(" * Generate Nn = ",Nn," points from the alpha*Nn points lying in F(",i-1,") with MH algorithm\n",sep="")}
			U$Nn = generateWithlrmM(seeds=U$Nn[,G$Nn<y0],seeds_eval=G$Nn[G$Nn<y0],N=Nn,limit_f=meta_fun[[i-1]])$points
		}
	
		#assessment of g on these points
		if(verbose>0){cat(" * Assessment of the LSF on these points\n")}
		G$Nn = lsf(U$Nn);Ncall = Ncall + Nn
	
		#Determination of y[i] as alpha-quantile of Nn points Un=U$Nn
		if(verbose>0){cat(" * Determination of y[",i,"] as alpha-quantile of these samples\n",sep="")}
		y0 <- y[i] <- getQuantile(data=G$Nn,alpha=alpha_quantile)

	
		#Add points U$Nn to the learning database
		if(verbose>0){cat(" * Add points U$Nn to the learning database\n\n")}
		if(i==1) {
			learn_db = cbind(seq(0,0,l=dimension),U$Nn)
			g0 = lsf(as.matrix(seq(0,0,l=dimension)));Ncall = Ncall + 1;
			G$g = c(g0,G$Nn)
		}
		else {
			learn_db = cbind(learn_db,U$Nn)
			G$g = c(G$g,G$Nn)
		}
		rownames(learn_db) <- rep(c('x', 'y'), length.out = dimension)
		
		if(plot==TRUE){
		  p <- p + geom_point(data = data.frame(t(learn_db), z = G$g), aes(color = z))
		  print(p)
		}

		cat("===========================================================================\n")
		cat(" STEP 2 : Metamodel algorithm part \n")
		cat("===========================================================================\n\n")
		if(i==1){
			arg = list(dimension=dimension,
					lsf=lsf,
					failure = y0,
					learn_db = learn_db,
					lsf_value = G$g,
					plot = plot,
					z_lsf = z_lsf,
					add = TRUE,
					output_dir = output_dir,
          verbose = verbose,...)
			    meta_step = do.call(SMART,arg)
		}
		else{
			arg = list(dimension=dimension,
					lsf = lsf,
					failure = y0,
					learn_db = learn_db,
					lsf_value = G$g,
					seeds = seeds,
					seeds_eval = seeds_meta,
					limit_fun_MH = meta_fun[[i-1]],
					z_lsf = z_lsf,
					plot = plot,
					add = TRUE,
					output_dir = output_dir,
          verbose = verbose,...)
				  meta_step = do.call(SMART,arg)
		}

		if(!is.null(output_dir)){
		  if(verbose>0){cat("\n * 2D PLOT : CLOSE DEVICE \n")}
			dev.off()
		}

		P = P*meta_step$proba
		delta2 = tryCatch(delta2 + (meta_step$cov)^2,error = function(cond) {return(NA)})
		Ncall = Ncall + meta_step$Ncall
		learn_db = meta_step$learn_db
		G$g = meta_step$lsf_value
		meta_fun[[i]] = meta_step$meta_fun
		meta_model = meta_step$meta_model
		seeds = meta_step$points
		seeds_meta = meta_step$meta_eval
		z_meta[[i]] = meta_step$z_meta

		if(y0>0) {cat("\n * Current threshold =",y0,"> 0 => start a new subset\n")
			  cat("   - Current probability =",P,"\n")
			  cat("   - Current number of call =",Ncall,"\n")}
		else {cat("\n * Current threshold =",y0,"=> end of the algorithm\n")
		      cat("   - Final probability =",P,"\n")
		      cat("   - Total number of call =",Ncall,"\n")}
	}

	#' @return   An object of class \code{list} containing the failure probability
	#' and some more outputs as described below:
	res = list(p=P,
	           #' \item{p}{The estimated failure probability.}
             cov=sqrt(delta2),
	           #' \item{cov}{The coefficient of variation of the Monte-Carlo probability
	           #' estimate.}
             Ncall=Ncall,
	           #' \item{Ncall}{The total number of calls to the \code{lsf}.}
             learn_db=learn_db,
	           
	           #' \item{learn_db}{The final learning database, ie. all points where \code{lsf}
	           #' has been calculated.}
             lsf_value=G$g,
	           
	           #' \item{lsf_value}{The value of the \code{lsf} on the learning database.}
             meta_model=meta_model
	           #' \item{meta_model}{The final metamodel. An object from \pkg{e1071}.}
             );
	
	return(res)
}
