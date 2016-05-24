#' @title Crude Monte Carlo method
#' 
#' @description Estimate a failure probability using a crude Monte Carlo method.
#' 
#' @author Clement WALTER \email{clement.walter@cea.fr}
#' 
#' @details This implementation of the crude Monte Carlo method works with evaluating
#' batchs of points sequentialy until a given precision is reached on the final
#' estimator
#' 
#' @note Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}}
#' to do so. Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow}
#' = \code{dimension} and \sQuote{ncol} = number of vector to be consistent with
#' \code{as.matrix} transformation of a vector.
#' 
#' Algorithm calls lsf(X) (where X is a matrix as defined previously) and expects a vector
#' in return. This allows the user to optimise the computation of a batch of points,
#' either by vectorial computation, or by the use of external codes (optimised C or
#' C++ codes for example) and/or parallel computation.
#' 
#' @references
#'   \itemize{
#'    \item
#'      R. Rubinstein and D. Kroese:\cr
#'      \emph{Simulation and the Monte Carlo method} \cr
#'      Wiley (2008)\cr
#'  }
#' 
#' @seealso
#' \code{\link{SubsetSimulation}}
#' \code{\link{foreach}}
#' 
#' @examples
#' #First some considerations on the usage of the lsf. 
#' #Limit state function defined by Kiureghian & Dakessian :
#' # Remember you have to consider the fact that the input will be a matrix ncol >= 1
#' lsf_wrong = function(x, b=5, kappa=0.5, e=0.1) {
#'   b - x[2] - kappa*(x[1]-e)^2 # work only with a vector of lenght 2
#' }
#' lsf_correct = function(x){
#'   apply(x, 2, lsf_wrong)
#' }
#' lsf = function(x, b=5, kappa=0.5, e=0.1) {
#'   x = as.matrix(x)
#'   b - x[2,] - kappa*(x[1,]-e)^2 # vectorial computation, run fast
#' }
#' 
#' y = lsf(X <- matrix(rnorm(20), 2, 10))
#' #Compare running time
#' \dontrun{
#'   require(microbenchmark)
#'   X = matrix(rnorm(2e5), 2)
#'   microbenchmark(lsf(X), lsf_correct(X))
#' }
#' 
#' #Example of parallel computation
#' require(doParallel)
#' lsf_par = function(x){
#'  foreach(x=iter(X, by='col'), .combine = 'c') %dopar% lsf(x)
#' }
#' 
#' #Try Naive Monte Carlo on a given function with different failure level
#' \dontrun{
#'   res = list()
#'   res[[1]] = MonteCarlo(2,lsf,q = 0,plot=TRUE)
#'   res[[2]] = MonteCarlo(2,lsf,q = 1,plot=TRUE)
#'   res[[3]] = MonteCarlo(2,lsf,q = -1,plot=TRUE)
#'   
#' }
#' 
#' 
#' #Try Naive Monte Carlo on a given function and change number of points.
#' \dontrun{
#'   res = list()
#'   res[[1]] = MonteCarlo(2,lsf,N_max = 10000)
#'   res[[2]] = MonteCarlo(2,lsf,N_max = 100000)
#'   res[[3]] = MonteCarlo(2,lsf,N_max = 500000)
#' }
#' 
#' @import ggplot2
#' 
#' @export
MonteCarlo = function(dimension, 
                      #' @param dimension the dimension of the input space.
                      lsf,
                      #' @param lsf the function defining safety/failure domain.
                      N_max 	= 500000,
                      #' @param N_max maximum number of calls to the \code{lsf}.
		                  N_batch 	= 1000,
		                  #' @param N_batch number of points onte evalutae the \code{lsf} at each iteration.
                      q	= 0,
		                  #' @param q the quantile
		                  lower.tail = TRUE,
		                  #' @param lower.tail as for pxxxx functions, TRUE for estimating P(lsf(X) < q), FALSE
		                  #' for P(lsf(X) > q)
		                  precision	= 0.05,
		                  #' @param precision a targeted maximum value for the coefficient of variation.
                      plot 	= FALSE,
		                  #' @param plot to plot the contour of the \code{lsf} as well as the generated samples.
                      output_dir = NULL,
		                  #' @param output_dir to save a copy of the plot in a pdf. This name will be
		                  #' pasted with
		                  #' "_Monte_Carlo_brut.pdf"
		                  verbose 	= 0){
                      #' @param verbose to control the level of outputs in the console; either 0 or 1 or 2 for
                      #' almost no outputs to a high level output.

cat("========================================================================\n")
cat("                 Beginning of Monte-Carlo algorithm\n")
cat("========================================================================\n\n")
  
# Step 0 : Initialization
cov = Inf;
Ncall = 0;


# Fix NOTE issue with R CMD check
x <- y <- z <- ..level.. <- NULL

if(lower.tail==FALSE){
  lsf_dec = lsf
  lsf = function(x) -1*lsf_dec(x)
  q <- -q
}

if(verbose>0){cat(" * STEP 1 : FIRST SAMPLING AND ESTIMATION \n")}

if(verbose>1){cat("   - Generate N_batch = ",N_batch," standard samples\n")}
U = matrix(rnorm(dimension*N_batch),dimension,N_batch, dimnames = list(rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]))

if(verbose>1){cat("   - Evaluate LSF on these samples\n")}
G = lsf(U); Ncall = Ncall + N_batch

if(verbose>1){cat("   - Evaluate q probability and corresponding CoV \n")}
P = mean(G < q)
cov = sqrt((1-P)/(Ncall*P))

if(plot==TRUE){
  
  if(!is.null(output_dir)) {
    fileDir = paste(output_dir,"_Monte_Carlo_brut.pdf",sep="")
    pdf(fileDir)
  }
  xplot <- yplot <- c(-80:80)/10
  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
  p <- ggplot(data = df_plot, aes(x,y))
  indCol = (G>q)+1
  p <- p + geom_contour(aes(z=z, color = ..level..), breaks = 0) +
    geom_point(data = data.frame(t(U), z = G, ind = indCol), colour = factor(indCol)) +
    theme(legend.position = "none") +
    xlim(-8, 8) + ylim(-8, 8)
  print(p)
}

Nrem = min(N_max - Ncall, N_batch)

if( (verbose>0) && (cov > precision) && (Nrem > 0) ) {cat(" * STEP 2 : LOOP UNTIL COV < PRECISION \n")}

while( (cov > precision) && (Nrem > 0) ) {

	if(verbose>0){cat(" * cov =",cov,">",precision,"and",N_max - Ncall,"remaining calls to the LSF \n")}
	if(verbose>1){cat("   - Generate N =",Nrem,"standard samples\n")}
	U = matrix(rnorm(dimension*Nrem),dimension,Nrem, dimnames = list(rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]))
	
	if(verbose>1){cat("   - Evaluate LSF on these samples\n")}
	G = c(G,lsf(U)); Ncall = Ncall + Nrem
	
	if(plot==TRUE){
	  if(verbose>1){cat("   - Plot these samples\n")}
	  indCol = c(tail(G, Nrem)>q)+1
	  p <- p + geom_point(data = data.frame(t(U), z = tail(G, Nrem), ind = indCol), colour= factor(indCol))
	  print(p)
	}
	
	if(verbose>1){cat("   - Evaluate q probability and corresponding CoV \n")}
	P = mean(G < q)
	cov = sqrt((1-P)/(Ncall*P))

	if(verbose>1){cat("   - P =",P,"\n")}
	if(verbose>1){cat("   - cov =",cov,"\n")}

	Nrem = min(N_max - Ncall, N_batch)
}

if(plot == TRUE) {
	if(!is.null(output_dir)) {dev.off()}
}

cat("========================================================================\n")
cat("                    End of Monte-Carlo algorithm\n")
cat("========================================================================\n\n")

cat("   - p =",P,"\n")
cat("   - q =", q, "\n")
cat("   - 95% confidence intervalle :",P*(1-2*cov),"< p <",P*(1+2*cov),"\n")
cat("   - cov =",cov,"\n")
cat("   - Ncall =",Ncall,"\n")

ecdf_MC <- local({
  lower.tail <- lower.tail
  G <- G*(-1)^(!lower.tail)
  function(q) {
    if(lower.tail==TRUE){
      p <- mean(G<q)
    }
    else{
      p <- mean(G>q)
    }
    p
  }
})

#' @return An object of class \code{list} containing the failure probability and some
#' more outputs as described below:
res = list(p=P,
           #' \item{p}{the estimated probabilty.}
           ecdf_MC = ecdf_MC,
           #' \item{ecdf_MC}{the empiracal cdf got with the generated samples.}
           cov=cov,
           #' \item{cov}{the coefficient of variation of the Monte Carlo estimator.}
           Ncall=Ncall,
           #' \item{Ncall}{the total numnber of calls to the \code{lsf}, ie the total
           #' number of generated samples.}
           X = U,
           #' \item{X}{the generated samples.}
           Y = G*(-1)^(!lower.tail)
           #' \item{Y}{the value \code{lsf(X)}.}
           )
return(res)
}