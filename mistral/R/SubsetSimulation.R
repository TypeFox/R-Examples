#' @title Subset Simulation Monte Carlo
#' 
#' @description Estimate a probability of failure with the Subset Simulation algorithm (also known as
#' Multilevel Splitting or Sequential Monte Carlo for rare events).
#' 
#' @aliases ss subset
#' 
#' @author Clement WALTER \email{clement.walter@cea.fr}
#' 
#' @details This algorithm uses the property of conditional probabilities on nested subsets
#' to calculate a given probability defined by a limit state function.
#' 
#' It operates iteratively on \sQuote{populations} to estimate the quantile
#' corresponding to a probability of \code{p_0}. Then, it generates
#' samples conditionnaly to this threshold, until found threshold be lower
#' than 0.
#' 
#' Finally, the estimate is the product of the conditional probabilities.
#' 
#' @note Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}}
#' to do so. Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow}
#' = \code{dimension} and \sQuote{ncol} = number of vector to be consistent with
#' \code{as.matrix} transformation of a vector.
#' 
#' Algorithm calls lsf(X) (where X is a matrix as defined previously) and expects a vector
#' in return. This allows the user to optimise the computation of a batch of points,
#' either by vectorial computation, or by the use of external codes (optimised C or
#' C++ codes for example) and/or parallel computation; see examples in \link{MonteCarlo}.
#' 
#' @seealso 
#' \code{\link{IRW}}
#' \code{\link{MP}}
#' \code{\link{MonteCarlo}}
#' 
#' @references 
#'   \itemize{
#'    \item
#'      S.-K. Au, J. L. Beck:\cr
#'      \emph{Estimation of small failure probabilities in high dimensions by Subset Simulation} \cr
#'      Probabilistic Engineering Mechanics (2001)\cr
#'    \item A. Guyader, N. Hengartner and E. Matzner-Lober:\cr
#'     \emph{Simulation and estimation of extreme quantiles and extreme
#'     probabilities}\cr
#'     Applied Mathematics \& Optimization, 64(2), 171-196.\cr
#'    \item F. Cerou, P. Del Moral, T. Furon and A. Guyader:\cr
#'    \emph{Sequential Monte Carlo for rare event estimation}\cr
#'    Statistics and Computing, 22(3), 795-808.\cr
#'  }
#'  
#' @examples 
#' #Try Subset Simulation Monte Carlo on a given function and change number of points.
#'  
#' \dontrun{
#'  res = list()
#'  res[[1]] = SubsetSimulation(2,kiureghian,N=10000)
#'  res[[2]] = SubsetSimulation(2,kiureghian,N=100000)
#'  res[[3]] = SubsetSimulation(2,kiureghian,N=500000)
#' }
#' 
#' # Compare SubsetSimulation with MP
#' \dontrun{
#' p <- res[[3]]$p # get a reference value for p
#' p_0 <- 0.1 # the default value recommended by Au \& Beck
#' N_mp <- 100
#' # to get approxumately the same number of calls to the lsf
#' N_ss <- ceiling(N_mp*log(p)/log(p_0))
#' comp <- replicate(50, {
#' ss <- SubsetSimulation(2, kiureghian, N = N_ss)
#' mp <- MP(2, kiureghian, N = N_mp, q = 0)
#' comp <- c(ss$p, mp$p, ss$Ncall, mp$Ncall)
#' names(comp) = rep(c("SS", "MP"), 2)
#' comp
#' })
#' boxplot(t(comp[1:2,])) # check accuracy
#' sd.comp <- apply(comp,1,sd)
#' print(sd.comp[1]/sd.comp[2]) # variance increase in SubsetSimulation compared to MP
#' 
#' colMeans(t(comp[3:4,])) # check similar number of calls
#' }
#' 
#' @import ggplot2
#' @import foreach
#' @export
SubsetSimulation = function(dimension,
                            #' @param dimension the dimension of the input space.
                            lsf,
                            #' @param lsf the function defining failure/safety domain.
                            p_0 = 0.1,
                            #' @param p_0 a cutoff probability for defining the subsets.
                            N = 10000,
                            #' @param N the number of samples per subset, ie the population size for the Monte
                            #' Carlo estimation of each conditional probability.
                            q = 0,
                            #' @param q the quantile defining the failure domain.
                            lower.tail = TRUE,
                            #' @param lower.tail as for pxxxx functions, TRUE for estimating P(lsf(X) < q), FALSE
                            #' for P(lsf(X) > q)
                            K,
                            #' @param K a transition Kernel for Markov chain drawing in the regeneration step.
                            #' K(X) should propose a matrix of candidate sample (same dimension as X) on which
                            #' \code{lsf} will be then evaluated and transition accepted of rejected. Default
                            #' kernel is the one defined K(X) = (X + sigma*W)/sqrt(1 + sigma^2) with W ~ N(0, 1).
                            burnin = 20,
                            #' @param burnin a burnin parameter for the the regeneration step.
                            save.all = FALSE,
                            #' @param save.all if TRUE, all the samples generated during the algorithms are saved
                            #' and return at the end. Otherwise only the working population is kept at each
                            #' iteration.
                            plot = FALSE,
                            #' @param plot to plot the contour of the \code{lsf} and the generated sample.
                            output_dir = NULL,
                            #' @param output_dir to save the plot into a pdf file. This variable will
                            #' be paster with
                            #' "_Subset_Simulation.pdf"
                            plot.lab = c('x', 'y'),
                            #' @param plot.lab the x and y labels for the plot
                            verbose = 0) {
  #' @param verbose Either 0 for almost no output, 1 for medium size output and 2 for all outputs
  
  cat("==========================================================================================\n")
  cat("                              Beginning of Subset Simulation algorithm \n")
  cat("==========================================================================================\n\n")
  
  # Fix NOTE issue with R CMD check
  x <- y <- z <- ..level.. <- NULL
  
  if(lower.tail==FALSE){
    lsf_dec = lsf
    lsf = function(x) -1*lsf_dec(x)
    q <- -q
  }
  
  if(verbose>0){cat("  * Generate the first N =",N,"samples by crude MC \n")}
  Utot <- U <- matrix(rnorm(dimension*N), nrow = dimension, dimnames = list(rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]))
  
  if(verbose>0){cat("  * Evaluate the limit state function \n")}
  Gtot <- G <- lsf(U);
  Ncall = N
  
  if(verbose>0){cat("  * Find the quantile q_0 verifying P[lsf(U) < q_0] = p_0 \n")}
  q_0 = max(quantile(G, probs = p_0), q)
  cat("   - q_0 =",q_0,"\n")
  
  if(verbose>0){cat("  * Evaluate actual probability \n")}
  indG = G<q_0
  P = mean(indG)
  cov = (1-P)/N/P
  
  n_subset = 1
  cat("   - P =",P,"\n")
  
  # Set sigma for transition Kernel K
  sigma.hist <- sigma <- 0.3
  
  if(plot==TRUE){
    
    if(!is.null(output_dir)) {
      fileDir = paste(output_dir,"_Subset_Simulation.pdf",sep="")
      pdf(fileDir)
    }
    xplot <- yplot <- c(-80:80)/10
    df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
    indCol = factor(indG+n_subset, levels = 1:20)
    p <- ggplot(data = df_plot, aes(x,y)) +
      geom_contour(aes(z=z, color=..level..), breaks = q_0) +
      geom_point(data = data.frame(t(U), z = G, ind = indCol), color=factor(indCol)) +
      theme(legend.position = "none") +
      xlim(-8, 8) + ylim(-8, 8) + xlab(plot.lab[1]) + ylab(plot.lab[2])
    print(p)
  }
  
  while(q_0>0){
    
    n_subset = n_subset + 1
    
    if(verbose>0){cat("  * Resample N(1-p_0) =",N*(1-p_0),"samples with MCMC \n")}
    U_tmp <- U[,indG]
    G_tmp <- G[indG]
    
    if(missing(K)) {
      K = function(x){
        W = array(rnorm(x), dim = dim(x))
        # sigma = apply(x, 1, sd)*2
        # sigma <- 0.3
        y = (x + sigma*W)/sqrt(1 + sigma^2)
        return(y)
      }
    }
    
    U <- G <- NULL
    acceptance.rate <- 0
    foreach::times(ceiling(N/sum(indG))) %do% {
      foreach::times(burnin) %do% {
        U_star <- K(U_tmp) # generate candidate
        G_star <- lsf(U_star) # calculate lsf
        if(save.all==TRUE){
          Utot <- cbind(U, U_star)
          Gtot <- c(Gtot, G_star)
        }
        Ncall <- Ncall + length(G_star) # update Ncall
        indG_star <- G_star<q_0 # get sample in the right domain
        acceptance.rate <- acceptance.rate + sum(indG_star)
        U_tmp[,indG_star] <- U_star[,indG_star] # update sample in the right domain
        G_tmp[indG_star] <- G_star[indG_star]
      }
      U <- cbind(U, U_tmp)
      G <- c(G, G_tmp)
    }
    acceptance.rate <- acceptance.rate/burnin/N
    if(acceptance.rate>0.4) sigma <- sigma*1.1
    if(acceptance.rate<0.2) sigma <- sigma*0.9
    sigma.hist <- c(sigma.hist, sigma)
    
    #     replicate(burnin, {
    #       U_star <- K(U)
    #       G_star <- lsf(U_star)
    #       Ncall <<- Ncall + N
    #       indG_star <- G_star<q_0
    #       U[,indG_star] <<- U_star[,indG_star]
    #       G[indG_star] <<- G_star[indG_star]
    #     })
    
    if(verbose>1){cat("  * Find the quantile q_0 verifying P[lsf(U) < q_0] = p_0 \n")}
    q_0 = max(quantile(G, probs = p_0), q)
    if(verbose>0){cat("  * Evaluate actual probability \n")}
    indG = G<q_0
    P_ss = mean(indG)
    P = P*P_ss
    cov = cov + (1-P_ss)/P_ss/N
    cat("   - q_0 =",q_0,"\n")
    cat("   - P =",P,"\n")
    
    if(plot==TRUE){
      indCol = factor(indG+n_subset, levels = 1:20)
      p <- p + geom_contour(aes(z=z, color=..level..), breaks = q_0) +
        geom_point(data = data.frame(t(U), z = G, ind = indCol), color= factor(indCol))
      print(p)
    }
    
  }
  
  if(!is.null(output_dir)) {dev.off()}
  
  cat("==========================================================================================\n")
  cat("                              End of Subset Simulation algorithm \n")
  cat("==========================================================================================\n\n")
  
  cat("   - p =",P,"\n")
  cat("   - q =", q, "\n")
  cat("   - 95% confidence intervalle :",P*(1-2*cov),"< p <",P*(1+2*cov),"\n")
  cat("   - cov =",cov,"\n")
  cat("   - Ncall =",Ncall,"\n")
  
  
  #return the result
  #' @return   An object of class \code{list} containing the failure probability and
  #' some more outputs as described below:
  result = list(p = P,
                #' \item{p}{the estimated failure probability.}
                cov = cov,
                #' \item{cov}{the estimated coefficient of variation of the estimate.}
                Ncall = Ncall,
                #' \item{Ncall}{the total number of calls to the \code{lsf}.}
                X = U,
                #' \item{X}{the working population.}
                Y = G*(-1)^(!lower.tail),
                #' \item{Y}{the value lsf(X).}
                Xtot = Utot,
                #' \item{Xtot}{if \code{save.list==TRUE}, all the \code{Ncall} samples generated by
                #' the algorithm.}
                Ytot = Gtot*(-1)^(!lower.tail),
                #' \item{Ytot}{the value lsf(Xtot).}
                sigma.hist = sigma.hist
                #' \item{sigma.hist}{if default kernel is used, sigma is initialized with 0.3 and
                #' then further adaptively updated to have an average acceptance rate of 0.3}
  )
  return(result)
  
}
