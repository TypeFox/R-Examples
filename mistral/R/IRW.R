#' @title Increasing Randow Walk
#' 
#' @description Simulate the increasing random walk associated with a real-valued continuous
#' random variable.
#' 
#' @aliases NestedSampling TPA
#' 
#' @author Clement WALTER \email{clement.walter@cea.fr}
#' 
#' @details This function lets generate the increasing random walk associated with a continous
#' real-valued random variable of the form \code{Y = lsf(X)} where \code{X} is
#' vectorial random variable.
#' 
#' This random walk can be associated with a Poisson process with parameter
#' \code{N} and hence the number of iterations before a given threshold \code{q}
#' is directly related to P[ lsf(X) > q]. It is the core tool of algorithms
#' such as nested sampling, Last Particle Algorithm or Tootsie Pop Algorithm.
#' 
#' Bascially for \code{N = 1}, it generates a sample \eqn{Y = lsf(X)} and iteratively
#' regenerates greater than the sought value: \eqn{Y_{n+1} \sim \mu^Y( \cdot \mid Y > Y_n}. This
#' regeneration step is done with a Metropolis-Hastings algorithm and that is why it is usefull
#' to consider generating several chains all together (\code{N > 1}).
#' 
#' The algorithm stops when it has simulated the required number of events \code{Nevent} or when
#' it has reached the sought threshold \code{q}.
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
#' \item C. Walter:\cr
#' \emph{Moving Particles: a parallel optimal Multilevel Splitting method
#' with application in quantiles estimation and meta-model based algorithms}\cr
#' Structural Safety, 55, 10-25.\cr
#' 
#' \item C. Walter:\cr
#' \emph{Point Process-based Monte Carlo estimation}\cr
#' arXiv preprint arXiv:1412.6368.\cr
#' 
#' \item J. Skilling:\cr
#' \emph{Nested sampling for general Bayesian computation}\cr
#' Bayesian Analysis, 1(4), 833-859.\cr
#' 
#' \item M. Huber \& S. Schott:\cr
#' \emph{Using TPA for Bayesian inference}\cr
#' Bayesian Statistics 9, 9, 257.\cr
#' 
#' \item A. Guyader, N. Hengartner and E. Matzner-Lober:\cr
#' \emph{Simulation and estimation of extreme quantiles and extreme
#' probabilities}\cr
#' Applied Mathematics \& Optimization, 64(2), 171-196.
#' }
#' 
#' @seealso
#' \code{\link{MP}}
#' 
#' @examples
#' # Get faililng samples for the kiureghian limit state function
#' # Failure is defined as lsf(X) < 0 so we have to invert the lsf
#' lsf <- function(x) -1*kiureghian(x)
#' \dontrun{
#' fail.samp <- IRW(2, lsf, q = 0, N = 10, plot = TRUE)
#' }
#' @import ggplot2
#' @export
IRW = function(dimension,
               #' @param dimension dimension of the input space.
                lsf,
               #' @param lsf limit state function.
                N = 10,             
               #' @param N number of particules.
                q = Inf,      
               #' @param q level until which the randow walk is to be generated.
                Nevent = Inf, 
               #' @param Nevent the number of desired events.
                particles,          
               #' @param particles to start with some given particles.
                LSF_particles = lsf(particles),      
               #' @param LSF_particles value of the \code{lsf} on these particles.
                K,                  
               #' @param K kernel transition for conditional generations.
                burnin = 20,             
               #' @param burnin burnin parameter.
                sigma = 0.3,        
               #' @param sigma radius parameter for \code{K}.
                last.return = TRUE,
               #' @param last.return if the last event should be returned.
                use.potential = TRUE,
               #' @param use.potential tu use a \sQuote{potential} matrix to select starting point not
               #' directly related to the sample to be moved with the MH algorithm.
                ## plot parameter
                plot = FALSE,          
               #' @param plot if \code{TRUE}, the algorithm plots the evolution of the particles. This
               #' requieres to evaluate the \code{lsf} on a grid and is only for visual purpose.
                print_plot = FALSE,       
               #' @param print_plot if TRUE, print the updated plot after each iteration. This might
               #' be slow; use with a small \code{N}. Otherwise it only prints the final plot.
                output_dir = NULL     
               #' @param output_dir if plots are to be saved in pdf in a given directory. This will
               #' be pasted with \sQuote{_IRW.pdf}. Together with \code{print_plot==TRUE} this will
               #' produce a pdf with a plot at each iteration, enabling \sQuote{video} reconstitution
               #' of the algorithm.
){

  cat("==========================================================================================\n")
  cat("                              Beginning of MP \n")
  cat("==========================================================================================\n\n")
  
  
  # Fix NOTE issue with R CMD check
  x <- z <- ..level.. <- NULL
  
  ## STEP 0 : INITIALISATION
  
  if(missing(K)) {
    K = function(x){
      W = rnorm(dimension)
      y = (x + sigma*W)/sqrt(1 + sigma^2)
      return(y)
    }
  }
  
  Ncall = 0;
  m = 1
  Ndep = 0;
  Ndup = 0
  acceptance = c()
  
  ## Setup potential matrix
  potentiel = matrix(1, nrow=N, ncol=N) - diag(1,N)
  
  ## Plotting part
  if(plot==TRUE){
    
    cat(" * 2D PLOT : SET-UP \n")
    if(!is.null(output_dir)){
      output_d = paste(output_dir,"_IRW.pdf",sep="")
      pdf(output_d)
    }
    xplot <- yplot <- c(-80:80)/10
    df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
    p <- ggplot(data = df_plot, aes(x,y, z = z)) +
      geom_contour(aes(color=..level..), breaks = q) +
      theme(legend.position = "none") +
      xlim(-8, 8) + ylim(-8, 8)
    if(print_plot) print(p)
    else{
      if(!is.null(output_dir)) list_plot <- list(p)
    }
  }
  
  cat(" =================================== \n")
  cat(" STEP 1 : FIRST SAMPLING AND MINIMUM \n")
  cat(" =================================== \n\n")
  
  if(missing(particles)){ 
    cat("## Draw an iid N =",N,"sample (X_i) \n")
    X = matrix(rnorm(dimension*N),ncol=N, dimnames = list(rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]))
    
    cat("## Calculate lsf \n")
    y = lsf(X); Ncall = Ncall + N
  }
  else{
    cat("## Restart chain from the particules \n")
    X = particles
    row.names(X) <- rep(c('x', 'y'), ceiling(dimension/2))[1:dimension]
    
    cat("## Get the LSF values \n")
    if(missing(LSF_particles)) {
      y = lsf(X)
    }
    else{
      y = LSF_particles
    }
    N = length(y)
  }
  
  if(plot==TRUE) {
    cat(" * 2D PLOT : UPDATE \n")
    p <- p + geom_point(data = data.frame(t(X), z = y), aes(color=z))
    if(print_plot) print(p)
    else{
      if(!is.null(output_dir)) list_plot <- c(list_plot, list(p))
    }
  }
  
  cat("## Find minimum \n\n")
  ind = which.min(y)
  L = y[ind]
  
  cat(" ================== \n")
  cat(" STEP 2 : CORE LOOP \n")
  cat(" ================== \n\n")
  
  while((L[Ndep+1] < q) && (Ndep <= Nevent)){
    cat("#### ITERATIION",m,"\n")
    cat("## MOVE nbr",Ndep+1,"\n\n")
    
    cat(" * Select randomly a particle to start sampling from \n")
    sel = tryCatch(sample(c(1:N)[potentiel[,ind]==1],1),error = function(cond) ind)      
    cat("   - mov. particle : ",ind,"; y =", y[ind],"\n")
    cat("   - sel. particle : ",sel,"; y =", y[sel],"\n")
    
    X_from  <-  X[,sel] 
    y_from  <-  y[sel]
    
    cat(" * Markov chain drawing \n")
    acceptance[Ndep + 1] = 0
    for(i in 1:burnin){
      ## Draw a new particule
      X_star = K(X_from)
      tryCatch(
        {y_star <- lsf(X_star); Ncall = Ncall + 1},
        error = function(cond) {
          message("Unable to evaluate the model at proposed point, transition refused \n");
          y_star <<- L[Ndep+1]-1
          return(y_star)
        }
      )
      
      ## Select it if in the right space
      if(y_star>L[Ndep+1]) {
        X_from = X_star; y_from = y_star;
        acceptance[Ndep + 1] = acceptance[Ndep + 1] +1
      }
    }
    
    if(y_from > L[Ndep+1]) {
      cat(" * New particle accepted \n")
      
      if(plot==TRUE) {
        cat(" * 2D PLOT : UPDATE \n")
        p <- p + geom_line(data = df_tmp <- data.frame(x = c(X[1,ind],X_from[1]), y = c(X[2,ind],X_from[2]), z = c(0,0)), color = "green", linetype = 4) +
          geom_point(data = data.frame(x = X_from[1], y = X_from[2], z = y_from), aes(color=z))
        if(print_plot) print(p)
        else{
          if(!is.null(output_dir)) list_plot <- c(list_plot, list(p))
        }
      }
      
      X[,ind] = X_from;
      y[ind] = y_from
      Ndep = Ndep + 1;
      
      cat(" * Refresh potential matrix \n")
      if(use.potential == TRUE){ 
        if(acceptance[Ndep]==0){
          potentiel[,ind] = potentiel[,sel]
          potentiel[ind,] = potentiel[sel,]
          Ndup = Ndup + 1
        }
        else{
          potentiel[,ind] = seq(1,1,l=N)
          potentiel[ind,] = seq(1,1,l=N)
          potentiel[ind,sel] = 0
          potentiel[ind,ind] = 0
        }
      }
      else{
        potentiel = matrix(1, nrow=N, ncol=N) - diag(1,N)
      }
      
      cat(" * Find new minimum \n")
      ind = which.min(y)
      L[Ndep+1] = y[ind]
      cat("   - current threshold :",L[Ndep+1],"\n\n")
    }
    else{
      cat(" * New particule rejected \n")
    }
    m = m + 1;
    
  }

  if(plot==TRUE) {
    cat(" * 2D PLOT : FINAL SAMPLING \n\n")
    
    if(!is.null(output_dir)){
      if(print_plot) print(p)
      else{
        list_plot <- c(list_plot, list(p))
        lapply(list_plot, print)
      }
      dev.off()
      output_d = paste(output_dir,"_MP_final_db.pdf",sep="")
      pdf(output_d)
    }
    print(p)
    print(ggplot(df_plot, aes(x,y)) +
            geom_contour(aes(z=z, color=..level..), breaks = q) +
            geom_point(data = data.frame(t(X), z = y) , aes(color = z)) +
            theme(legend.position = "none") +
            xlim(-8, 8) + ylim(-8, 8))
    if(!is.null(output_dir)){
      dev.off()
    }
  }
  
  cat("==========================================================================================\n")
  cat("                              End of MP \n")
  cat("==========================================================================================\n\n")

    cat("   - Number of iterations =",m-1,"\n")
    cat("   - Number of moves =",Ndep,"\n")
    cat("   - Number of wrong moves =",m-1-Ndep,"\n")
    cat("   - Total number of calls =",Ncall,"\n")
    
  if(last.return == FALSE){L = L[-length(L)]}
  
#' @return An object of class \code{list} containing the following data:
    res = list(L = L,
#' \item{L}{the events of the random walk.}
                              M = m-1,
#' \item{M}{the total number of iterations.}
               Ncall = Ncall,
#' \item{Ncall}{the total number of calls to the \code{lsf}.}
               particles = X,
#' \item{particles}{a matrix containing the final particles.}
               LSF_particles = y,
#' \item{LSF_particles}{the value of \code{lsf} on the \code{particles}.}
               q = q,
#' \item{q}{the threshold considered when generating the random walk.}
               Nevent = Nevent,
#' \item{Nevent}{the target number of events when generating the random walk.}
               Nwmoves = m-1-Ndep,
#' \item{Nwmoves}{the number of rejected transitions, ie when the proposed point was not stricly
#' greater/lower than the current state.}
               acceptance = acceptance/T)
#' \item{acceptance}{a vector containing the acceptance rate for each use of the MH algorithm.}
}