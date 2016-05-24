#' Markov Basis Metropolis-Hastings Algorithm
#'
#' Given a starting table (as a vector) and a loglinear model matrix A, compute the Markov basis of A with 4ti2 and then run the Metropolis-Hastings algorithm starting with the starting table.  
#'
#' See Algorithm 1.1.13 in LAS, the reference below.
#' 
#' @param init the initial step
#' @param moves the markov basis (the negatives will be added).  see ?markov
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param engine C++ or R? (C++ yields roughly a 20-25x speedup)
#' @return a list
#' @export metropolis
#' @author David Kahle
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
#' @examples
#'
#' \dontrun{
#'
#' 
#' 
#' data(handy)
#' 
#' exp   <- loglin(handy, as.list(1:2), fit = TRUE)$fit
#' e <- unname(tab2vec(exp))
#' h <- t(t(unname(tab2vec(handy))))
#' chisq <- algstat:::computeChisqsCpp(h, e)
#' 
#' out <- hierarchical(~ Gender + Handedness, data = handy)
#' chisqs <- algstat:::computeChisqsCpp(out$steps, e)
#' 
#' mean(chisqs >= chisq)
#' fisher.test(handy)$p.value
#' 
#' 
#' 
#' 
#' 
#' A <- hmat(c(2,2), as.list(1:2))
#' moves <- markov(A)
#' outC <- metropolis(tab2vec(handy), moves, 1e4, engine = "Cpp")
#' str(outC)
#' outR <- metropolis(tab2vec(handy), moves, 1e4, engine = "R", thin = 20)
#' str(outR)
#' 
#' # showSteps(out$steps)
#' 
#' 
#' library(microbenchmark)
#' microbenchmark(
#'   metropolis(tab2vec(handy), moves, engine = "Cpp"),
#'   metropolis(tab2vec(handy), moves, engine = "R")
#' )
#' 
#' # cpp ~ 20-25x faster
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' showSteps <- function(steps){
#'   apply(steps, 2, function(x){
#'     x <- format(x)
#'     tab <- vec2tab(x, dim(handy))
#'     message(
#'       paste(
#'         apply(tab, 1, paste, collapse = " "),
#'         collapse = " "
#'       )
#'     )
#'     message("
#' ", appendLF = F)
#'   })
#'   invisible()
#' }
#' # showSteps(out$steps)
#' 
#' 
#' 
#' 
#' 
#'
#'
#' 
#' }
#'
#' 
metropolis <- function(init, moves, iter = 1E3, burn = 1000, thin = 10,
  engine = c("Cpp","R")
){

  ## preliminary checking
  ##################################################
  engine <- match.arg(engine)

  if(engine == "R"){
  ## in R
  ##################################################

  nMoves <- ncol(moves)
  state  <- matrix(nrow = nrow(moves), ncol = iter)

  ## run burn-in

  current <- unname(init)
  unifs <- runif(burn)
  
  message("Running chain... ", appendLF = FALSE)
  
  for(k in 1:burn){

    move      <- sample(c(-1,1), 1) * moves[,sample(nMoves,1)]
    propState <- current + move
    
    if(any(propState < 0)){
      prob <- 0
    } else {
      prob <- exp( sum(lfactorial(current)) - sum(lfactorial(propState)) )
    }
    
    if(unifs[k] < prob) current <- propState # else current
    
  }
  state[,1] <- current


  ## run main sampler

  totalRuns <- 0
  probTotal <- 0
  unifs <- runif(iter*thin)  
  
  for(k in 2:iter){
  	
  	for(j in 1:thin){

      move      <- sample(c(-1,1), 1) * moves[,sample(nMoves,1)]
      propState <- current + move
    
      if(any(propState < 0)){
        prob <- 0
      } else {
        prob <- exp( sum(lfactorial(current)) - sum(lfactorial(propState)) )
      }
      probTotal <- probTotal + min(1, prob)

      if(unifs[k*(thin-1)+j] < prob) current <- propState # else current
      
      totalRuns <- totalRuns + 1        
    }

    state[,k] <- current    
  }
  message("done.")  
  
  ## format output
  out <- list(steps = state, moves = moves, 
    acceptProb = probTotal / totalRuns
  )
  
  
  

  } else if(engine == "Cpp"){


  ## in Cpp
  ##################################################
    
  current   <- unname(init)  
  allMoves  <- cbind(moves, -moves)  
  message("Running chain... ", appendLF = FALSE)  
  current   <- metropolisCpp(current, allMoves, burn, 1)$steps[,burn]
  out       <- metropolisCpp(current, allMoves, iter, thin)
  out$moves <- moves
  message("done.")

  }


  ## return output
  ##################################################  

  out[c("steps", "moves", "acceptProb")]
}

