#' Analyze a Rank Dataset
#'
#' \code{spectral} analyzes a rank dataset for order interactions; see examples for details.  
#' 
#' @param data a vector in the frequency format (see examples)
#' @param n the number of objects to select from
#' @param k the number of objects selected
#' @param levels the names of the outcomes, in the proper order
#' @param iter iterations in metropolis
#' @param burn burn in
#' @param thin thinning
#' @return a list containing named elements
#' \itemize{
#'   \item \code{effects}: the pure ith order effects as a data frame computed by projecting the data onto the isotypic subspaces.  the length of each is the same as the data, choose(n, k).
#'   \item \code{effectsNorms}: the l2 norms of each effect. 
#'   \item \code{statsMatrices}: the lower order statistics calculating matrices, made with Tmaker.
#'   \item \code{moves}: the markov moves for moving around each V, computed with markov on the statsMatrices.  only the positive moves are given.
#'   \item \code{samps}: the samples from each space conditioned on each level of statistics.  this is the output of metropolis.
#'   \item \code{obs}: a list of the observed data and its lower level summaries.
#'   \item \code{exp}: a list of the expected number of samples at each level given the summary statistics at the previous (lower) level.  these are computed from the samples from metropolis by (1) summarizing them with Tmaker and then (2) averaging the samples.  the expected V2 samples (for example), are determined by taking the samples with the same V1 statistics (in V3, say), summarizing them to V2 with Tmaker, and then averaging every cell across the samples.
#'   \item \code{fullExp}: this is the result of taking each of the samples with the same lower-order statistics and averaging them.  see exp, which is a reduction of fullExp.
#'   \item \code{residuals}: obs - exp
#'   \item \code{isotypicBases}: a list of the basis vectors of each isotypic subspace; computed as the eigenvalues of the result of Amaker, and grouped by eigenvalue.
#'   \item \code{sampsEffects}: the effects determined by each of the samples, projected onto the isotypic subspaces.
#'   \item \code{sampsEffectsNorms }: the norms of the effects of the samples.
#'   \item \code{sampsEffectsNormSummary }: a summary of the norms of the effects of the samples.
#'   \item \code{showStages}: a function that prints out the observed, expected, and residuals of sequentially conditioning on the sample size, first order statistics, second order statistics, and so on.
#'   \item \code{showFit}: a function that prints out a summary of the fit of the model.
#'   \item \code{decompose}: a function that takes a vector of the same length of the table given and summarizes it to its lower level statistics.
#'   \item \code{sampsDecomposed}: every sample decomposed.
#'   \item \code{statistic}: the pearson's chi-squared (X2), likelihood ratio (G2), Freeman-Tukey (FT), Cressie-Read (CR), and Neyman modified chi-squared (NM) statistics computed using the observed data (obs) and expected data (exp) at each level.
#'   \item \code{sampsStats}: the statistics computed on each of the samples.
#'   \item \code{p.value}: the exact p-values of individual tests, accurate to Monte-Carlo error.  these are computed as the proportion of samples with statistics equal to or larger than the oberved statistic.
#'   \item \code{p.value.se}: the standard errors of the p-values computed using the standard asymptotic formula of sqrt(p(1-p)/n).  a measure of the Monte-Carlo error.
#' }	

#' @export spectral
#' @examples
#'
#' \dontrun{
#'
#'
#'
#' ## voting statistics at different levels
#' ############################################################
#'
#' # load the cookies dataset:
#' data(cookie)
#' cookie$freq
#' cookie$cookies
#'
#'
#' # performing the spectral analysis
#' (out <- spectral(cookie$freq, 6, 3, cookie$cookies))
#'  
#'
#' out$obs # the original observations, and the summary statistics
#'
#' out$exp # each level is conditional on the previous level's statistics
#'         # (e.g. what you would expect for 1st order effects given sample size)
#'         # these are generated using 10k markov bases based mcmc samples
#'
#' out$p.value # these are approximate exact test p-values using various 
#'             # popular test statistics.  the approximations are good to
#'             # monte carlo error
#'
#' out$p.value.se # these are the standard errors using the sqrt(p*(1-p)/n)
#'                # asymptotic formula, known to have poor performance
#'                # for small/large p; see package binom for better
#'
#' out$statistic # the individual statistics are also available
#'               # the values are not comprable across Vi levels (the rows)
#'               # as they have asymptotic chi-squared distributions with 
#'               # different degrees of freedom
#'
#' out$fullExp # you can also get the expected number of samples at each scale
#'             # for tables with the same ith order statistics, i = 0, ..., k-1
#'
#'
#' # these can be seen to (re)construct an expected picture of the 
#' # complete data given each successive collection of statistics
#' cbind(
#'   obs = cookie$freq, 
#'   as.data.frame(lapply(out$fullExp, function(x) round(x[[4]],1)))
#' )[c(2:4,1)]
#' # notice that the reconstruction given only the first order statistics
#' # (the number of individual cookies selected) is quite good
#'
#' 
#' 
#' # instead of using the reconstructions from the exp coming from 
#' # the samples, you could reconstruct the summaries of the observed
#' # data using bump; it's not quite as good :
#' V0 <- bump(cookie$freq, 6, 3, 3, 0)
#' V1 <- bump(cookie$freq, 6, 3, 3, 1)
#' V2 <- bump(cookie$freq, 6, 3, 3, 2)
#' 
#' cbind(
#'   obs = cookie$freq, 
#'   round(data.frame(
#'     V0 = bump(V0, 6, 3, 0, 3),
#'     V1 = bump(V1, 6, 3, 1, 3),    
#'     V2 = bump(V2, 6, 3, 2, 3)
#'   ), 2)
#' )[c(2:4,1)]
#' 
#' 
#'
#'
#' # you can see the model step-by-step with showStages() :
#' out$showStages()
#' # notice (1) the significant reduction in the residuals after conditioning
#' # on the first order statistics and also (2) the powdery noise after 
#' # conditioning on the second order statistics.
#' # the p-values reflect the same: 
#' #   * the residuals from conditioning on the sample size show the first
#' #     order effects are strongly significant (in out$p.value V1 = 0)
#' #   * the residuals from conditioning on the first order effects suggest
#' #     the second order effects might be significant (V2 ~ .04-.13ish)
#' #   * the residuals from conditioning on the second order effects indicate
#' #     the third order effects are entirely insignificant (V3 > .2)
#' 
#' 
#' # the isotypic subpaces can be used to determine the pure order effects :
#'
#' out$isotypicBases # bases of the isotypic subspaces (here 4)
#'
#' out$effects # pure ith order effects; cookie$freq projected onto the bases
#'             # these are their effects at the data level, so they all have
#'             # the same length as the original dataset: choose(n, k)
#'
#' zapsmall(rowSums(out$effects)) # the effects sum to the data
#'
#'
#' # if the Vk effects are 0, then the conclusion is that Vk is perfectly
#' # predicted with the (k-1)st level statistics.  this may lead to the 
#' # conclusion that the l2 norms (say) of the effects might be used to 
#' # gauge the relative strength of effects :
#' out$effectsNorms # = apply(out$effects, 2, lpnorm)
#' 
#' 
#' # the natural (not full-dimensional) residuals can be seen with the summary
#' out
#' # or with
#' out$residuals
#' # these are the residuals (obs ith level stats) - (exp ith level stats)
#' # given the (i-1)st statistics
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
#' # bump is a useful function :
#' out$obs
#' bump(cookie$freq, 6, 3, 3, 0) # the 0 level is the number of voters, not votes
#' bump(cookie$freq, 6, 3, 3, 1) 
#' bump(cookie$freq, 6, 3, 3, 2)
#' bump(cookie$freq, 6, 3, 3, 3)
#' 
#' V1 <- out$obs$V1 # = bump(cookie$freq, 6, 3, 3, 1)
#' bump(V1, 6, 3, 1, 0)
#' bump(V1, 6, 3, 1, 1)
#' bump(V1, 6, 3, 1, 2) # cbind(bump(V1, 6, 3, 1, 2), out$exp$V2)
#' bump(V1, 6, 3, 1, 3) # cbind(bump(V1, 6, 3, 1, 3), out$fullExp$V1[[4]])
#' # the differences here are between an observation and an expectation
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
#' out$obs$V1 - out$exp$V1
#' out$residuals$V1
#' out$decompose(out$effects$V1)$V1
#' 
#' out$obs$V2 - out$exp$V2
#' out$residuals$V2
#' 
#'   out$decompose(out$effects$V0)$V2 + 
#'   out$decompose(out$effects$V1)$V2 +
#'   out$decompose(out$effects$V2)$V2 -
#'   out$exp$V2
#'   
#'
#'
#'
#' 
#' # this is how to reconstruct the observation given the effects
#' # the cols of out$effects are the Vk order effects reconstructed
#' # from the lower level effects
#' out$obs$V0
#' zapsmall(
#'   out$decompose(out$effects$V0)$V0
#' )
#' 
#' out$obs$V1
#' zapsmall(
#'   out$decompose(out$effects$V0)$V1 + 
#'   out$decompose(out$effects$V1)$V1
#' )
#' 
#' out$obs$V2
#' zapsmall(
#'   out$decompose(out$effects$V0)$V2 + 
#'   out$decompose(out$effects$V1)$V2 +
#'   out$decompose(out$effects$V2)$V2
#' )
#' 
#' out$obs$V3
#' zapsmall(
#'   out$decompose(out$effects$V0)$V3 + 
#'   out$decompose(out$effects$V1)$V3 +
#'   out$decompose(out$effects$V2)$V3 +
#'   out$decompose(out$effects$V3)$V3
#' )
#' zapsmall(rowSums(out$effects))
#' 
#' all(cookie$freq == zapsmall(rowSums(out$effects)))
#' 
#' 
#' 
#' out$effects$V0
#' out$effects$V0 + out$effects$V1
#' out$effects$V0 + out$effects$V2
#' out$effects$V0 + out$effects$V3
#' 
#' 
#' 
#' str(out$sampsDecomposed)
#' as.data.frame(lapply(out$sampsDecomposed, function(l) rowMeans(l$V3)))
#' 
#' eff0 <- rowMeans(out$sampsDecomposed$V0$V3)
#' cbind(eff0, out$effects$V0)
#' 
#' eff1 <- rowMeans(out$sampsDecomposed$V1$V3 - eff0)
#' cbind(eff1, out$effects$V1)
#' 
#' eff2 <- rowMeans(out$sampsDecomposed$V2$V3 - eff0 - eff1)
#' cbind(eff2, out$effects$V2)
#' 
#' sum(eff0)
#' sum(eff1)
#' sum(eff2)
#' 
#' 
#' 
#' 
#' 
#'
#'
#' str(out$sampsEffectsNorms)
#' 
#' data <- out$sampsEffectsNorms$V0$V3
#' plot(density(data))
#' curve(dnorm(x, mean(data), sd(data)), col = "red", add = TRUE)
#' 
#' data <- out$sampsEffectsNorms$V0$V2
#' plot(density(data))
#' curve(dnorm(x, mean(data), sd(data)), col = "red", add = TRUE)
#' 
#' data <- out$sampsEffectsNorms$V0$V1
#' plot(density(data))
#' curve(dnorm(x, mean(data), sd(data)), col = "red", add = TRUE)
#' 
#' 
#' data <- out$sampsEffectsNorms$V1$V3
#' plot(density(data))
#' curve(dnorm(x, mean(data), sd(data)), col = "red", add = TRUE)
#' 
#' data <- out$sampsEffectsNorms$V1$V2
#' plot(density(data))
#' curve(dnorm(x, mean(data), sd(data)), col = "red", add = TRUE)
#' 
#' 
#' data <- out$sampsEffectsNorms$V2$V3
#' plot(density(data))
#' curve(dnorm(x, mean(data), sd(data)), col = "red", add = TRUE)
#' 
#' 
#' 
#' 
#' 
#' 
#' ## how to convert data into the right format
#' ############################################################
#' # this essentially just uses some clever indexing tricks
#' # to reorder the data in the way you want
#' 
#' data <- cookie$raw       # an example raw, unordered dataset
#' levels <- cookie$cookies # the order of the objects you want
#' levsNndcs <- 1:length(levels)
#' names(levsNndcs) <- levels
#' 
#' 
#' # arrange selections within rows (order of selection doesn't matter)
#' data <- t(apply(data, 1, function(x) x[order(levsNndcs[x])] ))
#'
#' 
#' # arrange rows (order of selectors doesn't matter)
#' for(k in ncol(data):1) data <- data[order(levsNndcs[data[,k]]),]
#' 
#' 
#' # check that you've done the right thing
#' all( data == cookie$sorted )
#' 
#' # the data frequency order should match that of subsets:
#' subsets(levels, 1)
#'
#' subsets(levels, 2)
#' sapply(subsets(levels, 2), paste, collapse = ", ")
#' 
#' subsets(levels, 3)
#' sapply(subsets(levels, 3), paste, collapse = ", ")
#' 
#' names(cookie$freq)
#' names(cookie$freq) == sapply(subsets(levels, 3), paste, collapse = ", ")
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
#' 
#' 
#'
#' ## other examples
#' ############################################################
#' 
#' # rvotes provides uniform samples
#'
#' n <- 4
#' k <- 2
#' 
#' raw <- rvotes(250, n, k)
#' rawTogether <- apply(raw, 1, paste, collapse = " ")
#' levels <- sapply(subsets(n, k), paste, collapse = " ")
#' freq <- table( factor(rawTogether, levels = levels) )
#' (out <- spectral(freq, n, k))
#'
#' out$p.value
#' out$showStages()
#' 
#' out$obs
#' out$exp
#'
#'
#'
#'
#'
#' n <- 6
#' k <- 3
#' raw <- rvotes(250, n, k)
#' rawTogether <- apply(raw, 1, paste, collapse = " ")
#' levels <- sapply(subsets(n, k), paste, collapse = " ")
#' freq <- table( factor(rawTogether, levels = levels) )
#' (out <- spectral(freq, n, k))
#' 
#' 
#' n <- 7
#' k <- 3
#' raw <- rvotes(250, n, k)
#' rawTogether <- apply(raw, 1, paste, collapse = " ")
#' levels <- sapply(subsets(n, k), paste, collapse = " ")
#' freq <- table( factor(rawTogether, levels = levels) )
#' (out <- spectral(freq, n, k))
#'
#'
#'
#' n <- 8
#' k <- 3
#' raw <- rvotes(250, n, k)
#' rawTogether <- apply(raw, 1, paste, collapse = " ")
#' levels <- sapply(subsets(n, k), paste, collapse = " ")
#' freq <- table( factor(rawTogether, levels = levels) )
#' # out <- spectral(freq, n, k) # breaks
#'
#' 
#' 
#' 
#' }
spectral <- function(data, n, k, levels, iter = 1e4, burn = 1e3, thin = 10){
  
  if(missing(levels)) levels <- 1:n     
  if(!missing(levels) && length(levels) != n){
  	stop("if specified, levels must be of length n.", call. = FALSE)
  }
  
  
  # make list of moves
  message("Computing moves... ", appendLF = FALSE)
  listOfMoves <- as.list(rep(NA, k))
  names(listOfMoves) <- paste0("V", 0:(k-1))
  for(d in 0:(k-1)){
    listOfMoves[[d+1]] <- markov(Tmaker(n, k, d))
  }
  message("done.")



  # run metropolis-hastings
  samps <- as.list(rep(NA, k))
  names(samps) <- paste0("V", 0:(k-1))  
  for(d in 0:(k-1)){
    samps[[d+1]] <- metropolis(unname(data), listOfMoves[[d+1]],
      iter = iter, burn = burn, thin = thin
    )
  }
  
  
  

  # compute the A matrix
  A <- Amaker(n, k)
  
  # compute eigen stuff  
  e <- eigen(A)
  
  # clean eigen
  e$values <- zapsmall(e$values)
  e$vectors <- e$vectors
  
  # number of evals
  nEvals <- length(unique(e$values))  
  
  # split evecs based on eval
  listOfBasisVectors <- split(
    split(e$vectors, col(e$vectors)),
    e$values
  )
  
  # list from V0 to Vk
  listOfBasisVectors <- rev(listOfBasisVectors)
  
  # convert to matrices
  listOfBasisVectors <- lapply(listOfBasisVectors, function(x){
    if(is.list(x)) return(unname(Reduce(cbind, x)))
    t(t(x))
  })
  
  # assign names
  names(listOfBasisVectors) <- paste0("V", 0:k)
  
  
 
  
  # compute projections of data  
  qrs <- lapply(listOfBasisVectors, qr)  
  effects <- lapply(qrs, function(x) as.numeric(qr.fitted(x, data)) )
  effects <- as.data.frame(effects)
  row.names(effects) <- sapply(subsets(levels, k), paste, collapse = ", ")

  # norm projections of data
  effectsNorms <- sapply(effects, lpnorm)




  # compute projections of samps and their effectsNorms
  relative2Data <- sampsEffects <- sampsEffectsNorms <- as.list(rep(NA, k))
  names(sampsEffects) <- names(sampsEffectsNorms) <- paste0("V", 0:(k-1))
  ph <- as.list(rep(NA, k+1))  
  names(ph) <- paste0("V", 0:k)
  message("Projecting and norming... ", appendLF = FALSE)       

  for(d in 0:(k-1)){ # cycle through the same ith-order effects

  	sampsEffects[[d+1]] <- sampsEffectsNorms[[d+1]] <- relative2Data[[d+1]] <- ph

    for(dd in 0:k){ # cycle through dimensions to project onto    	

      # project      
      sampsEffects[[d+1]][[dd+1]] <- qr.fitted(qrs[[dd+1]], samps[[d+1]]$steps)        
      
      # norm
      sampsEffectsNorms[[d+1]][[dd+1]] <- apply(sampsEffects[[d+1]][[dd+1]], 2, lpnorm)      
      
      # compare to data
      relative2Data[[d+1]][[dd+1]] <- 
        mean( sampsEffectsNorms[[d+1]][[dd+1]] <= effectsNorms[dd+1] )
    
    }
    
  }
  message("done.")
  
  
  

  # summarize norms
  summarize <- function(v){
    if(length(v) == 1 && is.na(v)) return(NA)
    c(
        #"2.5%" = unname(quantile(v, .025)),
        "Mean" = mean(v),
       #"97.5%" = unname(quantile(v, .975)),
           "Sd" = sd(v)
    )
  }
  sampsEffectsNormsSummary <- lapply(sampsEffectsNorms, lapply, summarize)
  sampsEffectsNormsSummary <- lapply(sampsEffectsNormsSummary, as.data.frame)
  sampsEffectsNormsSummary <- lapply(sampsEffectsNormsSummary, function(df) df[,1:(k+1)])  
  relative2Data <- lapply(relative2Data, as.data.frame)
  for(i in 1:length(sampsEffectsNormsSummary)){
    sampsEffectsNormsSummary[[i]]   <- rbind(
      sampsEffectsNormsSummary[[i]], 
      unname(effectsNorms),      
      relative2Data[[i]]
    )
    row.names(sampsEffectsNormsSummary[[i]]) <- 
      c("Mean", "Sd", "data", "% <= data")
  }
  names(sampsEffectsNormsSummary) <- paste(names(sampsEffectsNormsSummary), "samples")




  # name listOfBasisVectors, make into matrices, zapsmall
  listOfBasisVectors <- lapply(listOfBasisVectors, function(x){
    if(is.list(x)) return(unname(Reduce(cbind, x)))
    t(t(x))
  })
  listOfBasisVectors <- lapply(listOfBasisVectors, zapsmall)
  

  
  
  # generate summary statistics
  observed <- as.list(rep(NA, k+1))
  for(d in 0:k){
  	observed[[d+1]] <- Tmaker(n, k, d) %*% unname(data)
    row.names(observed[[d+1]]) <- 
  	  sapply(subsets(levels, d), paste, collapse = ", ")
    colnames(observed[[d+1]]) <- "N"
  }
  names(observed) <- paste0("V", 0:k)
  
  

  # compute expected number of samples at each level
  # given the observed number at each level
  fullExpected <- as.list(rep(NA, k))
  names(fullExpected) <- paste0("V", 0:(k-1))
  Ts <- lapply(as.list(0:k), function(x) Tmaker(n, k, x)) # list of T's  
  for(d in 0:(k-1)){
  	fullExpected[[d+1]] <- lapply(Ts, function(Tmat){
  	  rowMeans(Tmat %*% samps[[d+1]]$steps)
  	})    
  }
  
  
  
  # summarize the above to the expected number of samples
  # at each level given the observed number of samples
  # at the previous level
  expected <- observed
  for(d in 1:k){
    expected[[d+1]] <- t(t(fullExpected[[d]][[d+1]]))
    row.names(expected[[d+1]]) <- 
      sapply(subsets(levels, d), paste, collapse = ", ")
    colnames(expected[[d+1]]) <- "N"      
  }
  
  
  
  
  # compute residuals
  residuals <- observed
  for(d in 1:k){
    residuals[[d+1]] <- observed[[d+1]] - expected[[d+1]]
  }  
  
  
  
  # showing stages
  showStages <- function(){
  	count <- 0
  	postScript <- c(
  	  "the sample size :\n",
  	  "first order statistics :\n",
  	  "second order statistics :\n",
  	  "third order statistics :\n",
  	  "fourth order statistics :\n",
  	  "fifth order statistics :\n",
  	  "sixth order statistics :\n"
  	)
    lapply(fullExpected, function(l){
      resMat <- obsMat <- expMat <- 
        matrix(numeric(0), nrow = choose(n,k), ncol = k+1)  
      obsMat[2,1] <- expMat[2,1] <- observed[[1]]          
      for(d in 0:k){
        expMat[1:choose(n, d), d+1] <- l[[d+1]]
        obsMat[1:choose(n, d), d+1] <- observed[[d+1]][,1]
        resMat[1:choose(n, d), d+1] <- observed[[d+1]][,1] - l[[d+1]]
      }
      mat <- cbind(obsMat, expMat, resMat)
      mat <- round(mat, 1)
      mat2 <- apply(mat, 2, format)
      mat2[is.na(mat)] <- ""
      mat2[1,0*(k+1)+1] <- "Obs: "
      mat2[1,1*(k+1)+1] <- "Exp: "    
      mat2[1,2*(k+1)+1] <- "Resid:"      
      mat2[3,0*(k+1)+1] <- paste0("(x", k, ")")
      mat2[3,1*(k+1)+1] <- paste0("(x", k, ")")      
      mat2 <- apply(mat2, 2, format)
      cat(paste("Conditioning on", postScript[count+1]))
      cat(apply(mat2, 1, paste, collapse = "  "), sep = "\n")      
      cat("\n")
      count <<- count + 1
    })  
    invisible()
  } 
  
  
  
  
  # showing expected
  showFit <- function(){
    resMat <- obsMat <- expMat <- 
      matrix(numeric(0), nrow = choose(n,k), ncol = k+1)  
    obsMat[2,1] <- expMat[2,1] <- observed[[1]]    
    for(d in 1:k){
      expMat[1:choose(n, d), d+1]  <- fullExpected[[d]][[d+1]]
      obsMat[1:choose(n, d), d+1] <- observed[[d+1]][,1]
      resMat[1:choose(n, d), d+1] <- residuals[[d+1]][,1]
    }
    mat <- cbind(obsMat, expMat, resMat)
    mat <- round(mat, 1)
    mat2 <- apply(mat, 2, format)
    mat2[is.na(mat)] <- ""
    mat2[1,0*(k+1)+1] <- "Obs: "
    mat2[1,1*(k+1)+1] <- "Exp: "    
    mat2[1,2*(k+1)+1] <- "Resid:"
    mat2[3,0*(k+1)+1] <- paste0("(x", k, ")")
    mat2[3,1*(k+1)+1] <- paste0("(x", k, ")")
    mat2 <- apply(mat2, 2, format)    
    cat(apply(mat2, 1, paste, collapse = "  "), sep = "\n")
    cat("\n")
    invisible()  
  } 
  
  
  
  
  # decompose an observation into the various order effects
  decompose <- function(v){
    out <- as.list(rep(NA, k+1))
    names(out) <- paste0("V", 0:k)
    for(d in 0:k){
      out[[d+1]] <- Ts[[d+1]] %*% unname(v)
      row.names(out[[d+1]]) <- 
        sapply(subsets(levels, d), paste, collapse = ", ")
      colnames(out[[d+1]]) <- "N"
    }
    out
  }   
    



  # decompose samps
  sampsDecomposed <- lapply(samps, function(x){
    out <- as.list(rep(NA, k+1))
    names(out) <- paste0("V", 0:k)  	
    for(d in 0:k) out[[d+1]] <- Ts[[d+1]] %*% x$steps
    out <- lapply(out, function(mat){
      matrix(as.integer(mat), nrow = nrow(mat))
    })
    out
  })
  #for(d in 1:k) sampsDecomposed[[d]] <- sampsDecomposed[[d]][[d+1]]
  #names(sampsDecomposed) <- paste0("V", 1:k)
  #lapply(sampsDecomposed, rowMeans) # = expected
  
  


  # compute statistics for each combo of observed/expected  
  X2 <- vector(length = k)
  names(X2) <- paste0("V", 1:k)
  PR <- G2 <- CR <- FT <- NM <- X2
  for(d in 1:k){
    PR[d] <- computeUProbsCpp(observed[[d+1]])  	
    X2[d] <- computeX2sCpp(observed[[d+1]], expected[[d+1]])
    G2[d] <- computeG2sCpp(observed[[d+1]], expected[[d+1]])    
    FT[d] <- computeCRsCpp(observed[[d+1]], expected[[d+1]], -.5)
    CR[d] <- computeCRsCpp(observed[[d+1]], expected[[d+1]], 2/3)    
    NM[d] <- computeNMsCpp(observed[[d+1]], expected[[d+1]])        
  }
  statistic <- rbind(PR, X2, G2, FT, CR, NM)
  
  
  
  
  
  # compute statistics for each combo of samps/expected  
  X2s <- as.list(rep(NA, length = k))
  names(X2s) <- paste0("V", 1:k)
  PRs <- G2s <- CRs <- FTs <- NMs <- X2s
  for(d in 1:k){
    PRs[[d]] <- computeUProbsCpp(sampsDecomposed[[d]][[d+1]])  	
    X2s[[d]] <- computeX2sCpp(sampsDecomposed[[d]][[d+1]], expected[[d+1]])
    G2s[[d]] <- computeG2sCpp(sampsDecomposed[[d]][[d+1]], expected[[d+1]])    
    FTs[[d]] <- computeCRsCpp(sampsDecomposed[[d]][[d+1]], expected[[d+1]], -.5)
    CRs[[d]] <- computeCRsCpp(sampsDecomposed[[d]][[d+1]], expected[[d+1]], 2/3)    
    NMs[[d]] <- computeNMsCpp(sampsDecomposed[[d]][[d+1]], expected[[d+1]])        
  }
  sampsStats <- list(PRs = PRs, X2s = X2s, G2s = G2s, FTs = FTs, CRs = CRs, NMs = NMs)

  
  
  
  # compute p values and their standard errors
  p.value <- matrix(nrow = nrow(statistic), ncol = k)
  colnames(p.value) <- paste0("V", 1:k)
  row.names(p.value) <- c("PR", "X2", "G2" ,"FT", "CR", "NM")
  mid.p.value <- p.value
  for(i in 1:nrow(statistic)){
    for(d in 1:k){
      p.value[i,d] <- mean( sampsStats[[i]][[d]] >= statistic[i,d] )
      mid.p.value[i,d] <- 
        mean( sampsStats[[i]][[d]] == statistic[i,d] )/2 +
        mean( sampsStats[[i]][[d]] > statistic[i,d] )
    }
  }
  p.value[1,] <- 1 - p.value[1,] 
  p.value.se <- sqrt(p.value*(1-p.value)/iter)
  mid.p.value.se <- sqrt(mid.p.value*(1-mid.p.value)/iter)  
  
    



  # return
  out <- list(
    effects = zapsmall(effects), effectsNorms = effectsNorms, 
    statsMatrices = Ts, moves = listOfMoves,
    samps = samps, 
    obs = observed, exp = expected, fullExp = fullExpected,
    residuals = residuals,    
    isotypicBases = listOfBasisVectors,    
    sampsEffects = sampsEffects, sampsEffectsNorms = sampsEffectsNorms,
    sampsEffectsNormsSummary = sampsEffectsNormsSummary,
    showStages = showStages, showFit = showFit,
    decompose = decompose, sampsDecomposed = sampsDecomposed,
    statistic = statistic, sampsStats = sampsStats,
    p.value = p.value, p.value.se = p.value.se,
    mid.p.value = mid.p.value, mid.p.value.se = mid.p.value.se,
    call = match.call()    
  )
  
  class(out) <- "spectral"
  out
}
