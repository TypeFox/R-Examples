optCliq <- function(cliqMat, tol = 10^-10, inner_loops = 100, outer_loops = 20){
  c_ans <- .Call('optCliq', cliqMat, tol, as.integer(inner_loops), as.integer(outer_loops), FALSE)
  ans <- cliqOptInfo(c_ans, tol)
  return(ans)
}

initICNPMLE <- function(L, R, Lin = NULL, Rin = NULL, A = NULL, 
                        max.inner = 100, max.outer = 100, tol = 1e-10){
    hm <- NULL
#    if(is.null(A) ){
      if(is.null(Lin))  Lin <- 1
      if(is.null(Rin))  Rin <- 1
    
      n <- length(L)
      if(n != length(R) ) stop('length(L) != length(R)')
    
      timeMatrix <- cbind(L, R, 0, 1)
      if(length(Lin) == 1 & length(Rin)) B <- c(Lin, Rin)
      else{
        if(length(Lin) != n | length(Rin) != n) stop('length of Lin or Rin not equal to 1 or n')
        B <- as.matrix( cbind(Lin, Rin, 1, 1) )
      }
    hm <- MLEcens::reduc(timeMatrix, B, cm = TRUE)
    A <- hm$cm
#    }
#  else
#  storage.mode(A) <- 'integer'
  cliqFit <- optCliq(A, tol, inner_loops = max.inner, outer_loops = max.outer)
  
  ans <- list(pf = cliqFit$pvec, intmap = t( hm$rects[,1:2]  ) )
  class(ans) <- 'icfit'
  return(ans)
}

ICNPMLE <- function(times, B = c(1,1), max.inner = 100, max.outer = 100, tol = 1e-10){
  times <- as.matrix(times)
  univariate = FALSE
  if(ncol(times) == 2){
    times <- cbind(times, 0, 1)
    univariate = TRUE
  }
  
  if(ncol(times) != 4){
    stop('number of columns of times not equal to 2 or 4')
  }
  hmInfo <- MLEcens::reduc(times, B = B, cm = T)
  cliqs <- hmInfo$cm
  
  optCliqInfo <- optCliq(cliqs, tol = tol, inner_loops = max.inner, outer_loops = max.outer)
  ans <- IC_NPMLE(optCliqInfo, hmInfo$rects, isUni = univariate, B = B)
  return(ans)
}


simBVCen <- function(n = 1000){
  t1 <- runif(n)
  t2 <- runif(n)
  
  l1 <- rep(0, n)
  l2 <- rep(0, n)
  r1 <- rep(1, n)
  r2 <- rep(1, n)
  
  c1.1 <- runif(n, 0, 2/3)
  c1.2 <- runif(n, c1.1, 1)
  
  btwn0_1.1 <- t1 < c1.1
  btwn1_2.1 <- c1.1 <= t1 & t1 < c1.2
  above2.1 <- !(btwn0_1.1 | btwn1_2.1)

  r1[btwn0_1.1] <- c1.1[btwn0_1.1]
    
  l1[btwn1_2.1] <- c1.1[btwn1_2.1]
  r1[btwn1_2.1] <- c1.2[btwn1_2.1]
  
  l1[above2.1] <- c1.2[above2.1]

  c2.1 <- runif(n, 0, 2/3)
  c2.2 <- runif(n, c2.1, 1)

  btwn0_1.2 <- t2 < c1.1
  btwn1_2.2 <- c2.1 <= t2 & t2 < c2.2
  above2.2 <- !(btwn0_1.2 | btwn1_2.2)
  
  r2[btwn0_1.2] <- c2.1[btwn0_1.2]
  
  l2[btwn1_2.2] <- c2.1[btwn1_2.2]
  r2[btwn1_2.2] <- c2.2[btwn1_2.2]
  
  l2[above2.2] <- c2.2[above2.2]
  
  ans <- data.frame(l1, r1, l2, r2)
  return(as.matrix(ans))
}

cliqOptInfo <- setRefClass('cliqOptInfo',
                           fields = c('pvec', 'llh', 'error', 'tot_iters', 'outer_iters', 'conv'),
                           methods = list(
                             show = function(){
                               cat('Clique Optimizer Object\n')
                               cat('Final log likelihood = ', llh, '\nNumeric Error = ',
                                   error, '\nTotal iterations = ', tot_iters, 
                                   '\nOuter iterations = ', outer_iters, '\n')
                             },
                             initialize = function(cList, tol){
                               pvec <<- cList[[1]]
                               llh <<- cList[[2]]
                               tot_iters <<- cList[[3]]
                               outer_iters <<- cList[[4]]
                               error <<- cList[[5]]
                               conv <<- error < tol
                               
                             }
                           ))

IC_NPMLE <- setRefClass('IC_NPMLE',
                         fields = c('p', 'rects', 'bounds', 'conv', 'llh', 'cliqOptInfo', 'isUni'),
                         methods = list(
                           show = function(){
                             if(!isUni) cat('Bivariate Interval Censored Optimization Results\n')
                             else       cat('Interval Censored Optimization Result\n')
                             cat('Final log likelihood = ', cliqOptInfo$llh, '\nNumeric Error = ',
                                 cliqOptInfo$error, '\nTotal iterations = ', cliqOptInfo$tot_iters, 
                                 '\nOuter iterations = ', cliqOptInfo$outer_iters, '\n')
                             cat('Number of maximal intersections = ', length(cliqOptInfo$pvec), 
                                 '\nNumber of maximal intersections w/ positive mass = ', 
                                 length(p), '\n')
                           },
                           initialize = function(cliqInfo, inRects, isUni, B){
                             cliqOptInfo <<- cliqInfo
                             isPos <- cliqInfo$pvec > 0
                             p <<- cliqInfo$pvec[isPos]
                             if(isUni){
                               inRects <-inRects[,1:2]
                             }
                             isUni <<- isUni
                             rects <<- inRects[isPos,]
                             llh <<- cliqOptInfo$llh
                             conv <<- cliqInfo$conv
                             bounds <<- B
                           }
                         ))