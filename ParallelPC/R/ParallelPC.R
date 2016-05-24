#########################################################################################
### Author: Thuc Duy Le, Tao Hoang and Shu Hu
### Date: 26th March 2015
### Please cite the following paper when using the code
### Paper: A fast PC algorithm for high dimensional causal discovery with multi-core PCs
#########################################################################################


if(getRversion() >= "2.15.1")  utils::globalVariables(c("gaussCItest","idaFast","detectCores","makeCluster","clusterEvalQ","stopCluster","pc.cons.intern","pdsep","udag2pag","udag2pdag","udag2pdagSpecial","udag2pdagRelaxed","rfci.vStruc","udag2apag","getNextSet","clusterApply","clusterSplit","mclapply","getNextSet","find.unsh.triple"))
#' Estimate (Initial) Skeleton of a DAG using the PC_stable Algorithm
#' @importFrom methods as new
#' @importFrom stats cor cov median qnorm
#' @importFrom utils read.csv
#' @description 
#' This is the skeleton (stable) function in the pcalg package. It is copied here to localise the parallel functions.
#' @param suffStat Sufficient statistics: List containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest Predefined function for testing conditional independence. The function is internally called as indepTest(x,y,S,suffStat) and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list containing all relevant elements for the conditional independence decisions. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha significance level (number in (0,1) for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param method Character string specifying method; the default, "stable" provides an order-independent skeleton, see 'Details' below.
#' @param m.max Maximal size of the conditioning sets that are considered in the conditional independence tests.
#' @param fixedGaps logical symmetric matrix of dimension p*p. If entry [i,j] is true, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges a logical symmetric matrix of dimension p*p. If entry [i,j] is true, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete logical needed for the case indepTest(*) returns NA. If it is true, the corresponding edge is deleted, otherwise not.
#' @param verbose if TRUE, detailed output is provided.
#' @return An object of class "pcAlgo" (see pcAlgo in the pcalg package) containing an estimate of the skeleton of the underlying DAG, the conditioning sets (sepset) that led to edge removals and several other parameters.
#' @examples
#' ##########################################
#' ## Using skeleton_stable
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' skeleton_stable(suffStat, indepTest=gaussCItest, p=p, method="stable", alpha=0.01)
#' @export
skeleton_stable <- function(suffStat, indepTest, alpha, labels, p,
                            method = c("stable", "original", "stable.fast"), m.max = Inf,
                            fixedGaps = NULL, fixedEdges = NULL,
                            NAdelete = TRUE, verbose = FALSE)
{ 
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
  #  if (method == "stable.fast" && .Platform$OS.type == "windows") {
  #    method <- "stable"
  #    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
  #  }
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (fixedEdges != t(fixedEdges))
    stop("fixedEdges must be symmetric")
  
  start_total <- proc.time()
  
  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose), 
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options, PACKAGE = "ParallelPC");
    G <- res$amat
    # sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1
  }
  else {
    ## Original R version
    
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remainingEdgeTests <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remainingEdgeTests,"\n",sep="")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remainingEdgeTests) {
        if (verbose && i%%100 == 0) cat("|i=", i, "|iMax=", nrow(ind), "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(pval >= alpha) { # independent
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }
  
  total_t = proc.time()-start_total
  
  # write results
 # cat('n=', suffStat[[2]], ',p=', p, '\n', sep="")
 # cat('Num CI Tests=', n.edgetests, ',Total CI Tests=', sum(unlist(n.edgetests)), ',Total Time=', total_t[3], '\n', sep=" ")
  
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }
  
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  
}## end{ skeleton }

#### Parallelized skeleton estimation ####
#' Estimate (Initial) Skeleton of a DAG. 
#' @importFrom methods as new
#' @description 
#' This is the parallelised version of the skeleton function in the pcalg package.
#' @param suffStat Sufficient statistics: List containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest Predefined function for testing conditional independence. The function is internally called as indepTest(x,y,S,suffStat) and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list containing all relevant elements for the conditional independence decisions. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha significance level (number in (0; 1) for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param method Character string specifying method; the default, "parallel" provides an efficient skeleton, see skeleton_parallel.
#' @param mem.efficient Uses less amount of memory at any time point while running the algorithm
#' @param workers Creates a set of copies of R running in parallel and communicating over sockets.
#' @param num_workers The numbers of cores CPU as numbers of workers to run the algorithm
#' @param m.max Maximal size of the conditioning sets that are considered in the conditional independence tests.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete  logical needed for the case indepTest(*) returns NA. If it is true, the corresponding edge is deleted, otherwise not.
#' @param verbose if TRUE, detailed output is provided.
#' @return An object of class "pcAlgo" (see pcAlgo in the pcalg package) containing an estimate of the skeleton of the underlying DAG, the conditioning sets (sepset) that led to edge removals and several other parameters.
#' @examples
#' ##########################################
#' ## Using skeleton_parallel without mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' skeleton_parallel(suffStat,indepTest=gaussCItest,p=p,method="parallel",alpha=0.01,num_workers=2)
#'
#' ##########################################
#' ## Using skeleton_parallel with mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' skeleton_parallel(suffStat,indepTest=gaussCItest,p=p,method="parallel",
#' alpha=0.01,num_workers=2,mem.efficient=TRUE)
#' @export
skeleton_parallel <- function(suffStat, indepTest, alpha, labels, p,
                              method = c("parallel"),  mem.efficient=FALSE, workers, num_workers,
                              m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                              NAdelete = TRUE, verbose = FALSE)
{
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
  } else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  } else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (fixedEdges != t(fixedEdges))
    stop("fixedEdges must be symmetric")
  
  pval <- NULL
  sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
  ## save maximal p value
  pMax <- matrix(-Inf, nrow = p, ncol = p)
  diag(pMax) <- 1
  done <- FALSE
  ord <- 0
  n.edgetests <- numeric(1)# final length = max { ord}
  
  # edge test function, conditioning on x's neighbours
  edge_test_xy <- function(x, y) {
    G_xy <- TRUE
    num_tests_xy <- 0
    pMax_xy <- pMax[x, y]
    sepset_xy <- NULL
    done_xy <- TRUE
    if (G_xy && !fixedEdges[y, x]) {
      nbrsBool <- G.l[[x]]
      nbrsBool[y] <- FALSE
      nbrs <- seq_p[nbrsBool]
      #rm(nbrsBool)
      length_nbrs <- length(nbrs)
      if (length_nbrs >= ord) {
        if (length_nbrs > ord) done_xy <- FALSE
        S <- seq_len(ord)
        repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
          num_tests_xy <- num_tests_xy + 1
          pval <- indepTest(x, y, nbrs[S], suffStat)
          if(is.na(pval)) pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
          if (pMax_xy < pval) pMax_xy <- pval
          if(pval >= alpha) { # independent
            G_xy <- FALSE
            sepset_xy <- nbrs[S]
            break
          } else {
            nextSet <- getNextSet(length_nbrs, ord, S)
            if (nextSet$wasLast)
              break
            S <- nextSet$nextSet
            #rm(nextSet)
          } ## if (pval >= alpha)
        } ## repeat
        #rm(S)
      } ## if (length_nbrs >= ord)
    } ## if(!done)
    list(G_xy, sepset_xy, num_tests_xy, pMax_xy, done_xy)
  }
  
  # edge test function
  edge_test <- function(i) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    num_tests_i <- 0
    G_i <- TRUE
    pMax_xy <- pMax[x, y]
    pMax_yx <- pMax[y, x]
    sepset_xy <- NULL
    sepset_yx <- NULL
    done_i <- TRUE
    
    # conditioning on neighbors of x
    res_x <- edge_test_xy(x, y)
    G_i <- res_x[[1]]
    sepset_xy <- res_x[[2]]
    num_tests_i <- num_tests_i + res_x[[3]]
    pMax_xy <- res_x[[4]]
    done_i <- done_i & res_x[[5]]
    
    if (G_i) {
      if (ord == 0) {
        num_tests_i <- num_tests_i + 1
      } else {
        # conditioning on neighbors of y
        res_y <- edge_test_xy(y, x)
        G_i <- res_y[[1]]
        sepset_yx <- res_y[[2]]
        num_tests_i <- num_tests_i + res_y[[3]]
        pMax_yx <- res_y[[4]]
        done_i <- done_i & res_y[[5]]
      }
    }
    
    # cleanup
    #rm(x)
    #rm(y)
    #rm(res_x)
    #rm(res_y)
    
    list(i, G_i, sepset_xy, sepset_yx, num_tests_i, pMax_xy, pMax_yx, done_i)
  }
  
  edge_tests <- function(l) {
    res <- vector("list",length(l))
    for (k in 1:length(l)) {
      res[[k]] <- edge_test(l[[k]])
    }
    res
  }
  
#   total_mem <- function() {
#     if (Sys.info()[["sysname"]] == "Linux") {
#       total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^MemFree:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Cached:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Inactive:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Buffers:' /proc/meminfo", intern = TRUE))))/1000
#       return(total)
#     } else if (Sys.info()[["sysname"]] == "Windows") {
#       #total <- as.numeric(memory.limit())
#       total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("wmic OS get FreePhysicalMemory /Value", intern=TRUE))[3]))/1000
#       return(total)
#     } else { # Mac OS X
#       total <- 4096*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages free'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages inactive'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages speculative'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages purgeable'", intern = TRUE))))/1000000
#       return(total)
#     }
#   }
  
total_mem <- function() {
  tryCatch({
    if (Sys.info()[["sysname"]] == "Linux") {
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^MemFree:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Cached:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Inactive:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Buffers:' /proc/meminfo", intern = TRUE))))/1000
      return(total)
    } else if (Sys.info()[["sysname"]] == "Windows") {
      #total <- as.numeric(memory.limit())
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("wmic OS get FreePhysicalMemory /Value", intern=TRUE))[3]))/1000
      return(total)
    } else if (Sys.info()[["sysname"]] == "Darwin") { # Mac OS X
      total <- 4096*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages free'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages inactive'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages speculative'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages purgeable'", intern = TRUE))))/1000000
      return(total)
    }
    else { # other OS, i.e. Solaris
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'free memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'inactive memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'buffer memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'swap cache'", intern = TRUE))))/1000
      return(total)
    }
  }, error=function(e){
    return(1024)
  }, warning=function(e){
    return(1024)
  })
}

  parallel_threshold <- 100
  if (mem.efficient) {
    mem_per_test <- 2 #MB
    tests_per_batch <- as.integer(total_mem() / mem_per_test)
  }
  
  start_total <- proc.time()
  
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord+1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ## For comparison with C++ sort according to first row
    ind <- ind[order(ind[, 1]), ]
    ## Consider only unique edge
    ind <- subset(ind, ind[, 1] < ind[, 2])
    remainingEdgeTests <- nrow(ind)
    ## Order-independent version: Compute the adjacency sets for any vertex
    ## Then don't update when edges are deleted
    G.l <- split(G, gl(p,p))
    
    if (!mem.efficient) {
      tests_per_batch <- remainingEdgeTests
    }
    
    for (j in seq(1, remainingEdgeTests, by=tests_per_batch)) {
      l <- min(remainingEdgeTests, j + tests_per_batch - 1)
      if (l - j + 1 < num_workers) {
        num_workers <- l - j + 1
      }
      res <- NULL
      if (l - j + 1 < parallel_threshold) {
        res <- lapply(j:l, edge_test)
      } else if (Sys.info()[['sysname']] == 'Windows') {
        res <- do.call("c", clusterApply(workers, clusterSplit(workers, j:l), edge_tests))
      } else {
        res <- mclapply(j:l, edge_test, mc.cores=num_workers, mc.set.seed=FALSE, mc.cleanup=TRUE, mc.allow.recursive=FALSE)
      }
      
      # synchronize
      for (p_obj in res) {
        i <- p_obj[[1]]
        x <- ind[i, 1]
        y <- ind[i, 2]
        n.edgetests[ord1] <- n.edgetests[ord1] + p_obj[[5]]
        pMax[x, y] <- p_obj[[6]]
        pMax[y, x] <- p_obj[[7]]
        G[x, y] <- G[y, x] <- p_obj[[2]]
        if (!p_obj[[2]]) {
          if (!is.null(p_obj[[3]])) sepset[[x]][[y]] <- p_obj[[3]]
          if (!is.null(p_obj[[4]])) sepset[[y]][[x]] <- p_obj[[4]]
        }
        done <- done & p_obj[[8]]
      }
    }
    
    # increase the nbrs size
    ord <- ord + 1
  } ## while()
  
  total_t = proc.time()-start_total
  
  # write results
  cat('n=', suffStat[[2]], ',p=', p, '\n', sep="")
  cat('Num CI Tests=', n.edgetests, ',Total CI Tests=', sum(unlist(n.edgetests)), ',Total Time=', total_t[3], '\n', sep=" ")
  
  for (i in 1:(p - 1)) {
    for (j in 2:p)
      pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
  }
  
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }
  
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  
}## end{ skeleton }

#########################################################################################
### pc_stable and pc_parallel
#########################################################################################


### This function is exactly the same as pc() in the pcalg package ###
#' Estimate the Equivalence Class of a DAG using the PC_stable Algorithm
#' @importFrom methods as new
#' @description 
#' Estimate the equivalence class of a directed acyclic graph (DAG) from observational data, using the PC_stable algorithm.
#' @param suffStat A list of sufficient statistics, containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest A function for testing conditional independence. It is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha significance level (number in (0,1) for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param verbose If TRUE, detailed output is provided.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max Maximal size of the conditioning sets that are considered in the conditional independence tests.
#' @param u2pd String specifying the method for dealing with conflicting information when trying to orient edges (see pcalg for details).
#' @param skel.method Character string specifying method; the default, "stable" provides an order-independent skeleton.
#' @param conservative Logical indicating if the conservative PC is used. In this case, only option u2pd = "relaxed" is supported. See pcalg for more information.
#' @param maj.rule Logical indicating that the triples shall be checked for ambiguity using a majority rule idea, which is less strict than the conservative PC algorithm. For more information, see the pcalg package.
#' @param solve.confl If TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations. In this case, only option u2pd = relaxed is supported. Note, that therefore the resulting object might not be a CPDAG because bi-directed edges might be present. See details for more information.
#' @return An object of class "pcAlgo" (see pcAlgo in the pcalg package) containing an estimate of the equivalence class of the underlying DAG.
#' @examples
#' ##########################################
#' ## Using pc_stable
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' pc_stable(suffStat, indepTest=gaussCItest, p=p, skel.method="stable", alpha=0.01)
#' @export
pc_stable <- function(suffStat, indepTest, alpha, labels, p,
                      fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
                      u2pd = c("relaxed", "rand", "retry"),
                      skel.method = c("stable", "original", "stable.fast"),
                      conservative = FALSE, maj.rule = FALSE,
                      solve.confl = FALSE, verbose = FALSE)
{ 
  
  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  
  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")
  
  ## Skeleton
  skel <- skeleton_stable(suffStat, indepTest, alpha, labels=labels, method = skel.method,
                          fixedGaps=fixedGaps, fixedEdges=fixedEdges,
                          NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  
  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose=verbose))
  }
  else { ## u2pd "relaxed" : conservative _or_ maj.rule
    
    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf=c(2,1), maj.rule=maj.rule, verbose=verbose)
    udag2pdagRelaxed(pc.$sk, verbose=verbose,
                     unfVect=pc.$unfTripl)
  }
}##{pc_stable}

############## Parallel-PC algorithm (based on pc() from pcalg package) #################
#' Estimate the Equivalence Class of a DAG  using the PC_parallel Algorithm
#' @importFrom methods as new
#' @description 
#' Estimate the equivalence class of a directed acyclic graph (DAG) from observational data, using the PC_parallel algorithm.
#' @param suffStat A list of sufficient statistics, containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest A function for testing conditional independence. It is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha significance level (number in (0,1) for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param verbose If TRUE, detailed output is provided.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max Maximal size of the conditioning sets that are considered in the conditional independence tests.
#' @param u2pd String specifying the method for dealing with conflicting information when trying to orient edges (see pcalg for details).
#' @param mem.efficient If TRUE, uses less amount of memory at any time point while running the algorithm.
#' @param skel.method Character string specifying method; the default, "parallel",  skeleton_parallel for learning the causal structure.
#' @param conservative Logical indicating if the conservative PC is used. In this case, only option u2pd = "relaxed" is supported. Note that therefore the resulting object might not be extendable to a DAG. See pcalg for details.
#' @param maj.rule Logical indicating that the triples shall be checked for ambiguity using a majority rule idea, which is less strict than the conservative PC algorithm. For more information, see pcalg.
#' @param solve.confl If TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations.See pcalg for details.
#' @param num.cores The numbers of cores CPU to run the algorithm.
#' @return An object of class "pcAlgo" (see pcAlgo in the pcalg package) containing an estimate of the equivalence class of the underlying DAG.
#' @examples
#' ##########################################
#' ## Using pc_parallel without mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' pc_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", alpha=0.01, num.cores=2)
#' 
#' ##########################################
#' ## Using pc_parallel with mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' pc_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", 
#' alpha=0.01, num.cores=2, mem.efficient=TRUE)
#' 
#' #################################################
#' ## Using pc_parallel with mutual information test
#' #################################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' #The first parameter is the dataset rather than suffStat
#' pc_parallel(gmG$x, indepTest=mig, p=p, skel.method="parallel", 
#' alpha=0.01, num.cores=2, mem.efficient=TRUE)
#' 
#' @export
pc_parallel <- function(suffStat, indepTest, alpha, labels, p,
                        fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
                        u2pd = c("relaxed", "rand", "retry"),
                        skel.method = c("parallel"), mem.efficient=FALSE,
                        conservative = FALSE, maj.rule = FALSE,
                        solve.confl = FALSE, verbose = FALSE, num.cores = detectCores())
{
  
  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  
  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")
  
  # prepare the workers
  num_workers <- num.cores
  if (num_workers < 2) {
    stop("The number of cores is insufficient to run parallel-PC")
  }
  workers <- NULL
  if (Sys.info()[['sysname']] == 'Windows') {
    workers <- makeCluster(num_workers, type="PSOCK")
    eval(suffStat)
    clusterEvalQ(workers, library(pcalg))
  }
  
  ## Skeleton
  skel <- skeleton_parallel(suffStat, indepTest, alpha, labels=labels, method = skel.method, workers=workers, num_workers=num_workers,
                            fixedGaps=fixedGaps, fixedEdges=fixedEdges, mem.efficient=mem.efficient,
                            NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  
  # stop workers
  if (Sys.info()[['sysname']] == 'Windows') {
    stopCluster(workers)
  }
  
  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose=verbose))
  } else { ## u2pd "relaxed" : conservative _or_ maj.rule
    
    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf=c(2,1), maj.rule=maj.rule, verbose=verbose)
    udag2pdagRelaxed(pc.$sk, verbose=verbose,
                     unfVect=pc.$unfTripl)
  }
}##{pc_parallel}

#########################################################################################
### fci_stable and fci_parallel
#########################################################################################

### This function is exactly the same as the FCI_Stable in the pcalg package ###
#' Estimate a PAG, using the FCI_stable algorithm
#' @importFrom methods as new
#' @description 
#' This is the FCI stable version in the pcalg package.
#' @param suffStat Sufficient statistics: List containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest Predefined function for testing conditional independence. The function is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list containing all relevant elements for the conditional independence decisions. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha Significance level for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param skel.method Character string specifying method; the default, "stable", provides an order-independent skeleton, see skeleton.
#' @param type Character string specifying the version of the FCI algorithm to be used. By default, it is "normal", and so the normal FCI algorithm is called. If set to "anytime", the 'Anytime FCI' is called and m.max needs to be specified. If set to "adaptive", the 'Adaptive Anytime FCI' is called and m.max is not used. For more information, see the FCI function in the pcalg package.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max Maximum size of the conditioning sets that are considered in the conditional independence tests.
#' @param pdsep.max Maximum size of Possible-D-SEP for which subsets are considered as conditioning sets in the conditional independence tests. See pcalg for more details.
#' @param rules Logical vector of length 10 indicating which rules should be used when directing edges. See pcalg for more details.
#' @param doPdsep If TRUE, Possible-D-SEP is computed for all nodes, and all subsets of Possible-D-SEP are considered as conditioning sets in the conditional independence tests, if not defined otherwise in pdsep.max. If FALSE, Possible-D-SEP is not computed, so that the algorithm simplifies to the Modified PC algorithm of Spirtes, Glymour and Scheines (2000, p.84).
#' @param biCC If TRUE, only nodes on paths between nodes x and y are considered to be in Possible-D-SEP(x) when testing independence between x and y. Uses biconnected components, biConnComp from RBGL.
#' @param conservative Logical indicating if the unshielded triples should be checked for ambiguity the second time when v-structures are determined.
#' @param maj.rule Logical indicating if the unshielded triples should be checked for ambiguity the second time when v-structures are determined using a majority rule idea, which is less strict than the standard conservative. For more information, see details.
#' @param verbose If true, more detailed output is provided.
#' @return An object of class fciAlgo (see fciAlgo in the pcalg package) containing the estimated graph (in the form of an adjacency matrix with various possible edge marks), the conditioning sets that lead to edge removals (sepset) and several other parameters.
#' @examples
#' ##########################################
#' ## Using fci_stable
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' fci_stable(suffStat, indepTest=gaussCItest, p=p, skel.method="stable", alpha=0.01)
#' @export
#' @references
#' 1. Diego Colombo, Marloes H Maathuis, Markus Kalisch, Thomas S Richardson, et al. Learning high-dimensional directed acyclic graphs with latent and selection variables. The Annals of Statistics, 40(1):294-321, 2012.
#' 
#' 2. Markus Kalisch, Martin Machler, Diego Colombo, Marloes H Maathuis, and Peter Buhlmann. Causal inference using graphical models with the r package pcalg.
#'  Journal of Statistical Software, 47(11):1-26, 2012.
                                        
                                        
fci_stable <- function(suffStat, indepTest, alpha, labels, p,
                       skel.method = c("stable", "original", "stable.fast"),
                       type = c("normal", "anytime", "adaptive"),
                       fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                       m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                       doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                       maj.rule = FALSE, verbose = FALSE)
{  
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")
  
  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")
  
  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")
  
  skel <- skeleton_stable(suffStat, indepTest, alpha, labels=labels, method = skel.method,
                          fixedGaps=fixedGaps, fixedEdges=fixedEdges,
                          NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL
  
  if (doPdsep) {
    if (verbose) cat("\nCompute PDSEP\n=============\n")
    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
                            alpha=alpha, version.unf = c(1,1),
                            maj.rule=FALSE, verbose=verbose)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- pdsep(skel@graph, suffStat, indepTest=indepTest, p=p,
                      sepset = pc.ci$sk@sepset, alpha=alpha, pMax=pMax,
                      m.max = if (type == "adaptive") max.ordSKEL else m.max,
                      pdsep.max=pdsep.max, NAdelete=NAdelete,
                      unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                      biCC=biCC, verbose=verbose)
    
    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                            verbose=verbose, version.unf = c(1, 1),
                            maj.rule=maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                                verbose=verbose, version.unf = c(2, 1),
                                maj.rule=maj.rule)
      tripleList <- nopdsep$unfTripl
      ##update the sepsets
      sepset <- nopdsep$sk@sepset
    }
  }
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
                  verbose = verbose)
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset=sepset, pMax=pMax, allPdsep=allPdsep)
  
}## {fci_stable}

############## Parallel-FCI algorithm (based on fci() from pcalg package) #################
#' Estimate a Partial Ancestral Graph using the FCI_parallel algorithm
#' @importFrom methods as new
#' @description 
#' Estimate a Partial Ancestral Graph (PAG) from observational data, using the FCI_parallel Algorithm. This is the parallelised version of the FCI algorithm in the pcalg package.
#' The parameters are consistent with the FCI algorithm in pcalg, except the parameter num.cores for specifying the number of cores CPU.
#' @param suffStat Sufficient statistics: List containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest Predefined function for testing conditional independence. The function is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list containing all relevant elements for the conditional independence decisions. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha Significance level for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param skel.method Character string specifying method; the default, "parallel", uses the parallelised method to build the skeleton of the graph, see skeleton_parallel.
#' @param mem.efficient Uses less amount of memory at any time point while running the algorithm.
#' @param type Character string specifying the version of the FCI algorithm to be used. By default, it is "normal", and so the normal FCI algorithm is called. If set to "anytime", the 'Anytime FCI' is called and m.max needs to be specified. If set to "adaptive", the 'Adaptive Anytime FCI' is called and m.max is not used. For more information, see Details.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max Maximum size of the conditioning sets that are considered in the conditional independence tests.
#' @param pdsep.max Maximum size of Possible-D-SEP for which subsets are considered as conditioning sets in the conditional independence tests. See pcalg for more details.
#' @param rules Logical vector of length 10 indicating which rules should be used when directing edges. See pcalg for more details.
#' @param doPdsep If TRUE, Possible-D-SEP is computed for all nodes, and all subsets of Possible-D-SEP are considered as conditioning sets in the conditional independence tests, if not defined otherwise in pdsep.max. If FALSE, Possible-D-SEP is not computed, so that the algorithm simplifies to the Modified PC algorithm of Spirtes, Glymour and Scheines (2000, p.84).
#' @param biCC If TRUE, only nodes on paths between nodes x and y are considered to be in Possible-D-SEP(x) when testing independence between x and y. Uses biconnected components, biConnComp from RBGL.
#' @param conservative Logical indicating if the unshielded triples should be checked for ambiguity the second time when v-structures are determined. 
#' @param maj.rule Logical indicating if the unshielded triples should be checked for ambiguity the second time when v-structures are determined using a majority rule idea, which is less strict than the standard conservative. For more information, see details.
#' @param verbose If true, more detailed output is provided.
#' @param num.cores Numbers of cores CPU to run the algorithm
#' @return An object of class fciAlgo (see fciAlgo in the pcalg package) containing the estimated graph (in the form of an adjacency matrix with various possible edge marks), the conditioning sets that lead to edge removals (sepset) and several other parameters.
#' @examples
#' ##########################################
#' ## Using fci_parallel without mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' fci_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", alpha=0.01, num.cores=2)
#'
#' ##########################################
#' ## Using fci_parallel with mem.efficeient
#' ##########################################
#' 
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' fci_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", 
#' alpha=0.01, num.cores=2, mem.efficient=TRUE)
#' 
#' #################################################
#' ## Using fci_parallel with mutual information test
#' #################################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' #' # The first parameter is the dataset
#' fci_parallel(gmG$x, indepTest=mig, p=p, skel.method="parallel", 
#' alpha=0.01, num.cores=2, mem.efficient=TRUE)
#' @export
#' @references
#' 1. Diego Colombo, Marloes H Maathuis, Markus Kalisch, Thomas S Richardson, et al. Learning high-dimensional directed acyclic graphs with latent and selection variables. The Annals of Statistics, 40(1):294-321, 2012.
#' 
#' 2. Markus Kalisch, Martin Machler, Diego Colombo, Marloes H Maathuis, and Peter Buhlmann. Causal inference using graphical models with the r package pcalg.
#'  Journal of Statistical Software, 47(11):1-26, 2012.
fci_parallel <- function(suffStat, indepTest, alpha, labels, p,
                         skel.method = c("parallel"),mem.efficient=FALSE,
                         type = c("normal", "anytime", "adaptive"),
                         fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                         m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                         doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                         maj.rule = FALSE, verbose = FALSE, num.cores = detectCores())
{
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")
  
  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")
  
  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")
  
  # prepare the workers
  num_workers <- num.cores
  if (num_workers < 2) {
    stop("The number of cores is insufficient to run parallel-PC")
  }
  workers <- NULL
  if (Sys.info()[['sysname']] == 'Windows') {
    workers <- makeCluster(num_workers, type="PSOCK")
    eval(suffStat)
    clusterEvalQ(workers, library(pcalg))
  }
  
  ## Skeleton
  skel <- skeleton_parallel(suffStat, indepTest, alpha, labels=labels, method = skel.method, workers=workers, num_workers=num_workers,
                            fixedGaps=fixedGaps, fixedEdges=fixedEdges, mem.efficient=mem.efficient,
                            NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  
  # stop workers
  if (Sys.info()[['sysname']] == 'Windows') {
    stopCluster(workers)
  }
  
  #
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL
  
  if (doPdsep) {
    if (verbose) cat("\nCompute PDSEP\n=============\n")
    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
                            alpha=alpha, version.unf = c(1,1),
                            maj.rule=FALSE, verbose=verbose)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- pdsep(skel@graph, suffStat, indepTest=indepTest, p=p,
                      sepset = pc.ci$sk@sepset, alpha=alpha, pMax=pMax,
                      m.max = if (type == "adaptive") max.ordSKEL else m.max,
                      pdsep.max=pdsep.max, NAdelete=NAdelete,
                      unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                      biCC=biCC, verbose=verbose)
    
    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                            verbose=verbose, version.unf = c(1, 1),
                            maj.rule=maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                                verbose=verbose, version.unf = c(2, 1),
                                maj.rule=maj.rule)
      tripleList <- nopdsep$unfTripl
      ##update the sepsets
      sepset <- nopdsep$sk@sepset
    }
  }
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
                  verbose = verbose)
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset=sepset, pMax=pMax, allPdsep=allPdsep)
  
}##{fci_parallel}

#########################################################################################
### rfci_stable and rfci_parallel
#########################################################################################

### This function is exactly the same as rfci() in the pcalg package ###
#' Estimate a PAG using the RFCI_stable Algorithm
#' @importFrom methods as new
#' @description 
#' This is the RFCI stable version in the pcalg package.
#' @param suffStat Sufficient statistics: List containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest Predefined function for testing conditional independence. The function is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list containing all relevant elements for the conditional independence decisions. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha significance level (number in (0,1) for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param skel.method Character string specifying method; the default, "stable" provides an order-independent skeleton, see skeleton.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max Maximum size of the conditioning sets that are considered in the conditional independence tests.
#' @param rules Logical vector of length 10 indicating which rules should be used when directing edges. See the pcalg package for details.
#' @param conservative Logical indicating if the unshielded triples should be checked for ambiguity after the skeleton has been found, similar to the conservative PC algorithm.
#' @param maj.rule Logical indicating if the unshielded triples should be checked for ambiguity after the skeleton has been found using a majority rule idea, which is less strict than the conservative.
#' @param verbose If true, more detailed output is provided.
#' @return An object of class fciAlgo (see fciAlgo in the pcalg package) containing the estimated graph (in the form of an adjacency matrix with various possible edge marks), the conditioning sets that lead to edge removals (sepset) and several other parameters.
#' @examples
#' ##########################################
#' ## Using rfci_stable
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' rfci_stable(suffStat, indepTest=gaussCItest, p=p, skel.method="stable", alpha=0.01)
#' @export
#' @references
#' 1. Diego Colombo, Marloes H Maathuis, Markus Kalisch, Thomas S Richardson, et al. Learning high-dimensional directed acyclic graphs with latent and selection variables. The Annals of Statistics, 40(1):294-321, 2012.
#' 
#' 2. Markus Kalisch, Martin Machler, Diego Colombo, Marloes H Maathuis, and Peter Buhlmann. Causal inference using graphical models with the r package pcalg.
#'  Journal of Statistical Software, 47(11):1-26, 2012.
rfci_stable <- function(suffStat, indepTest, alpha, labels, p,
                        skel.method = c("stable", "original", "stable.fast"),
                        fixedGaps = NULL, fixedEdges = NULL,
                        NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10),
                        conservative = FALSE, maj.rule = FALSE,
                        verbose = FALSE)
{
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  if (conservative && maj.rule)
    stop("Can only choose one of conservative or majority rule RFCI")
  if (verbose) cat("Compute Skeleton\n================\n")
  
  skel <- skeleton_stable(suffStat, indepTest, alpha, labels=labels, method = skel.method,
                          fixedGaps=fixedGaps, fixedEdges=fixedEdges,
                          NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  sk.A <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  ##the list of all ordered unshielded triples (the graph g does not change it is just a search!)
  u.t <- find.unsh.triple(sk.A, check=FALSE)
  
  ## check and orient v-structures recursively
  r.v. <- rfci.vStruc(suffStat, indepTest, alpha, sepset, sk.A,
                      unshTripl = u.t$unshTripl, unshVect = u.t$unshVect,
                      conservative = (conservative || maj.rule),
                      version.unf=c(1,1), maj.rule=maj.rule, verbose=verbose)
  A <- r.v.$amat
  sepset <- r.v.$sepset
  
  ## orient as many edge marks as possible
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), "\n")
  
  res <- udag2apag(A, suffStat, indepTest, alpha, sepset,
                   rules=rules, unfVect = r.v.$unfTripl, verbose=verbose)
  Amat <- res$graph
  colnames(Amat) <- rownames(Amat) <- labels
  new("fciAlgo", amat = Amat, call = cl, n = integer(0),
      max.ord = as.integer(skel@max.ord), max.ordPDSEP = 0L,
      n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0,
      sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", p))
  
}## {rfci_stable}

############## Parallel-RFCI algorithm (based on rfci() from pcalg package) #################
#' Estimate a PAG fast using the RFCI_parallel Algorithm
#' @importFrom methods as new
#' @description 
#' This is the parallelised version of the RFCI algorithm in the pcalg package.
#' @param suffStat Sufficient statistics: List containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest Predefined function for testing conditional independence. The function is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list containing all relevant elements for the conditional independence decisions. The return value of indepTest is the p-value of the test for conditional independence.
#' @param alpha Significance level for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param skel.method Character string specifying method; the default, "parallel" provides an efficient skeleton, see skeleton_parallel.
#' @param mem.efficient Uses less amount of memory at any time point while running the algorithm
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max Maximum size of the conditioning sets that are considered in the conditional independence tests.
#' @param rules Logical vector of length 10 indicating which rules should be used when directing edges. See the pcalg package for details.
#' @param conservative Logical indicating if the unshielded triples should be checked for ambiguity the second time when v-structures are determined. For more information, see pcalg.
#' @param maj.rule Logical indicating if the unshielded triples should be checked for ambiguity the second time when v-structures are determined using a majority rule idea, which is less strict than the standard conservative. For more information, see pcalg.
#' @param verbose If true, more detailed output is provided.
#' @param num.cores The numbers of cores CPU to run the algorithm
#' @return An object of class fciAlgo (see fciAlgo in the pcalg package) containing the estimated graph (in the form of an adjacency matrix with various possible edge marks), the conditioning sets that lead to edge removals (sepset) and several other parameters.
#' @examples
#' ##########################################
#' ## Using rfci_parallel without mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' rfci_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", alpha=0.01, num.cores=2)
#' 
#' ##########################################
#' ## Using rfci_parallel with mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel) 
#' data("gmG")
#' p<-ncol(gmG$x)
#' suffStat<-list(C=cor(gmG$x),n=nrow(gmG$x))
#' rfci_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", 
#' alpha=0.01, num.cores=2, mem.efficient=TRUE)
#' 
#' #################################################
#' ## Using fci_parallel with mutual information test
#' #################################################
#' library(pcalg)
#' library(parallel)
#' data("gmG")
#' p<-ncol(gmG$x)
#' 
#' # The first parameter is the dataset
#' rfci_parallel(gmG$x, indepTest=mig, p=p, skel.method="parallel", 
#' alpha=0.01, num.cores=2, mem.efficient=TRUE)
#' @export
#' @references
#' 1. Diego Colombo, Marloes H Maathuis, Markus Kalisch, Thomas S Richardson, et al. Learning high-dimensional directed acyclic graphs with latent and selection variables. The Annals of Statistics, 40(1):294-321, 2012.
#' 
#' 2. Markus Kalisch, Martin Machler, Diego Colombo, Marloes H Maathuis, and Peter Buhlmann. Causal inference using graphical models with the r package pcalg.
#'  Journal of Statistical Software, 47(11):1-26, 2012.
rfci_parallel <- function(suffStat, indepTest, alpha, labels, p,
                          skel.method = c("parallel"),mem.efficient=FALSE,
                          fixedGaps = NULL, fixedEdges = NULL,
                          NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10),
                          conservative = FALSE, maj.rule = FALSE,
                          verbose = FALSE, num.cores = detectCores())
{
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  if (conservative && maj.rule)
    stop("Can only choose one of conservative or majority rule RFCI")
  if (verbose) cat("Compute Skeleton\n================\n")
  
  # prepare the workers
  num_workers <- num.cores
  if (num_workers < 2) {
    stop("The number of cores is insufficient to run parallel-rfci")
  }
  workers <- NULL
  if (Sys.info()[['sysname']] == 'Windows') {
    workers <- makeCluster(num_workers, type="PSOCK")
    eval(suffStat)
    clusterEvalQ(workers, library(pcalg))
  }
  
  ## Skeleton
  skel <- skeleton_parallel(suffStat, indepTest, alpha, labels=labels, method = skel.method, workers=workers, num_workers=num_workers,
                            fixedGaps=fixedGaps, fixedEdges=fixedEdges, mem.efficient=mem.efficient,
                            NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  sk.A <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  ##the list of all ordered unshielded triples (the graph g does not change it is just a search!)
  u.t <- find.unsh.triple(sk.A, check=FALSE)
  
  ## check and orient v-structures recursively
  r.v. <- rfci.vStruc(suffStat, indepTest, alpha, sepset, sk.A,
                      unshTripl = u.t$unshTripl, unshVect = u.t$unshVect,
                      conservative = (conservative || maj.rule),
                      version.unf=c(1,1), maj.rule=maj.rule, verbose=verbose)
  A <- r.v.$amat
  sepset <- r.v.$sepset
  
  # stop workers
  if (Sys.info()[['sysname']] == 'Windows') {
    stopCluster(workers)
  }
  
  #
  
  ## orient as many edge marks as possible
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), "\n")
  
  res <- udag2apag(A, suffStat, indepTest, alpha, sepset,
                   rules=rules, unfVect = r.v.$unfTripl, verbose=verbose)
  Amat <- res$graph
  colnames(Amat) <- rownames(Amat) <- labels
  new("fciAlgo", amat = Amat, call = cl, n = integer(0),
      max.ord = as.integer(skel@max.ord), max.ordPDSEP = 0L,
      n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0,
      sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", p))
  
}## {rfci_parallel}

#########################################################################################
### IDA_stable and IDA_parallel
#########################################################################################

#' Estimate  Total Causal Effects
#' 
#' @description 
#' This the stable version (using stable-PC for structure learning) of the IDA algorithm in the pcalg package.
#' @param datacsv The dataset in csv format with rows are samples and columns are variables
#' @param cause The number of integer positions of the cause variables in the dataset
#' @param effect The number of integer positions of the target variables  in the dataset.
#' @param pcmethod Character string specifying method; the default, "stable", provides an order-independent skeleton. See Colombo, 2014.
#' @param alpha significance level (number in (0; 1) for the individual conditional independence tests.
#' @return A matrix that shows the causal effects (minimum of all possible effects) of the causes (columns) on the effects (rows).
#' @examples
#' ##########################################
#' ## Using IDA_stable
#' ##########################################
#' library(pcalg)
#' data("gmI")
#' datacsv <- cov(gmI$x)
#' IDA_stable(datacsv,1:2,3:4,"stable",0.01) 
#' @export
#' @references
#' 1. Marloes H Maathuis, Markus Kalisch, Peter Buhlmann, et al. Estimating high-dimensional
#' intervention effects from observational data. The Annals of Statistics, 37(6A):3133-3164,2009.
#' 
#' 2. Diego Colombo and Marloes H Maathuis. Order-independent constraint-based causal structure learning. The Journal of Machine Learning Research, 15(1):3741-3782, 2014.

IDA_stable <- function(datacsv, cause, effect, pcmethod, alpha){
  #result=IDA("EMT-35.csv", 1:35, 36:end, "stable", 0.01)
  
  
  if(is.character(datacsv)){
    data=read.csv(datacsv)
    #data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
  } else {
    data=datacsv #Assume there is no samplenames column and this is a data.frame.
  }  						#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
  data=scale(data) #standardise the data
  #print(data[1:5,])
  allnames=colnames(data)
  causenames=allnames[cause]
  effectnames=allnames[effect]
  
  multiset=character(0)
  result=matrix(nrow=length(effect), ncol=length(cause))
  suffStat=list(C=cor(data), n=nrow(data))
  indepTest=gaussCItest
  
  start_total_ida <- proc.time()
  pcFit <- pc_stable(suffStat, indepTest, p=ncol(data), alpha=alpha, skel.method=pcmethod)
  #pcFit <-pc_parallel(suffStat, indepTest=gaussCItest, p=ncol(data), skel.method=pcmethod, alpha=alpha, num.cores=4)
  for (l in cause){
    
    #Inferring causal effects
    caef<-idaFast(l,effect,cov(data), pcFit@graph )
    
    #min of absolute values.
    caef1<-matrix(nrow=length(effect),ncol=1)
    for (k in 1:length(effect)){
      caefabs<-abs(caef)
      index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
      pos<-index[1,2]
      caef1[k,]<-caef[k,pos]
    }
    result[,l]<-caef1
  }
  #total_t_ida = proc.time()-start_total_ida
  #cat('Total Time ida=', total_t_ida[3], '\n', sep=" ")
  colnames(result)=causenames
  rownames(result)=effectnames
  return(result)	
}##{IDA_stable}

#' Estimate  Total Causal Effects  using the IDA_parallel Algorithm
#' 
#' @description 
#' This is the parallelised version of the IDA (stable) algorithm in the pcalg package.
#' 
#' @param datacsv The dataset in csv format.
#' @param cause The number of integer positions of the cause variables in the dataset.
#' @param effect The number of integer  positions of the target variables  in the dataset.
#' @param pcmethod Character string specifying method; the default, "parallel", will use the parallelised method for learning the skeleton of the graph, see skeleton_parallel.
#' @param alpha significance level (number in (0; 1) for the individual conditional independence tests.
#' @param num.cores The numbers of cores CPU to run the algorithm
#' @param mem.efficient If TRUE, uses less amount of memory at any time point while running the algorithm
#' @return A matrix that shows the causal effects (minimum of all possible effects) of the causes (columns) on the effects (rows)
#' @examples
#' ##########################################
#' ## Using IDA_parallel without mem.efficeient
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' library(parallel)
#' data("gmI")
#' datacsv <- cov(gmI$x)
#' IDA_parallel(datacsv,1:2,3:4,"parallel",0.01, 2)
#' 
#' ##########################################
#' ## Using IDA_parallel with mem.efficeient
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' library(parallel)
#' data("gmI")
#' datacsv <- cov(gmI$x)
#' IDA_parallel(datacsv,1:2,3:4,"parallel",0.01, 2, TRUE)
#'
#' @export
#' @references
#' Marloes H Maathuis, Markus Kalisch, Peter Buhlmann, et al. Estimating high-dimensional
#' intervention effects from observational data. The Annals of Statistics, 37(6A):3133-3164,2009.
IDA_parallel <- function(datacsv, cause, effect, pcmethod, alpha, num.cores, mem.efficient=FALSE){
  #result=IDA("EMT-35.csv", 1:35, 36:end, "parallel",0.01, 4, TRUE)
  
  if(is.character(datacsv)){
    data=read.csv(datacsv)
    #data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
  } else {
    data=datacsv #Assume there is no samplenames column and this is a data.frame.
  }  						#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
  data=scale(data) #standardise the data
  #print(data[1:5,])
  allnames=colnames(data)
  causenames=allnames[cause]
  effectnames=allnames[effect]
  
  multiset=character(0)
  result=matrix(nrow=length(effect), ncol=length(cause))
  suffStat=list(C=cor(data), n=nrow(data))
  indepTest=gaussCItest
  
  start_total_ida <- proc.time()
  #pcFit <- pc_stable(suffStat, indepTest, p=ncol(data), alpha=alpha, skel.method=pcmethod)
  #pcFit <-pc_parallel(suffStat, indepTest=gaussCItest, p=ncol(data), skel.method=pcmethod, alpha=alpha, num.cores=num.cores, mem.efficient=mem.efficient)
  pcFit <-pc_parallel(suffStat, indepTest, p=ncol(data), skel.method=pcmethod, alpha=alpha, num.cores=num.cores, mem.efficient=mem.efficient)
  for (l in cause){
    
    #Inferring causal effects
    caef<-idaFast(l,effect,cov(data), pcFit@graph )
    
    #min of absolute values.
    caef1<-matrix(nrow=length(effect),ncol=1)
    for (k in 1:length(effect)){
      caefabs<-abs(caef)
      index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
      pos<-index[1,2]
      caef1[k,]<-caef[k,pos]
    }
    result[,l]<-caef1
  }
  #total_t_ida = proc.time()-start_total_ida
  #cat('Total Time ida=', total_t_ida[3], '\n', sep=" ")
  colnames(result)=causenames
  rownames(result)=effectnames
  return(result)	
}##{IDA_parallel}


################################################################

##Different indepTest methods

################################################################

###############################################################
#gaussCItest=bnlearn::zf
#' Gaussian conditional independence test
#' 
#' @description 
#' Gaussian conditional independence test. See the zf function in the bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat the data matrix with rows are samples and columns are the variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using zf
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' zf(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
zf=function(x,y,S, suffStat){
  
  #x, y, and S are passed from the main algorithm. We need to specify the paramater to 
  #satify the pre-defined tests
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStatis the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="zf")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}

###############################################################
#' The Monte Carlo permutation test for Gaussian conditional independence test
#' 
#' @description 
#' The Monte Carlo permutation test for Gaussian conditional independence test. See the mc-zf function in the bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The dataset in matrix format with rows are samples and columns are variables.
#' @return the p-value of the test.
#' @examples
#' ##########################################
#' ## Using mczf
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' mczf(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
mczf=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="mc-zf")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}
##############################################################
#' The sequential Monte Carlo permutation test for Gaussian conditional independence test.
#' 
#' @description 
#'The sequential Monte Carlo permutation test for Gaussian conditional independence test. See the smc-zf function in the bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The data matrix with rows are samples and columns are variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using smczf
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' smczf(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
smczf=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="smc-zf")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}
################################################################

################################################################

#mutual information: bnlearn::mi-g
#' Mutual information test
#' 
#' @description 
#' Mutual information test. See function mi-g in bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The dataset in matrix format with rows are samples and columns are variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using mig
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' mig(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 

mig=function(x,y,S, suffStat){
  
  #x, y, and S are passed from the main algorithm. We need to specify the paramater to 
  #satify the pre-defined tests
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStatis the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="mi-g")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}

################################################################
#mutual information Monte Carlo: bnlearn::mc-mi-g
#' The Monte Carlo permutation test (mc-mi-g)
#' 
#' @description 
#' The Monte Carlo permutation test for mutual information. See bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The dataset in matrix format with rows are samples and columns are variables.
#' @return the p-value of the test.
#' @examples
#' ##########################################
#' ## Using mcmig
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' mcmig(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
mcmig=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="mc-mi-g")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}

#mutual information Sequential Monte Carlo: bnlearn:: smc-mi-g
################################################################
#' The sequential Monte Carlo permutation test (smc-mi-g)
#' 
#' @description 
#' The sequential Monte Carlo permutation test. See bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The data matrix with rows are samples and columns are variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using smcmig
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' smcmig(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
smcmig=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="smc-mi-g")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}

#Shrinkage estimator for mutual information
################################################################
#' Shrinkage estimator for the mutual information (mi-g-sh)
#' 
#' @description 
#' Shrinkage estimator for the mutual information. See bnlearn package for more details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The dataset in matrix format with rows are samples and columns are variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using migsh
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' migsh(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
migsh=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="mi-g-sh")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}
#Pearson's chi-square
################################################################
#' The Pearson's correlation test
#' 
#' @description 
#' Linear correlation: Pearson's linear correlation test. 
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat the dataset with rows are samples and columns are variables.
#' @return the p-value of the test.
#' @examples
#' ##########################################
#' ## Using cor2 as a conditional independence test
#' ##########################################
#' library(pcalg)
#' library(bnlearn)
#' data("gmG")
#' suffStat<-gmG$x
#' cor2(1,2,3,suffStat)
#' ##Use cor2 with a causal discovery algorithm, e.g. PC
#' pc_stable(gmG$x, indepTest=cor2, p=ncol(gmG$x), alpha=0.01)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
cor2=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="cor")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}
################################################################
#Monte Carlo Pearson's chi-square
#' The Monte Carlo permutation test (mc-cor)
#' 
#' @description 
#' The Monte Carlo permutation test for Pearson's chi-square. See bnlearn package for details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The dataset in matrix format with rows are samples and columns are variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using mccor
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' mccor(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
mccor=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="mc-cor")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}
################################################################
#Sequential Monte Carlo Pearson's chi-square
#' The sequential Monte Carlo permutation test (smc-cor)
#' 
#' @description 
#' The sequential Monte Carlo permutation test. See bnlearn package for details.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of
#'        the remaining nodes. x, y, S all are integers, corresponding to variable or node
#'        numbers.
#' @param suffStat The dataset in matrix format with rows are samples and columns are variables.
#' @return The p-value of the test.
#' @examples
#' ##########################################
#' ## Using smccor
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' data("gmG")
#' suffStat<-gmG$x
#' smccor(1,2,3,suffStat)
#' @export
#' @references
#' Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. 
smccor=function(x,y,S, suffStat){
  
  if(!is.data.frame(suffStat)) suffStat=data.frame(suffStat)#suffStat is the dataset
  Namex=colnames(suffStat)[x]
  Namey=colnames(suffStat)[y]
  Namez=colnames(suffStat)[S]
  
  test=bnlearn::ci.test(x=Namex, y=Namey, z=Namez, data=suffStat, test="smc-cor")
  pval=test$p.value
  cat("x=",x," y=",y," S=", S, "pvalue=", pval, "\n")
  pval
  
  
}





#' Estimate  Total Causal Effects  of Joint Interventions 
#' 
#' @description 
#' This is the parallelised version of the jointIDA (stable) algorithm in the pcalg package.
#' 
#' @param datacsv The dataset in the csv format with rows are samples and columns are the variables.
#' @param cause The number of integer positions of the intervention variables in the dataset.
#' @param effect the  integer position of the target variable in the dataset.
#' @param method the method of calculating the final effect from multiple possible effects, e.g. min, max, median
#' @param pcmethod Character string specifying the method of the PC algorithm, e.g. stable for stable-PC, and parallel for parallel-PC.
#' @param alpha significance level (number in (0; 1) for the conditional independence tests.
#' @param num.cores  The numbers of cores CPU to run the algorithm
#' @param mem.efficient If TRUE, uses less amount of memory at any time point while running the algorithm
#' @param technique The character string specifying the technique that will be used to estimate the total joint causal effects in the pcalg package. 
#' RRC for Recursive regression for causal effects
#' MCD for Modifying the Cholesky decomposition
#' @return A matrix that shows the direct causal effects (minimum of all possible effects) of the (first) cause (columns) on the effects (rows)
# @examples
# ##########################################
# ## Using IDA_parallel without mem.efficeient
# ##########################################
# library(bnlearn)
# library(pcalg)
# library(parallel)
# data("gmI")
# datacsv <- cov(gmI$x)
# jointIDA_parallel(datacsv,1:2,3:4,method="min", pcmethod="parallel",0.01, 2, technique="RRC")
# 
# ##########################################
# ## Using IDA_parallel with mem.efficeient
# ##########################################
# library(bnlearn)
# library(pcalg)
# library(parallel)
# data("gmI")
# datacsv <- cov(gmI$x)
# jointIDA_direct(datacsv,1:2,3:4,method="min", pcmethod="parallel",0.01, 2, TRUE, technique="RRC")

jointIDA_direct=function(datacsv, cause, effect,method=c("min","max","median"), pcmethod="stable", alpha, num.cores=1, mem.efficient=FALSE,technique = c("RRC","MCD")){
  if(is.character(datacsv)){
    data=read.csv(datacsv)
    #data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
  } else {
    data=datacsv #Assume there is no samplenames column and this is a data.frame.
  }      				#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
  data=scale(data) #standardise the data
  #print(data[1:5,])
  allnames=colnames(data)
  causenames=allnames[cause]
  effectnames=allnames[effect]
  
  multiset=character(0)
  result=matrix(nrow=length(effect), ncol=length(cause))
  suffStat=list(C=cor(data), n=nrow(data))
  indepTest=gaussCItest
  
  start_total_jointida <- proc.time()
  
  if (pcmethod == "stable"){
    pcFit <- pc_stable(suffStat, indepTest, p=ncol(data), alpha=alpha, skel.method=pcmethod)
  }else {
    pcFit <- pc_parallel(suffStat, indepTest=gaussCItest, p=ncol(data), skel.method=pcmethod, alpha=alpha, num.cores=num.cores, mem.efficient=mem.efficient)
  }
  #pcFit<-pc(suffStat, indepTest=gaussCItest, p=ncol(data),alpha=alpha)
  #return(jointIda(cause,effect,cov(data),pcFit@graph,technique=technique))
  for(k in 1:length(effect)){
    
    caef<-pcalg::jointIda(cause,effect[k],cov(data),pcFit@graph,technique=technique)
    
    caefabs<-abs(caef)
    
    for(l in 1:length(cause)){
      
      if(method=="min"||method=="max"){
        if(method=="min"){
         # cat("hushu")
          index<-which(caefabs==min(caefabs[l,],na.rm = TRUE), arr.ind=TRUE)
        }else{
       #   cat("hushu1")
          index<-which(caefabs==max(caefabs[l,],na.rm = TRUE), arr.ind=TRUE)
        }
       # cat("index",index,"\n")
      #  cat("index[1,2]",index[1,2],"\n")
        
        pos<-index[1,2]
        
        result[k,l]<-caef[l,pos]
      }else if(method=="median"){
        result[k,l]<-median(caef[l,],na.rm = TRUE)
        
      }
      
    }
    
  }
  total_t_jointida = proc.time()-start_total_jointida
  #cat('Total Time jointida=', total_t_jointida[3], '\n', sep=" ")
  colnames(result)=causenames
  rownames(result)=effectnames
  return(result)
  
}#jointIDA



#' Estimate  Total Causal Effects  of Joint Interventions 
#' 
#' @description 
#' This is the parallelised version of the IDA (stable) algorithm in the pcalg package.
#' @param datacsv The dataset in csv format with rows are samples and columns are variables.
#' @param cause The number of integer positions of the intervention variables in the dataset.
#' @param effect the integer position of the target variable in the dataset.
#' @param pcmethod Character string specifying the method of the PC algorithm, e.g. stable for stable-PC, and parallel for parallel-PC.
#' @param alpha significance level (number in (0; 1) for the conditional independence tests.
#' @param num.cores  The numbers of cores CPU to run the algorithm
#' @param mem.efficient If TRUE, uses less amount of memory at any time point while running the algorithm
#' @param technique The character string specifying the technique that will be used to estimate the total joint causal effects in the pcalg package. 
#' RRC for Recursive regression for causal effects
#' MCD for Modifying the Cholesky decomposition
#' @return A matrix that shows the causal effects  of the causes (rows) on the effect. Different columns show different possible causal effect values.
#' @examples
#' ##########################################
#' ## Using IDA_parallel without mem.efficeient
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' library(parallel)
#' data("gmI")
#' datacsv <- cov(gmI$x)
#' jointIDA_parallel(datacsv,1:2,3, pcmethod="parallel",0.01, 2, technique="RRC")
#' 
#' ##########################################
#' ## Using IDA_parallel with mem.efficeient
#' ##########################################
#' library(bnlearn)
#' library(pcalg)
#' library(parallel)
#' data("gmI")
#' datacsv <- cov(gmI$x)
#' jointIDA_parallel(datacsv,1:2,3, pcmethod="parallel",0.01, 2, TRUE, technique="RRC")
#' @export
jointIDA_parallel=function(datacsv, cause, effect, pcmethod="stable", alpha, num.cores=1, mem.efficient=FALSE,technique = c("RRC","MCD")){
  if(is.character(datacsv)){
    data=read.csv(datacsv)
    #data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
  } else {
    data=datacsv #Assume there is no samplenames column and this is a data.frame.
  }        			#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
  data=scale(data) #standardise the data
  #print(data[1:5,])
  allnames=colnames(data)
  causenames=allnames[cause]
  effectnames=allnames[effect]
  
  #multiset=character(0)
  #result=matrix(nrow=length(effect), ncol=length(cause))
  suffStat=list(C=cor(data), n=nrow(data))
  indepTest=gaussCItest
  
  start_total_jointida <- proc.time()
  
  if (pcmethod == "stable"){
    pcFit <- pc_stable(suffStat, indepTest, p=ncol(data), alpha=alpha, skel.method=pcmethod)
  }else {
    pcFit <- pc_parallel(suffStat, indepTest=gaussCItest, p=ncol(data), skel.method=pcmethod, alpha=alpha, num.cores=num.cores, mem.efficient=mem.efficient)
  }
  #pcFit<-pc(suffStat, indepTest=gaussCItest, p=ncol(data),alpha=alpha)
  #return(jointIda(cause,effect,cov(data),pcFit@graph,technique=technique))
  #for(k in 1:length(effect)){
    
    result<-pcalg::jointIda(cause,effect,cov(data),pcFit@graph,technique=technique)
    
    
  total_t_jointida = proc.time()-start_total_jointida
  #cat('Total Time jointida=', total_t_jointida[3], '\n', sep=" ")
 # colnames(result)=causenames
  rownames(result)=causenames
  return(result)
  
}#jointIDA
############## Parallel-pcSelect algorithm (based on pcSelect() from pcalg package) #################
#' Estimate subgraph around a response variable using pcSelect_parallel.
#' @importFrom methods as new
#' @description 
#' This is the parallelised version of the pcSelect (stable) function in the pcalg package. Assume that we have a fixed target variable, the algorithm will test the
#' dependency between each variable and the target variable conditioning on combinations of other variables.
#' @param y The target (response) variable.
#' @param dm Data matrix with rows are samples and columns are variables.
#' @param method Character string specifying method; the default, "parallel" provides an parallelised method to implement all the conditional independence tests.
#' @param alpha Significance level of individual partial correlation tests.
#' @param corMethod "standard" or "Qn" for standard or robust correlation estimation
#' @param verbose Logical or in \{0,1,2\};
#' 
#'        FALSE, 0: No output,
#'        
#'        TRUE, 1: Little output,
#'        
#'        2: Detailed output.
#'        
#'        Note that such output makes the function very much slower.
#' @param directed Logical; should the output graph be directed?
#' @param num_workers The numbers of cores CPU to run the algorithm
#' @param mem.efficient If TRUE, uses less amount of memory at any time point while running the algorithm
#' @return 
#' G    A logical vector indicating which column of dm is associated with y.
#' 
#' zMin   The minimal z-values when testing partial correlations between y and each column of dm. The larger the number, the more consistent is the edge with the data.
#' @examples
#' ##########################################
#' ## Using pcSelect_parallel without mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' p <- 10
#' set.seed(101)
#' myDAG <- randomDAG(p, prob = 0.2)
#' n <- 1000
#' d.mat <- rmvDAG(n, myDAG, errDist = "normal")
#' pcSelect_parallel(d.mat[,10],d.mat[,-10], alpha=0.05,num_workers=2)
#' 
#' ##########################################
#' ## Using pcSelelct_parallel with mem.efficeient
#' ##########################################
#' library(pcalg)
#' library(parallel) 
#' p <- 10
#' set.seed(101)
#' myDAG <- randomDAG(p, prob = 0.2)
#' n <- 1000
#' d.mat <- rmvDAG(n, myDAG, errDist = "normal")
#' pcSelect_parallel(d.mat[,10],d.mat[,-10], alpha=0.05,mem.efficient=TRUE,num_workers=2)
#' @export

pcSelect_parallel <- function(y,dm,method = c("parallel"),
                              mem.efficient=FALSE,num_workers, alpha,
                              corMethod = "standard", verbose = FALSE, directed=FALSE)
{
  ## Purpose: Find columns in dm, that have nonzero parcor with y given
  ## any other set of columns in dm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - y: Response Vector (length(y)=nrow(dm))
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - verbose: 0-no output, 1-small output, 2-details
  ## ----------------------------------------------------------------------
  ## Value: List
  ## - G: boolean vector with connected nodes
  ## - zMin: Minimal z values
  ## ----------------------------------------------------------------------
  ## Author: LiangZHANG, Date: 25.08.15
  print("----------in the pc_select_parallel myself's----------")
  time0 <- Sys.time()
  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  print("vNms = ")
  print(vNms)
  cl <- match.call()
  
  zMin <- c(0,rep.int(Inf,p))
  C <- pcalg::mcor(cbind(y,dm), method = corMethod)
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- c(FALSE,rep.int(TRUE,p))
  seq_p <- 1:(p+1)
  done <- FALSE
  ord <- 0
  
  # prepare the workers
  if (num_workers < 2) {
    stop("The number of cores is insufficient to run parallel-PC")
  }
  workers <- NULL
  time1 <- Sys.time()
  if (Sys.info()[['sysname']] == 'Windows') {
    workers <- makeCluster(num_workers, type="PSOCK")
    eval(dm)
    clusterEvalQ(workers, library(pcalg))
  } 
  # the upper function cost a lot time!!! nearly 5 secs. why ??
  print("parallel prepare cost time:")
  print(Sys.time()-time1)
  time2 <- Sys.time()
  print("parallel perpare finished.")
  
  
  # edge test function, conditioning on x's neighbours
  edge_test_xy <- function(x,y){
    if (verbose>=1){
      cat("in_fun edge_test_xy,x=",x,"y=",y)
    }
    G_x <- TRUE
    G_y <- TRUE
    num_tests_xy <- 0
    zMin_x <- zMin[x]
    done_xy <- TRUE
    if (G_x){
      #nbrsBool <- G.l[[x]]
      nbrsBool <-as.logical(G.l)
      nbrsBool[x] <- FALSE
      nbrs <- seq_p[nbrsBool]
      ## neighbors of y without itself and x
      length_nbrs <- length(nbrs)
      if (length_nbrs >= ord){
        if (length_nbrs > ord) done_xy <- FALSE
        S <- seq(length = ord)
        
        ## now includes special cases (ord == 0) or (length_nbrs == 1):
        repeat{
          #cat("num_tests_xy=",num_tests_xy,"\n")
          num_tests_xy <- num_tests_xy + 1
          z <- pcalg::zStat(x,y, nbrs[S], C,n)
          #cat('x=',x,'y=',y,"\n")
          if(abs(z)<zMin_x) zMin_x <- abs(z)
          #cat('zMin_x=',zMin_x,"\n")
          if (verbose >= 2){
            #sprintf("x = %d",x)
            print(paste0("here,x = ",x))
            print(paste0("here,y = ",y))
            cat(paste("x:",vNms[x-1],"y:",(ytmp <- round((p+1)/2)),"S:"),
                c(ytmp,vNms)[nbrs[S]],paste("z:",z,"\n"))
          }
          if (abs(z) <= cutoff) {
            G_x <- FALSE
            #cat('do1',"\n")
            break
          }
          else {
            nextSet <- getNextSet(length_nbrs, ord, S)
            if(nextSet$wasLast){
              #cat('do2',"\n")
              break
            }
            S <- nextSet$nextSet
            #cat('do3',"\n")
          }
        } ## {repeat}
      } ## if (length_nbrs >= ord)
    } ## if(!done)
    list(G_x, num_tests_xy, zMin_x, done_xy)
  }
  
  # edge test function
  edge_test <- function(i) { # i for output actually is X
    y <- 1
    x <- ind[i]
    #print(paste0("-----in edge test no xy-----",i, x))
    if (verbose>=1){cat("-----in edge test no xy-----i, x =",i,x,"\n")}
    num_tests_i <- 0
    G_i <- TRUE
    zMin_x <- zMin[x]
    #zMin_y <- zMin[y]
    done_i <- TRUE
    # conditioning on neighbors of x
    if (verbose>=1){
      print("-----in edge_test(x, y)-----")
    }
    res_x <- edge_test_xy(x, y)   #########################
    #print("-----the res_x of the edge_test_xy-----")
    #print(res_x)
    G_i <- res_x[[1]]
    #cat("G_i=",G_i,"\n")
    #sepset_xy <- res_x[[2]]
    num_tests_i <- num_tests_i + res_x[[2]]
    # num_tests_i <- num_tests_i + 1
    #cat("num_tests_i_1=",num_tests_i,"\n")
    zMin_x <- res_x[[3]]
    done_i <- done_i & res_x[[4]]
    ####################### cut
    #     if (G_i) {
    #       if (ord == 0) {
    #         #         num_tests_i <- num_tests_i+1
    #       }else {
    #         #         # conditioning on neighbors of y
    #         print("-----in edge_test(y, x)-----")
    #         res_y <- edge_test_xy(y, x)              
    #         #         G_i <- res_y[[1]]
    #         #         num_tests_i <- num_tests_i + res_y[[2]]
    #         #         num_tests_i <- num_tests_i + 1
    #         #        zMin_y <- res_y[[3]]
    #         done_i <- done_i & res_y[[4]]
    #       }
    #     }
    ####################### cut
    # rm(x)
    # rm(y)
    # rm(res_x)
    # rm(res_y)
    list(i, G_i, num_tests_i, zMin_x, done_i)
  }
  ###############################################
  edge_tests <- function(l) {
    res <- vector("list",length(l))
    for (k in 1:length(l)) {
      res[[k]] <- edge_test(l[[k]])
    }
    res
  }
  ###############################################
#   total_mem <- function() {
#     if (Sys.info()[["sysname"]] == "Linux") {
#       total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^MemFree:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Cached:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Inactive:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Buffers:' /proc/meminfo", intern = TRUE))))/1000
#       return(total)
#     } else if (Sys.info()[["sysname"]] == "Windows") {
#       #total <- as.numeric(memory.limit())
#       total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("wmic OS get FreePhysicalMemory /Value", intern=TRUE))[3]))/1000
#       return(total)
#     } else { # Mac OS X
#       total <- 4096*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages free'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages inactive'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages speculative'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages purgeable'", intern = TRUE))))/1000000
#       return(total)
#     }
#   }
total_mem <- function() {
  tryCatch({
    if (Sys.info()[["sysname"]] == "Linux") {
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^MemFree:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Cached:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Inactive:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Buffers:' /proc/meminfo", intern = TRUE))))/1000
      return(total)
    } else if (Sys.info()[["sysname"]] == "Windows") {
      #total <- as.numeric(memory.limit())
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("wmic OS get FreePhysicalMemory /Value", intern=TRUE))[3]))/1000
      return(total)
    } else if (Sys.info()[["sysname"]] == "Darwin") { # Mac OS X
      total <- 4096*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages free'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages inactive'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages speculative'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages purgeable'", intern = TRUE))))/1000000
      return(total)
    }
    else { # other OS, i.e. Solaris
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'free memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'inactive memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'buffer memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'swap cache'", intern = TRUE))))/1000
      return(total)
    }
  }, error=function(e){
    return(1024)
  }, warning=function(e){
    return(1024)
  })
}
  ###############################################
  print("-----pepare to parallel------")
  parallel_threshold <- 4   # what is this ?
  if (mem.efficient) {
    mem_per_test <- 2 #MB
    tests_per_batch <- as.integer(total_mem() / mem_per_test)
  }
  time3 <- Sys.time()
  start_total <- proc.time()
  while (!done && any(G)){
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G)
    #     ## Consider only unique edge
    #     ind <- subset(ind, ind[1] < ind[2])
    remainingEdgeTests <- length(ind)
    if(verbose>=1)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    G.l <- split(G, gl(p+1,1))
    if (!mem.efficient) {
      tests_per_batch <- remainingEdgeTests
    }
    for (j in seq(1, remainingEdgeTests, by=tests_per_batch)){
      l <- min(remainingEdgeTests, j + tests_per_batch - 1)
      #cat("l - j + 1=",l - j + 1,"\n")
      #cat("parallel_threshold=",parallel_threshold,"\n")
      if (l - j + 1 < num_workers) {
        num_workers <- l - j + 1        
      }
      res <- NULL
      if (l - j + 1 < parallel_threshold) 
      {
        cat("l-j+1 = :",l-j+1,"parallel_threshold is :",parallel_threshold)
        print("*** system ***")
        res <- lapply(j:l, edge_test)
      } 
      else if (Sys.info()[['sysname']] == 'Windows') 
      {
        print("windows")
        res <- do.call("c", clusterApply(workers, clusterSplit(workers, j:l), edge_tests))
        #cat('res=',res,"\n")                                                # 
      } 
      else 
      {
        print("not windows")
        res <- mclapply(j:l, edge_test, mc.cores=num_workers, mc.set.seed=FALSE, 
                        mc.cleanup=TRUE, mc.allow.recursive=FALSE)
      }
      # synchronize
      for (p_obj in res) {
        #print("p_obj = ", p_obj)
        #cat("p_obj = ",p_obj)
        i <- p_obj[[1]]
        #x <- ind[i, 1]
        #y <- ind[i, 2]
        y <- 1
        x <- ind[i]
        n.edgetests[ord+1] <- n.edgetests[ord+1] + p_obj[[3]]
        #cat('p_obj[[3]]=',p_obj[[3]],"\n")
        #cat('n.edgetests[ord+1]=',n.edgetests[ord+1],"\n")
        #pMax[x, y] <- p_obj[[6]]
        zMin[x] <- p_obj[[4]]
        #zMin[y] <- p_obj[[5]]
        #G[x, y] <- G[y, x] <- p_obj[[2]]
        G[x] <-p_obj[[2]]
        done <- done & p_obj[[5]]
      }
    }
    # increase the nbrs size
    ord <- ord + 1
    
  }## end while
  
  #time4 <- Sys.time()
  print("parallel calculate cost:")
  print(Sys.time()-time3)
  total_t = proc.time()-start_total
  
  # write results
  cat('Num CI Tests=', n.edgetests, ',Total CI Tests=', 
      sum(unlist(n.edgetests)), ',Total Time=', total_t[3], '\n', sep=" ")
  cat("total is :",total_t,"\n")
  cat("cost time 1:",time1-time0,"cost time 2:",time2-time1,"\n")
  #print(paste0("---information about run time:---"),time2-time1, time3-time2)
  Gres <- G[-1]
  names(Gres) <- vNms
  list(G = Gres, zMin = zMin[-1])
}## pcSelect_parallel


#' Estimate subgraph around a response variable using pcSelect
#' 
#' @description 
#' This is the stable version (order independent version) of the pcSelect function (pc-Simple algorithm) in the pcalg package.
#' @param y The target (response) variable.
#' @param dm Data matrix with rows are samples and columns are variables.
#' @param alpha Significance level of individual partial correlation tests.
#' @param corMethod "standard" or "Qn" for standard or robust correlation estimation
#' @param method Character string specifying method; the default, "stable" provides an Order-independent version. 
#' @param verbose Logical or in \{0,1,2\};
#' 
#'        FALSE, 0: No output,
#'        
#'        TRUE, 1: Little output,
#'        
#'        2: Detailed output.
#'        
#'        Note that such output makes the function very much slower.
#' @param directed Logical; should the output graph be directed?
#' @return 
#' G    A logical vector indicating which column of dm is associated with y.
#' 
#' zMin   The minimal z-values when testing partial correlations between y and each column of dm. The larger the number, the more consistent is the edge with the data.
#' @examples
#' ##########################################
#' ## Using pcSelect_stable
#' ##########################################
#' library(pcalg)
#' library(parallel)
#' p <- 10
#' set.seed(101)
#' myDAG <- randomDAG(p, prob = 0.2)
#' n <- 1000
#' d.mat <- rmvDAG(n, myDAG, errDist = "normal")
#' pcSelect_stable(d.mat[,10],d.mat[,-10], alpha=0.05)
#' @export
pcSelect_stable <- function(y,dm, alpha, corMethod = "standard",method = "stable", verbose = FALSE, directed=FALSE)
{
  
  
  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  cl <- match.call()
  
  zMin <- c(0,rep.int(Inf,p))
  C <- pcalg::mcor(cbind(y,dm), method = corMethod)
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  #G <- c(FALSE,rep.int(TRUE,p))
  G <- matrix(c(FALSE,rep.int(TRUE,p)),nrow = p+1,ncol = 1)
  #G <- as.matrix(G)
  #G<-as.list(G)
  #View(G)
  #View(G1)
  seq_p <- 1:(p+1)
  #seq_p <- seq_len(p)
  
  done <- FALSE
  ord <- 0
  
  start_total <- proc.time()
  
  while (!done && any(G)) {
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G)
    remainingEdgeTests <- length(ind)
    if(verbose>=1)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    ##hushu
    if(method == "stable") {
      #View(G)
      #cat(typeof(G))
      ## Order-independent version: Compute the adjacency sets for any vertex
      ## Then don't update when edges are deleted
      G.l <- split(G, gl(p+1,1))
      #View(G.l)
      #G<-as.factor(G)
      #G.l <- split(G,G$factor)
    }
    #hushu
    for (i in 1:remainingEdgeTests) {
      if(verbose && i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n")
      y <- 1
      x <- ind[i]
      
      if (G[x]) {
        nbrsBool <- if(method == "stable") as.logical(G.l) else G 
        #nbrsBool <- G
        #View(G.l[[x]])
        #View(G.l[x])
        #View(G)
        nbrsBool[x] <- FALSE
        #View(seq_p[nbrsBool])
        #nbrs <- seq_p[nbrsBool[[i]]]
        nbrs <- seq_p[nbrsBool]
        ## neighbors of y without itself and x
        length_nbrs <- length(nbrs)
        
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq(length = ord)
          
          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat {
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            #cat('n.edgetests[ord+1]=',n.edgetests[ord+1],"\n")
            z <- pcalg::zStat(x,y, nbrs[S], C,n)
            #cat('x=',x,'y=',y,"\n")
            if(abs(z)<zMin[x]) zMin[x] <- abs(z)
            #cat('zMin_x=',zMin[x],"\n")
            if (verbose >= 2)
              cat(paste("x:",vNms[x-1],"y:",(ytmp <- round((p+1)/2)),"S:"),
                  c(ytmp,vNms)[nbrs[S]],paste("z:",z,"\n"))
            if (abs(z) <= cutoff) {
              G[x] <- FALSE
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast){
                break
              }
              S <- nextSet$nextSet
            }
          } ## {repeat}
        }
      } ## end if( G )
    } ## end for(i ..)
    ord <- ord+1
  } ## end while
  
  total_t = proc.time()-start_total
  
  
  # write results
  cat('Num CI Tests=', n.edgetests, ',Total CI Tests=', sum(unlist(n.edgetests)), ',Total Time=', total_t[3], '\n', sep=" ")
  
  Gres <- G[-1]
  names(Gres) <- vNms
  list(G = Gres, zMin = zMin[-1])
}## pcSelect_stable

