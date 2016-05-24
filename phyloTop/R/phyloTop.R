#' phyloTop: topological properties of phylogenies
#' 
#' Calculate a range of topological properties for one or more phylogenetic trees.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param treeList a \code{list} or \code{multiPhylo} object, or a single tree of class \code{phylo} or \code{phylo4}. All trees should be binary and rooted; if not they will be coerced into binary rooted trees using multi2di, if possible.
#' @param funcs a list of functions. The default is to apply all of the topological functions from the package, but a subset can be specified instead. The functions available are:
#' \itemize{
#' \item \code{\link{avgLadder}}
#' \item \code{\link{cherries}}
#' \item \code{\link{colless.phylo}}
#' \item \code{\link{ILnumber}}
#' \item \code{\link{maxHeight}}
#' \item \code{\link{pitchforks}}
#' \item \code{\link{sackin.phylo}}
#' \item \code{\link{stairs}} (note that this adds two columns to the output, "stairs1" and "stairs2")
#' }
#' @param normalise option to normalise the results of functions where possible. Default is \code{FALSE}

#' @return A matrix where rows correspond to trees and columns correspond to topological properties.
#' 
#' @import ape
#' 
#' @seealso \code{\link{avgLadder}}, \code{\link{cherries}}, \code{\link{colless.phylo}}, \code{\link{ILnumber}}, \code{\link{maxHeight}}, \code{\link{pitchforks}}, \code{\link{sackin.phylo}}, \code{\link{stairs}}
#'   
#' @examples
#' ## Apply all of the functions to a list of 10 random trees, each with 50 tips:
#' phyloTop(rmtree(10,50))
#' ## Normalising the results where possible:
#' phyloTop(rmtree(10,50), normalise=TRUE)
#' 
#' 
#' @export
phyloTop <- function(treeList,funcs="all", normalise=FALSE){
  
  # check input:
  if (!class(treeList) %in% c("list","multiPhylo")) stop("Please supply a list or multiPhylo object for treeList")

  # functions which return an integer tree statistic, and their dependencies:
  # avgLadder (ladderSizes)
  # cherries (nConfig)
  # colless.phylo (treeImb, nConfig)
  # ILnumber (treeImb, nConfig)
  # maxHeight (getDepths)
  # pitchforks (nConfig)
  # sackin.phylo (getDepths)
  # stairs(1 and 2) (treeImb, nConfig)
  
  allfuncs <- c("avgLadder","cherries","colless.phylo","ILnumber","maxHeight","pitchforks","sackin.phylo","stairs")
  
  if (funcs=="all") {funcs <- allfuncs}
  else {
    # check each func is recognised
    for (i in funcs) {
      if (!i %in% allfuncs) {
        stop(paste0("Function '",i,"' not recognised."))
      }
    }
    # put in alphabetical order, mainly because "stairs" needs to be at the end!
    funcs <- sort(funcs)
  }
  
  
  
  # initialise matrix
  # if "stairs" is requested, need two columns for it
  if ("stairs" %in% funcs) {treeSummaries <- matrix(nrow=length(treeList),ncol=(length(funcs)+1),0)
  colnames(treeSummaries) <- c(setdiff(funcs,"stairs"),"stairs1","stairs2")}
  
  else {treeSummaries <- matrix(nrow=length(treeList),ncol=length(funcs),0)
  colnames(treeSummaries) <- funcs }

  
  # go through each tree at a time
  for (i in 1:length(treeList)) {
    ntips <- length(treeList[[i]]$tip.label)
    
    # get the functions which are reused
    # can we be more efficient than this repeated calling of phyloCheck?
    if (any(is.element(c("maxHeight", "sackin.phylo"),funcs))) {depths <- getDepths(treeList[[i]]) }
    if (any(is.element(c("cherries", "colless.phylo","ILnumber","pitchforks","stairs"),funcs))) {
      
      nConfig <- nConfig(treeList[[i]])
      
      # and find treeImb without recalling nConfig:
      
      nn=treeList[[i]]$Nnode
      Ancs=(ntips+1):(ntips+nn) 
      
      # for each internal node, find its immediate children 
      Pointers=t(vapply(Ancs, function(x) treeList[[i]]$edge[treeList[[i]]$edge[,1]==x,2], FUN.VALUE=c(1,2))) 
      
      configs <- nConfig$cladeSizes
      imbalance <- t(sapply(c(1:ntips,Ancs), function(node) {
        if (node <= ntips) {return(c(0,0))}
        else {
          children <- Pointers[node-ntips,]
          left <- configs[[children[[1]]]]
          right <- configs[[children[[2]]]]
          return(c(left,right))
        }
      }))
      treeImb <- imbalance}
    
    
    # apply each function
    for (j in funcs) {
      if (j=="avgLadder") { l <- ladderSizes(treeList[[i]])$ladderSizes
                            if (length(l)==0) {treeSummaries[i,j] <- 0} 
                            else {
                              if (normalise==FALSE) {treeSummaries[i,j] <- mean(l) }
                              else {
                                if (ntips==2) {treeSummaries[i,j] <- 0 }
                                else {treeSummaries[i,j] <- mean(l)/(ntips-2) }
                              }
                            }
      }
      
      if (j=="cherries") { 
        if (normalise==FALSE) {treeSummaries[i,j] <- nConfig$numClades[[2]]}
        else {treeSummaries[i,j] <- 2*nConfig$numClades[[2]]/ntips}
      }
      
      else if (j=="colless.phylo") { if (ntips==2) {treeSummaries[i,j] <- 0 }
                                else { diffs <- abs(apply(treeImb,1,diff))
                                n <- ((ntips-1)*(ntips-2))/2
                                treeSummaries[i,j] <- sum(diffs)/n
                                }
        }
      else if (j=="ILnumber") {
        if (ntips==2) {treeSummaries[i,j] <- 0} # if N=2 the result is 0 (and we should not try to normalise it!)
        else {NDs <- treeImb[(ntips+1):(2*ntips-1),]
              if (normalise==FALSE) {treeSummaries[i,j] <- sum(apply(NDs,1, function(x) sum(x==1)==1))}
              else {treeSummaries[i,j] <- (sum(apply(NDs,1, function(x) sum(x==1)==1)))/(ntips-2)}
        }
      }
      
      else if (j=="maxHeight") {
        heights <- depths$tipDepths
        if (normalise==FALSE) {treeSummaries[i,j] <- max(heights)}   
        else {treeSummaries[i,j] <- max(heights)/(ntips - 1)}
      }
      
      else if (j=="pitchforks") { 
        if (ntips==2) { treeSummaries[i, j] <- 0 }
        else if (normalise==FALSE) {treeSummaries[i,j] <- nConfig$numClades[[3]]}
        else {treeSummaries[i,j] <- 3*nConfig$numClades[[3]]/ntips}
      }
      
      else if (j=="sackin.phylo") {tipDepths <- depths$tipDepths - 1
                                    if (normalise==FALSE) {
                                      treeSummaries[i,j] <- sum(tipDepths)}
                                    else {treeSummaries[i,j] <- (sum(tipDepths))/((1/2)*(ntips*(ntips+1)) - 1)}
                                   
      }
      
      else if (j=="stairs") {
        if(ntips==2){
          treeSummaries[i, "stairs1"] <- 0
          treeSummaries[i, "stairs2"] <- 1
        }
        else {NDs <- treeImb[(ntips + 1):(2 * ntips - 1),]
              stair1 <- (1/(ntips - 1)) * sum(abs(NDs[2, ] - NDs[1, ]))
              stair2 <- (1/(ntips - 1)) * sum(pmin(NDs[2, ], NDs[1, ])/pmax(NDs[2, ], NDs[1, ]))
              treeSummaries[i,"stairs1"] <- stair1
              treeSummaries[i,"stairs2"] <- stair2
        }
      }
    }
  }
  return(as.data.frame(treeSummaries))
}

