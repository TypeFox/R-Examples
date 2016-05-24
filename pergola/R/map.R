
#' Creates map object
#' 
#' Creates map object from matrix of pairwise recombination frequencies.
#'   
#' @param rf Matrix of pairwise recombination frequencies.
#' @param split Split object.
#' @param fun Function to space the markers on the map. 
#' Default is "haldane". Alternatives are "kosambi", "carter" and "none.
#' @param corr Corrector, if recombinations are overestimated.
#' Allows to multiply all spaces by a fixed factor.
#' @return Ordered vector of marker locations. Each marker has a name attribute.
#' @examples
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' pullMap(rfMat, split = split)   
#' @export
pullMap <- function(rf, split, fun = "haldane", corr = 1){
  fun <- match.arg(fun, c("haldane", "kosambi", "carter", "none"))
  rf <- rf[split$order, split$order]
  msplit <- max(split$split)
  globalMap <- vector(mode = "list", length = msplit)
  for(i in 1:msplit){
    subset <- which(split$split[split$order] == i)
    lss <- length(subset)
    map <- rep(0, lss)
    for(j in 2:lss){
      map[j] <- rf[subset, subset][j, j - 1]
    }
    map <- switch(fun, 
                  haldane = -log(1 - 2 * map) * 100,
                  kosambi = 0.25 * log((1 + 2 * map) / (1 - 2 * map)) * 100,
                  carter = 0.25 * (0.5 * (log(1 + 2 * map) - log(1 - 2 * map)) + atan(2 * map)) * 100,
                  none = map)
    map <- cumsum(map) * corr
    names(map) <- split$names[split$order][split$split[split$order] == i]
    if(any(is.infinite(map))) map <- map[!is.infinite(map)]
    globalMap[[i]] <- map
  }  
  return(globalMap)
}


#' Switch Chromosomes
#' 
#' Wrapper function to switch chromosomes for the whole map
#'  
#' @param map Map to switch.
#' @param comp Other map for comparison.
#' @return map
#' @examples 
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' map <- pullMap(rfMat, split = split)   
#' split <- sortLeafs(rfMat, split, method = "endlink")
#' map2 <- pullMap(rfMat, split = split)   
#' map <- switchChrs(map, map2)
#' 
#' @export
switchChrs <- function(map, comp){
  lapply(map, function(x) switchChr(x, comp))
}

#' Switch Chromosomes
#' 
#' Switches chromosomes if comparison map fits better.
#' This does improve the comparability without changing the order, 
#' as chromosomes are independent.
#'  
#' @param map Map to switch.
#' @param comp Other map for comparison.
#' @return map
#' @keywords internal
switchChr <- function(map, comp){
  markers <- names(map)
  matchchr <- findChr(markers, comp)
  nmark <- length(markers)
  cnm <- ceiling(nmark / 2)
  front <- markers[1:cnm]
  back <- markers[(cnm + 1):nmark]
  markers <- names(comp[[matchchr]])
  nmark <- length(markers)
  front2 <- markers[1:cnm]
  back2 <- markers[(cnm + 1):nmark]
  ff <- sum(front %in% front2)
  fb <- sum(front %in% back2)
  bf <- sum(back %in% front2)
  bb <- sum(back %in% back2)
  if(ff < fb && bb < bf){
    offs <- map[1] + map[length(map)]
    map <- offs - map
    return(rev(map))
  }
  map
}

#' Find most equal chromosome in other map
#' 
#' 
#'  
#' @param map Map to switch.
#' @param comp Other map for comparison.
#' @return map
#' @keywords internal
findChr <- function(map, comp){
  matches <- sapply(comp, function(x) sum(names(x) %in% map))
  which(matches == max(matches))[1]
}

#' Swap chromosomes
#' 
#' Finds best matching chromosome for each chromosome and brings them into the same order.
#'  
#' @param map Map to switch.
#' @param comp Other map for comparison.
#' @return map
#' @examples 
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' map <- pullMap(rfMat, split = split)   
#' split <- sortLeafs(rfMat, split, method = "endlink")
#' map2 <- pullMap(rfMat, split = split)   
#' map <- swapChrs(map, map2)
#' @export
swapChrs <- function(map, comp){
  chrs <- sapply(map, function(x) findChr(names(x), comp))
  if(any(duplicated(chrs))){
    warning("Not all chromosomes could be matched with their corresponding chromosome.")
  }else{
    map <- map[order(chrs)]
  }
  map
}



#' Transforming a map into a dendrogram
#' 
#' Create dendrogram object. The map specific distance are ignored and only the grouping and ordering is maintained.
#' Allows for comparison of whole map with package 'dendextend'
#'     
#' @param map One map. Required.
#' @param mergeoff Numeric, offset between chromosomes, to avoid equal heights in dendrogram.
#' Equal heights lead to problems in cor_bakers_gamma().
#' @return Dendogram object.
#' @import stats
#' @examples
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' map <- pullMap(rfMat, split = split)   
#' dend <- map2dend(map)  
#' plot(dend)
#' @export
map2dend<-function(map, mergeoff = 0L){
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("dendextend needed for this function to work. Please install it.",
         call. = FALSE)
  }
  nChr <- length(map)
  mergeval <- (1:nChr) * mergeoff + max(unlist(map)) * 1.2
  nMar<-sapply(map, length)
  out <- NULL
  for(i in 1:nChr){
    hclu <- hclust(dist(map[[i]]), method = "complete")
    dend <- as.dendrogram(hclu)
    dend <- dendextend::rotate(dend, names(map[[i]]))
    if(is.null(out)){
      out <- dend
    }else{
      out <- merge(out, dend, height = mergeval[i])
    }
  }
  return(out)  
}

#' Add offset
#' 
#' Add offset to zero distance markers to allow computation of correlation between maps.
#'     
#' @param map One map. Required.
#' @param offset Numeric value for offset.
#' 
#' @return Map object.
#' @examples
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' map <- pullMap(rfMat, split = split)   
#' map <- add_offset(map)
#' @export
add_offset<-function(map, offset = 0.1){
  for(i in 1:length(map)){
    while(any(duplicated(map[[i]]))){
      map[[i]][duplicated(map[[i]])]<-map[[i]][duplicated(map[[i]])] + offset
    }
  }
  map
}
