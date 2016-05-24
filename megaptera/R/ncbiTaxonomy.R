# RETRIEVE NCBI TAXONOMIES FOR CLADES
# PACKAGE: megaptera
# CALLED BY: USER
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-10-26

ncbiTaxonomy <- function(taxon, species.list = FALSE, kingdom, trim, parallel = FALSE, type = "SOCK"){
  
  kingdom <- match.arg(kingdom, c("Fungi", "Metazoa", "Viridiplantae"))
  
  if ( species.list ){
    spec <- taxon
    taxon <- unique(strip.spec(taxon))
  }
  
  id <- seq(from = 1, to= length(taxon), by = 50)
  id <- data.frame(from = id, to = c(id[-1] - 1, length(taxon)))
  x <- list()
  for ( i in 1:nrow(id) ){
    j <- seq(from = id$from[i], to = id$to[i])
    x <- c(x, ncbiLineage(taxon[j], species.list = species.list, 
                          parallel = parallel, kingdom = kingdom))
  }

  if ( length(x) == 0 ) return(NULL)
  
  ## delete ranks which are incertae sedis
  dis <- function(x){
    id <- grep("incertae sedis", x$name)
    if ( length(id) > 0 ){
      x <- x[-id, ]
    }
    x
  }
  x <- lapply(x, dis)
  
  ## check ranks and add columns if necessary
  ## ----------------------------------------
  ranks <- lapply(x, function(x) x$rank)
  rankSet <- unique(ranks)
  if ( length(rankSet) == 1 ){
    cat("\n.. unique rank set ..")
    ranks <- rankSet[[1]]
  } else {
    ranks <- sortRanks(rankSet)
    addRanks <- function(x, ranks){
      nr <- x$name[x$rank == "no rank"]
      id <- match(ranks, x$rank, incomparables = "no rank")
      xx <- x[id, ]
      xx$rank <- ranks
      
      nr <- x$name[x$rank == "no rank"]
      ff <- function(i, x){
        ff <- x$name[(which(x$name == i) + 1):nrow(x)]
        1:(min(which(x$name %in% ff)) - 1)
      }
      nrr <- lapply(nr, ff, x = x)
      names(nrr) <- nr
      nrl <- sapply(nrr, length)
      nid <- sort(nrl, FALSE)
      xx$name[nid] <- names(nid)
      
      ## fill remaining with "-"
      xx$name[is.na(xx$name)] <- "-"
      xx
    }
    x <- lapply(x, addRanks, ranks)
  } # end of IF clause
  
  x <- lapply(x, function(x) x$name)
  x <- do.call(rbind, x)
  colnames(x) <- ranks
  x <- data.frame(x)
  
  
  if ( species.list ){
    colnames(x)[colnames(x) == "species"] <- "genus"
    spec <- data.frame(strip.spec(spec), spec)
    x <- cbind(x[match(spec[, 1], x$genus), ], species = spec[, 2])
    
    no.info <- is.na(x$genus)
    if ( any(no.info) ){
      no.info <- which(no.info)
      warning("no classification found for: ", paste(x$species[no.info], collapse = ", "))
      x <- x[-no.info, ]
     
    }
    
  }
  cat("\n.. number of species found:", nrow(x), "..")
  
  ## delete internal, noninformative "no rank"-columns
  ## -------------------------------------------------
  id <- apply(x, 2, unique)
  id <- which(id == "-")
  if ( length(id) > 0 ){
    x <- x[, -id]
  }
  
  ## trimming taxonomy to required lower ranks
  ## -----------------------------------------
  if ( !missing(trim) ){
    trim <- match.arg(trim, c("auto", colnames(x)))
    if ( trim == "auto" ){
      id <- max(which(sapply(apply(x, 2, unique), length) == 1))
    }
    if ( trim %in% colnames(x) ){
      id <- which(colnames(x) %in% trim)
    }
    x <- x[, id:ncol(x)]
  } 
  
  ## order taxonomy
  ## --------------
  for ( i in rev(colnames(x)) ){
    x <- x[order(x[, i]), ]
  }
  x
}