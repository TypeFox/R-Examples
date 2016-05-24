


nodiv_data <- function(phylo, commatrix, coords, proj4string_in = CRS(as.character(NA)), type = c("auto", "grid", "points"), shape = NULL)
{
  type = match.arg(type)
  
  if(!(class(phylo) == "phylo")) stop ("phylo must be a phylogeny in the ape format") 
  
  if(!inherits(commatrix, "distrib_data"))
  {
    if(missing(coords)) 
      dist_dat <- distrib_data(commatrix, proj4string_in = proj4string_in, type = type, shape = shape) else
        dist_dat <- distrib_data(commatrix, coords, proj4string_in, type, shape)
  } else dist_dat <- commatrix
  
  if(is.null(dist_dat$species_stats))
    stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")
  
  nodiv_dat <- dist_dat[c("comm", "type", "coords", "species_stats", "hcom")]
  
  # TODO It should also be possible to give all of the below as function arguments, and also to the distrib_data functoin and integrate in subsample
  if(!is.null(dist_dat$shape))
    nodiv_dat$shape <- dist_dat$shape
  if(!is.null(dist_dat$node_stats))
    nodiv_dat$node_stats <- dist_dat$nodestats
  
  cat("Comparing taxon names in phylogeny and communities (using picante)\n")
  
  colnames(dist_dat$comm) <- match_speciesnames(phylo$tip.label, colnames(dist_dat$comm), do_not_match = TRUE)
  nodiv_dat$species_stats$species <- colnames(dist_dat$comm)
  temp <- capture.output(dat <- match.phylo.comm(phylo, dist_dat$comm))
  nodiv_dat$phylo <- dat$phy
  nodiv_dat$comm <- dat$comm
  if(!(is.data.frame(nodiv_dat$comm) & nrow(nodiv_dat$comm) > 1)) stop("The tip labels in the phylogeny do not match the names in the community matrix")
  phylodropped <- Ntip(phylo) - Ntip(nodiv_dat$phylo)
  commdropped <- ncol(dist_dat$comm) - ncol(nodiv_dat$comm)
  if(phylodropped > 0)
    cat(paste("  - removed ", phylodropped, " species from phylo that are not found in commatrix\n"))
  if(commdropped > 0)
    cat(paste("  - removed ", commdropped, " species from commatrix not found in phylo\n"))
  
  md <- match(colnames(nodiv_dat$comm), nodiv_dat$species_stats$species)
  nodiv_dat$species_stats <- subrow_data.frame(nodiv_dat$species_stats, md)
    
  nodiv_dat$coords <- dist_dat$coords[na.omit(match(rownames(nodiv_dat$comm), dist_dat$coords$sites)),]
  nodiv_dat$hcom <- matrix2sample(nodiv_dat$comm) # do I actually need this for anything?
  nodiv_dat$hcom[,1] <- as.character(nodiv_dat$hcom[,1])
  nodiv_dat$hcom[,3] <- as.character(nodiv_dat$hcom[,3])
  
  cat("Calculating which species descend from each node\n")
  nodiv_dat$node_species <- Create_node_by_species_matrix(nodiv_dat$phylo)
    
  class(nodiv_dat) <- c("nodiv_data","distrib_data")
  return(nodiv_dat)
}

is01line <- function(vec){
  num <- unique(as.numeric(vec))
  return(!sum(!num %in% 0:1))
}

BenHoltMatrix <- function(commatrix){
  last <- sapply(max(1, ncol(commatrix) - 3):ncol(commatrix), function(line) is01line(commatrix[, line]))
  if(!sum(last) == length(last))
    return(FALSE)
  first01 <- sapply(1:min(ncol(commatrix), 30), function(line) is01line(commatrix[, line]))  
  first01 <- min(which(first01))
  return(first01)
}

distrib_data <- function(commatrix, coords = NULL, proj4string_in = CRS(as.character(NA)), type = c("auto", "grid", "points"), shape = NULL)
{
  type = match.arg(type)
  cat("Checking input data\n")
  if(inherits(commatrix, "distrib_data")){
    if(is.null(commatrix$species_stats))
      stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")
    
    ret <- list()
    ret$comm <- commatrix$comm
    ret$coords <- commatrix$coords
    ret$type <- commatrix$type
    ret$shape <- commatrix$shape
    ret$species_stats <- commatrix$species_stats
    class(ret) <- "distrib_data"
    return(ret)
  }
  if(is.null(coords)){
    if(isWorldmapData(commatrix)){
      cat("Data format identified as Worldmap export file\n")
      coords <- data.frame(site = paste(commatrix[, 4], commatrix[, 5], sep = '_'), Long = commatrix[, 4], Lat = commatrix[, 5])
      commatrix <- data.frame(site = coords$site, abu = rep(1, nrow(commatrix)), species = commatrix[, 1])
      coords <- coords[!duplicated(coords$site), ]
      if(identical(proj4string_in, CRS(as.character(NA))))
        proj4string_in <- CRS("+proj=longlat +ellps=WGS84")
    } else {
      if((firstnumeric <- BenHoltMatrix(commatrix)) > 1){
        # this first bit should be moved up under isworldmapmatrix and be used to check for the Ben type of matrix perhaps
        cat(paste("Commatrix assumed to be a concatenation of coordinates (", firstnumeric - 1," columns) and community matrix\n", sep = ""))
        coords <- commatrix[, 1:(firstnumeric-1)]
        commatrix <- commatrix[, -(1:firstnumeric-1)]
      } else stop("If not commatrix is already of type distrib_data or nodiv_data, a worldmap matrix, or a concatenation of coords and community matrix, coords must be specified")       
    }     
  }
    
  

  ## Testing that input objects are all right
  if(class(coords) == "SpatialPointsDataFrame" | class(coords) == "SpatialPixelsDataFrame")
    if(!all.equal(proj4string_in,coords@proj4string))
    { 
      proj4string_in <- proj4string(coords)
      warning("specified proj4string overridden by the coords data")
    } 

  if(is.data.frame(commatrix) & ncol(commatrix) == 3 & !is.numeric(commatrix[,3])) #i.e. is the commatrix in phylocom format?
  {
    cat("Commatrix identified as phylocom format\n")
    commatrix[,1] <- as.character(commatrix[,1])
    commatrix[,3] <- as.character(commatrix[,3])
    commatrix <- sample2matrix(commatrix)   
  }

  if(is.data.frame(commatrix)) commatrix <- as.matrix(commatrix)
  if(!is.matrix(commatrix)) stop("commatrix must be a matrix of 0's and 1's, indicating presence or absence")
  if(!is.numeric(commatrix)) stop("commatrix must be a numeric matrix of 0's and 1's, indicating presence or absence")
  if(!sum(!unique(as.numeric(commatrix)) %in% 0:1) == 0) stop("commatrix must be a matrix of 0's and 1's, indicating presence or absence")
  
  temp <- floor(commatrix) #this is currently not necessary, due to the previous line
  if(sum(commatrix-temp) > 0)
    stop("commatrix had non-integer entries, please revise")
  
  if(is.matrix(coords)) coords <- as.data.frame(coords)
  cat("Transforming coords to spatial points\n")
  if(is.data.frame(coords)) coords <- toSpatialPoints(coords,proj4string_in, commatrix, type)

  if(class(coords) == "SpatialPixelsDataFrame") type <- "grid" else if (class(coords) == "SpatialPointsDataFrame") type <- "points" else stop("coords must be a data.frame of coordinates or an sp data.frame object")
  

  ## making sure that the points and the commatrix fit
  
  
  commatrix <- match_commat_coords(commatrix, coords$sites)  
  
  not.occurring.species <- which(colSums(commatrix) == 0)
  not.occupied.sites <- which(rowSums(commatrix) == 0)
  
  if(length(not.occupied.sites) > 0)
  {
    message(paste(length(not.occupied.sites), "sites where dropped because no species occupied them:\n", paste(coords$sites[not.occupied.sites], collapse = "\t")))
    coords <- coords[ - not.occupied.sites, ]
    commatrix <- commatrix[ - not.occupied.sites,  ]
  }
  
  if(length(not.occurring.species) > 0)
  {
    message(paste(length(not.occurring.species), "species where dropped because of 0 occurrences in the areas defined by coords:\n", paste(colnames(commatrix)[not.occurring.species], collapse = "\t")))
    commatrix <- commatrix[, - not.occurring.species]
  }
  
  ret <- list(comm = as.data.frame(commatrix), type = type, coords = coords)
  if(type == "grid")
    ret$grid <- summary(ret$coords)$grid 
  
  ret$species_stats <- data.frame(species = colnames(ret$comm), stringsAsFactors = FALSE)
  
  if(!is.null(shape)) ret$shape <- shape
  
  class(ret) <- "distrib_data"
  return(ret)

}


#internal functions

#TODO
#much of this testing can be done with try-catch phrases
#use the testthat library to test everything


match_commat_coords <- function(commatrix, sitenames)
{
  cat("Comparing sites in community data and spatial points\n")
  
  if(is.null(rownames(commatrix)))
    if(nrow(commatrix) == length(sitenames)) 
      rownames(commatrix) <- sitenames else stop("The number of sites in coords and the data matrix do not fit and there are no rownames in the community matrix to use for matching")  
      
  if(sum(sitenames %in% rownames(commatrix)) < 2)
    if(nrow(commatrix) == length(sitenames)) 
      rownames(commatrix) <- sitenames else
        stop("the coordinate names and the rownammes of the community matrix do not match")
  
  if(sum(rownames(commatrix) %in% sitenames) < length(sitenames)*0.8)
    if(nrow(commatrix) == length(sitenames)) {
      rownames(commatrix) <- sitenames
      cat("Rownames of the matrix ignored because not a sufficient match to sitenames\n")
    } else stop("The number of sites in coords and the data matrix do not fit and the sitenames did not match the matrix names sufficiently well for matching")  
  
  if(sum(sitenames %in% rownames(commatrix)) < length(rownames(commatrix)))
    cat(paste(length(rownames(commatrix)) - sum(sitenames %in% rownames(commatrix)), " sites removed from the dataset because the rownames of the data matrix did not match the site names. Make sure the rownames are correct\n"))
  
  sitenames <- sitenames[sitenames %in% rownames(commatrix)]

    
  commatrix <- commatrix[match(sitenames, rownames(commatrix)),]
  return(commatrix) 
}


toSpatialPoints <- function(coords, proj4string, commatrix, type)
{
  
    xcol <- 0
    ycol <- 0
    
    ret <- coords
    
    colnames(ret) <- tolower(colnames(ret))
    if('x' %in% colnames(ret) & 'y' %in% colnames(ret))
    {
      xcol <- which(colnames(ret) == 'x')
      ycol <- which(colnames(ret) == 'y')
      
    } else if('lon' %in% substr(colnames(ret),1,3) & 'lat' %in% substr(colnames(ret), 1, 3)) {
      colnames(ret) <- substr(colnames(ret), 1, 3)
      xcol <- which(colnames(ret) == 'lon')[1]
      ycol <- which(colnames(ret) == 'lat')[1]
    }

    names(coords)[xcol] = "myX"
    names(coords)[ycol] = "myY"
    
    cat("Identifying sites identifier\n")
    
    if (ncol(coords)==3 & !(xcol + ycol == 0) & isTRUE(all.equal(coords[,-c(xcol, ycol)], unique(coords[,-c(xcol, ycol)])))) names(coords)[!names(coords) %in% c("myX", "myY")] <- "sites" else 
      if(nrow(coords) == nrow(commatrix) & !is.null(rownames(commatrix)) & ncol(coords) == 2){
        if(is.null(rownames(coords)) | identical(rownames(coords), as.character(1:nrow(coords))))
          coords = data.frame(sites = rownames(commatrix)) else 
            if (!identical(rownames(commatrix), rownames(coords))) stop("Because the rownames of commatrix and coords differ, sitenames cannot be established unless they are included explicitly as a third column of coords")
      }  else {
        if(is.null(rownames(commatrix))){
          stop("There must be valid site names in the rownames of commatrix or in the coords data")
        } else {
          coords <- infer_sites_intern(rownames(commatrix), coords)
        }
      }
    
    if(sum(names(coords) == "sites") > 1)
      stop(paste("Could not match on the variable called sites, please rename"))
    
    coords$sites <- as.character(coords$sites)
    ids <- c(which(names(coords) == "myX"), which(names(coords) == "myY") )
    xy <- coords[, ids]    
    if(!ncol(xy) == 2) stop("coords should be a data.frame or spatial data.frame with 2 columns, giving the x/longitude, and y/latitude of all sites")
    
    xy <- SpatialPoints(xy, proj4string)
    type_auto <- ifelse(isGrid(xy), "grid", "points")
    
    if(type == "auto") type <- type_auto else 
      if(!type == type_auto)
        warning(paste("The specified type of data (", type, ") seems to conflict with the automatic setting. This may cause problems", sep = ""))
    
    if(ncol(coords) == 3)
      ret <- data.frame(sites = coords[, - ids], stringsAsFactors = FALSE) else
        ret <- data.frame(coords[, -ids], stringsAsFactors = FALSE)
    if(type == "grid") ret <- SpatialPixelsDataFrame(xy, ret) else
      ret <- SpatialPointsDataFrame(xy, ret)
    
    return(ret)  
}

isGrid <- function(coords)
  return(isGridVar(coordinates(coords)[,1]) & isGridVar(coordinates(coords)[,2]))

isGridVar <- function(gridVar)
{
  dists <- diff(sort(unique(gridVar)))
  distab <- table(dists)
  smallest <- as.numeric(names(distab[1]))
  most_common <- as.numeric(names(distab))[distab == max(distab)]
  return(isTRUE(all.equal(dists/smallest, floor(dists/smallest))) & smallest %in% most_common)
  #if all differences are a multiplum of the smallest, and the smallest distance is the most common, it is probably a grid
}


Create_node_by_species_matrix = function(tree)
{
  # create a matrix with 0s and 1s indicating which species descend from each node
  nodespecies <- matrix(0, nrow = Nnode(tree), ncol = Ntip(tree))
  colnames(nodespecies) <- tree$tip.label
  rownames(nodespecies) <- nodenumbers(tree)
 
  ntip <- Ntip(tree)
  .local <- function(tree, node)
  {
    if(node <= Ntip(tree))
      return(node)
    ret <- lapply(Descendants(node, tree), .local, tree = tree)
    ret <- do.call(c, ret)
    nodespecies[node - ntip, ret] <<- 1
    ret
  }
  .local(tree, basal_node(tree))
  
  return(nodespecies)
}

Node_spec <- function(tree, node, names = TRUE)
{
  .local <- function(tree, node)
  {
    if(node <= Ntip(tree))
      return(node)
    ret <- lapply(Descendants(node, tree), .local, tree = tree)
    do.call(c, ret)
  }
  
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")
  
  node <- identify_node(node, tree)
  ret <- .local(tree, node)
  
  if(names)
    ret <- tree$tip.label[ret]
  
  ret
}

isWorldmapData <- function(dat){
  if(is.data.frame(dat)){
    if(is.factor(dat[, 1]))
      dat[, 1] <- as.character(dat[, 1])
    if(ncol(dat) == 5)
      if(is.character(dat[, 1]))
        if(is.numeric(dat[, 4]))
          if(is.numeric(dat[, 5]))
            if(min(dat[, 4], na.rm = T) > -181 & max(dat[, 4], na.rm = T) < 181)
              if(min(dat[, 5], na.rm = T) > -91 & max(dat[, 5], na.rm = T) < 91)
                return(TRUE)
  }
  return(FALSE)
}

