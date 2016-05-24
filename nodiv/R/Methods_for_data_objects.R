
## TODO only generate the node_species object in the Node_analysis function, and remove it from the summary functions (if feasible - check for calls)

head.distrib_data <- function(x, ...) print.distrib_data(x, ...)
head.nodiv_data <- function(x, ...) print.nodiv_data(x, ...)
head.nodiv_result <- function(x, ...) print.nodiv_result(x, ...)

coords <- function(distrib_data){
  if(!inherits(distrib_data, "distrib_data"))
    stop("object must be of class distrib_data")
  ret <- coordinates(distrib_data$coords)
  rownames(ret) <- sites(distrib_data)
  ret
}

occurrences <- function(distrib_data, species, value = c("index", "names", "logical", "raw")){
  value = match.arg(value)
  if (!inherits(distrib_data, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  specs <- identify_species(species, distrib_data)
  
  ret <- distrib_data$comm[,specs]
  if(value == "logical")
    ret <- (ret > 0)
  
  if(is.null(dim(ret))){  #if it is just a vector
    if(value == "index" | value == "names")
      ret <- which(ret > 0)
    if(value == "names")
      ret <- sites(distrib_data)[ret]    
  } else {
    if(value == "index" | value == "names"){
      ret <- lapply(1:ncol(ret), function(x) which(ret[,x] > 0))
      names(ret) <- species
    }
    if(value == "names")
      ret <- lapply(ret, function(x) sites(distrib_data)[x])
  }
 
  ret  
}


#This function presently throws out the site statistics - instead it should be aggregating them!
gridData <- function(dist_data, cellsize_x = 1, cellsize_y = cellsize_x, xll_corner, yll_corner){
  if (!inherits(dist_data, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  if(missing(xll_corner))
    xll_corner <- cellsize_x * floor(min(coordinates(dist_data$coords)[, 1], na.rm = TRUE) / cellsize_x)
  if(missing(yll_corner))
    yll_corner <- cellsize_y * floor(min(coordinates(dist_data$coords)[, 2], na.rm = TRUE) / cellsize_y)
  if(dist_data$type == "grid"){
    cat('dist_data is already of type grid - resampling instead\n')
    dist_data$grid <- summary(dist_data$coords)$grid
    diffs <- c(cellsize_x/dist_data$grid$cellsize[1], cellsize_y/dist_data$grid$cellsize[2]) 
    if(min(diffs < 1))
      stop(paste('cannot resample grid to a finer scale - original cellsize is', paste(dist_data$grid$cellsize, collapse = ','), 'target cellsize is', paste(c(cellsize_x), collapse = ',' )))
       
     if(!identical(diffs, floor(diffs)))
       stop(paste('ratio between cellsizes must be an integer - original cellsize is', paste(dist_data$grid$cellsize, collapse = ','), 'target cellsize is', paste(c(cellsize_x), collapse = ',' )))
  }

  newx <- cellsize_x * floor((coordinates(dist_data$coords)[, 1] - xll_corner) / cellsize_x) + 0.5 * cellsize_x + xll_corner
  newy <- cellsize_y * floor((coordinates(dist_data$coords)[, 2] - yll_corner) / cellsize_y) + 0.5 * cellsize_y + yll_corner
  newsites <- paste(newx, newy, sep = '_')
  newcoords <- data.frame(site = newsites, X = newx, Y = newy)
  newcoords <- newcoords[!duplicated(newcoords$site), ]
  newhcom <- matrix2sample(dist_data$comm)
  newhcom$plot <- newsites[match(newhcom$plot, dist_data$coords$sites)]
  newcomm <- as.data.frame.matrix(table(newhcom$plot, newhcom$id))

  if(max(dist_data$comm) == 1)
    newcomm[newcomm > 0] <- 1

  newdist_data <- distrib_data(commatrix = newcomm, coords = newcoords, proj4string_in = CRS(proj4string(dist_data$coords)), type = "grid")

  dist_data$comm <- newdist_data$comm
  dist_data$coords <- newdist_data$coords
  dist_data$type <- "grid"
  dist_data
}

assemblage <- function(distrib_data, site, value = c("index", "names", "logical", "raw")){
  value = match.arg(value)
  if (!inherits(distrib_data, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  site <- identify_sites(site, distrib_data)
  
  ret <- distrib_data$comm[site, ]
  if(value == "logical")
    ret <- (ret > 0)
  if(value == "index" | value == "names")
    ret <- which(ret > 0)
  if(value == "names")
    ret <- species(distrib_data)[ret]
  
  ret   
}

Nspecies <- function(distrib_data)
{
  if (!inherits(distrib_data, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  length(species(distrib_data))
}

Nsites<- function(distrib_data)
{
  if (!inherits(distrib_data, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  length(distrib_data$coords)
}

print.distrib_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type," distributions of", Nspecies(x),"species in", Nsites(x),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(species(x)[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(sites(x)[1:printlen], collapse = ", "),", ...\n", sep = ""))
}

print.nodiv_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type,"distributions and phylogenetic relationships of", Nspecies(x),"species in", Nsites(x),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(species(x)[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(sites(x)[1:printlen], collapse = ", "),", ...\n", sep = ""))
}

identify.distrib_data <- function(x, ...)
  identify(coordinates(x$coords),  ...)

summary.distrib_data <- function(object, ...)
{
  if(object$type == "grid")
    richness <- suppressWarnings(SpatialPixelsDataFrame(SpatialPoints(object$coords), data.frame(richness = richness(object)))) else
      richness <- SpatialPointsDataFrame(SpatialPoints(object$coords), data.frame(richness = richness(object))) #Is it really necessary to make this a spatial data frame? 
      
  occupancy <- occupancy(object) 
  ret <- list(species = species(object), coords = object$coords, richness = richness, occupancy = occupancy, type = object$type)
  class(ret) <- "summary_distrib_data"
  ret
}

print.summary_distrib_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type,"distributions of", length(x$species),"species in", length(x$coords$sites),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("Species richness:  min", min(x$richness$richness), "max", max(x$richness$richness), "mean",mean(x$richness$richness),"\n"))
  cat(paste("Species occupancy:  min", min(x$occupancy, "max", max(x$occupancy), "mean",mean(x$occupancy),"\n")))
}

summary.nodiv_data <- function(object, ...)
{
  ret <- summary.distrib_data(object, ...)
  ret$nodes <- nodenumbers(object$phylo)
  ret$node.label = object$phylo.node.label
  class(ret) <- "summary_nodiv_data"
  ret
}

print.summary_nodiv_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type,"distributions and phylogenetic relationships of", length(x$species) ,"species in", length(x$coords$sites),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("Species richness:   min:", min(x$richness$richness), "\tmax:", max(x$richness$richness), "\tmean:", round(mean(x$richness$richness),2),"\n"))
  cat(paste("Species occupancy:  min:", min(x$occupancy), "\tmax:", max(x$occupancy), "\tmean:",round(mean(x$occupancy),2),"\n\n"))
  cat(paste("The phylogeny has", length(x$nodes), "internal nodes"))
  if (!is.null(x$node.label)) 
  {
    cat("Node labels:\n")
    if (x$nodes > printlen) 
    {
      cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "),", ...\n", sep = ""))
    } else print(x$node.label)
  }
}

add_shape <- function(distrib_data, shape)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("argument must be an object of types distrib_data, nodiv_data or nodiv_results")
  distrib_data$shape <- shape
  distrib_data
}

plot.distrib_data <- function(x, ...)
{
  if(is.null(x$shape)) shape <- NULL else shape <- x$shape
  if(x$type == "grid")
    plot_grid(richness(x), x$coords, shape = shape, ...) else
    plot_points(richness(x), x$coords, shape = shape, ...)
}  

plot_richness <- function(distrib_data, ...)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("argument must be an object of type distrib_data, nodiv_data or nodiv_results")
  plot.distrib_data(distrib_data, ...)
}

plot_node <- function(nodiv_data, node = basal_node(nodiv_data), sites = NULL, ...)
{
  if(!inherits(nodiv_data, "nodiv_data"))
    stop("argument must be an object of type nodiv_data or nodiv_result")
  node <- identify_node(node, nodiv_data)
  plot_richness(subsample.distrib_data(nodiv_data, species = Node_species(nodiv_data, node), sites = sites), ...)
}

plot.nodiv_data <- function(x,  ...)
{
  oldpar <- par()
  par(mfrow = c(1,2))
  plot.distrib_data(x, ...)
  plot(x$phylo, show.tip.label = isTRUE(Nspecies(x) < 40), cex = 0.7) #need to specify explicitly which

  par(mfrow = c(1,1))
}

subsample<- function(x, ...) UseMethod("subsample")

subsample.distrib_data <- function(x, sites = NULL, species = NULL, ...)
{
  #restorepoint::restore.point("subsample.distrib_data", TRUE)
  if(is.null(x$species_stats))
    stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")  
  
  if(inherits(sites, "SpatialPoints")) sites <- infer_sites(x, sites@data)$site
  
  keep_sites <- F
  if(is.null(sites)) sites <- 1:Nsites(x)
  if(is.logical(sites)) sites <- which(sites)
  if(is.character(sites)) 
  {
    if(length(sites) == 1)
    {
      if(sites == "all") keep_sites <- T 
    } else 
      sites <- match(sites, x$coords$sites)
  }  
  if(keep_sites) sites <- 1:Nsites(x)
  
  keep_species <- F
  if(is.null(species)) species <- 1:Nspecies(x)
  if(is.logical(species)) species <- which(species)
  if(is.character(species))
  {
    if(length(species) == 1)
    {
      if(species == "all") keep_species <- T 
    } else
    species <- match(species, x$species_stats$species)
  }
  
  if(keep_species) species <- 1:Nspecies(x)
  
  ret <- x[c("comm", "type", "coords", "species_stats")]
  if(!is.null(x$shape)) 
    ret$shape <- x$shape

  
  ret$comm <- ret$comm[sites, species]
  
  if(keep_sites) sites_keep <- rep_len(TRUE, nrow(ret$comm)) else sites_keep <- which(rowSums(ret$comm, na.rm = T) > 0)
  
  if(keep_species) species_keep <- rep_len(TRUE, ncol(ret$comm)) else species_keep <- which(colSums(ret$comm, na.rm = T) > 0)

  ret$comm <- ret$comm[sites_keep, species_keep]
  
  ret$species_stats <- subrow_data.frame(ret$species_stats, species[species_keep])
  ret$coords <- ret$coords[ret$coords$sites %in% rownames(ret$comm),]
  
  
  class(ret) <- "distrib_data"
  return(ret)
}

subsample.nodiv_data <- function(x, sites = NULL, species = NULL, node = NULL, ...)
{
  if(is.null(x$species_stats))
    stop("The nodiv_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")
  
  ret_phylo <- x$phylo
  ret_phylo$node.label <- nodenumbers(x)  #this line
  
  if(!is.null(node))
  {
    node <- identify_node(node, x)
    ret_phylo <- extract.clade(ret_phylo, node)
    species <- ret_phylo$tip.label
  } 
  
  ret <- subsample.distrib_data(x, sites, species)
  temp <- capture.output(dat <- match.phylo.comm(ret_phylo, ret$comm))
  
  new_phylo <- dat$phy
  old_nodes <- as.numeric(new_phylo$node.label)
  ret$phylo <- drop.tip(x$phylo, which(! species(x) %in% new_phylo$tip.label))  #this line
  ret$species_stats <- subrow_data.frame(x$species_stats, match(ret$phylo$tip.label, x$species_stats$species))
  
  ret$hcom <- subset(x$hcom, x$hcom$plot %in% ret$coords$sites & x$hcom$id %in% ret$species_stats$species)
  ret$node_species <- x$node_species[, colnames(x$node_species) %in% ret$species_stats$species]
  ret$node_species <- ret$node_species[rowSums(ret$node_species) > 0,]
  ret$node_species <- ret$node_species[as.numeric(rownames(ret$node_species)) %in% old_nodes,] # and this line are an ugly hack to make sure the node_species matrix does not get perverted
  if(!is.matrix(ret$node_species)) ret$node_species <- rbind(ret$node_species) #TODO this is a hack for when a node only has tips
  rownames(ret$node_species) <- nodenumbers(ret$phylo)
  class(ret) <- c("nodiv_data", "distrib_data")
  attr(ret, "old_nodes") <- old_nodes
  return(ret)
}

sites <- function(distrib_data){
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  return(distrib_data$coords@data$sites)
}


species <- function(distrib_data){
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  if(is.null(distrib_data$species_stats))
    stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")
  return(distrib_data$species_stats$species)
}

richness <- function(distrib_data, sites = NULL)
{  
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  
  if(!is.null(sites)){
    sites <- identify_sites(sites, distrib_data)
  } else sites <- 1:Nsites(distrib_data)

  if(length(sites) == 1)
    return(sum(distrib_data$comm[sites, ] > 0, na.rm = T))
  
  return(rowSums(distrib_data$comm[sites, ] > 0, na.rm = T))
}

occupancy <- function(distrib_data, species = NULL)
{  
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  
  if(!is.null(species)){
    species <- identify_species(species, distrib_data)
  } else species <- 1:Nspecies(distrib_data)
  
  if(length(species) == 1)
    return(sum(distrib_data$comm[, species] > 0, na.rm = T))
  
  return(colSums(distrib_data$comm[, species] > 0, na.rm = T))
}

plot_sitestat <- function(distrib_data, x, shape = NULL, type = c("auto", "points","grid"), ...)
{
  coords = NULL
  type = match.arg(type)
  if(inherits(distrib_data, "distrib_data"))
  {
    coords = distrib_data$coords
    if(!is.null(distrib_data$shape)) {
      if(is.null(shape))
        shape <- distrib_data$shape else
        warning("overriding the shape file associated with distrib_data")
    }
    
    if(is.character(x))
      if(length(x) == 1)
        x <- sitestat(distrib_data, x) 
      
    if(type == "auto")
      type <- distrib_data$type
    if(!type == distrib_data$type)
      warning(paste("Argument type has value",type,"but distrib_data is of type", distrib_data$type,"; this can cause problems or crashes"))
    nsit = Nsites(distrib_data)
  } 
  
  if(inherits(distrib_data, "SpatialPoints")){
    coords <- distrib_data
    nsit <- nrow(coordinates(coords))
  }

  if(is.matrix(distrib_data))
    distrib_data <- data.frame(distrib_data)
  if(is.data.frame(distrib_data))
    if(ncol(distrib_data) == 2){
      coords <- distrib_data
      nsit <- nrow(coords)
    }

  if(is.null(coords))
    stop("Wrong argument type for distrib_data - should be a distrib_data objects, spatial points or a two-column matrix of x and y values")
  
  if(type == "auto") type <- ifelse(isGrid(coords), "grid", "points")

  if(!length(x) == nsit)
    stop(paste("x must be a numeric vector of length", nsit, "or the name of a site variable in distrib_data"))

  if(type == "grid")
    plot_grid(x, coords, shape = shape, ...) else
      plot_points(x, coords, shape = shape, ...)
}

sitestat <- function(distrib_data, statname = NULL, site = NULL)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  
  if(is.null(site)) 
    site <- 1:Nsites(distrib_data) else site <- identify_sites(site, distrib_data)
  
  sitestatnames <- names(distrib_data$coord@data)
  
  if(is.null(statname))
    return(sitestatnames)
  
  if(length(statname) == 1)
    if(statname == "all")
      statname <- sitestatnames
  
  fitnames <- which(statname %in% sitestatnames)
  
  if(length(fitnames) == 0)
    stop(paste("Sitestat", paste(statname, collapse = ", ") , "not found in distrib_data!\nPotential sitestats are", paste(sitestatnames, collapse = ", ")))
  if(length(fitnames) < length(statname))
    warning(paste("Dropping sitestats:", paste(statname[-fitnames], collapse = ", ") , "not found in distrib_data"))
  
  ret <- distrib_data$coord@data[,sitestatnames %in% statname]
  if(is.vector(ret) | is.factor(ret)) {
    names(ret) <- sites(distrib_data)
    ret <- ret[site]
  } else ret <- ret[site,]
  ret
}

update_object <- function(distrib_data){
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  
  if(is.null(distrib_data$species_stats))
    distrib_data$species_stats <- data.frame(species = colnames(distrib_data$comm), stringsAsFactors = FALSE)
  if(!is.null(distrib_data$species))
    distrib_data$species <- NULL
  
  distrib_data
}

species_stat <- function(distrib_data, statname = NULL, specs = NULL)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  if(is.null(distrib_data$species_stats))
    stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")
  
  if(is.null(specs))
    specs <- 1:Nspecies(distrib_data) else specs <- identify_species(specs, distrib_data)
  
  if(length(specs) == 0 | (length(specs) == 1 & is.na(specs[1])))
    stop("Species not found in the dataset")
  
  species_statnames <- names(distrib_data$species_stats)
  
  if(is.null(statname))
    return(species_statnames)
  
  fitnames <- which(statname %in% species_statnames)
  
  if(length(fitnames) == 0)
    stop(paste("Species stat \"", paste(statname, collapse = ", ") , " \" not found in distrib_data!\nPotential species stats are:", paste(species_statnames, collapse = ", ")))
  if(length(fitnames) < length(statname))
    warning(paste("Dropping species_stats", paste(statname, collapse = ", ") , "not found in distrib_data"))
  
  ret <- distrib_data$species_stats[,species_statnames %in% statname]
  if(is.vector(ret)){
    names(ret) <- species(distrib_data)
    ret <- ret[specs]
  } else ret <- ret[specs,]
ret
}

plot_species <- function(distrib_data, species, col = c("darkgreen", "red"), ...)
{
  species <- identify_species(species, distrib_data)
  if (length(species) > 1)
  {
    warning("species had length > 1; only the first species is plotted")
  }
  if(is.null(list(...)$main)) main = distrib_data$species_stats$species[species]
  plot_sitestat(distrib_data, distrib_data$comm[,species], col = col, main = main, ...)
}

identify_species <- function(species, distrib_data, as.name = FALSE)
{
  if(inherits(distrib_data, "distrib_data")){
    specieslist <- species(distrib_data)
    speciesnumber <- Nspecies(distrib_data)
  } else if(inherits(distrib_data, "phylo")){
    specieslist <- distrib_data$tip.label
    speciesnumber <- Ntip(distrib_data)
  } else stop("distrib_data must be an object of type ape:::phylo, distrib_data, nodiv_data or nodiv_result")
  
  if(!is.vector(species)) 
    stop("species must be either numeric or character")
  
  if(is.character(species))
  {
    specs <- match(species, specieslist)
    omits <- which(is.na(specs))
    specs <- specs[!is.na(specs)]
    if(length(specs) == 0){
      warning("Species not found in the dataset") 
      return(NA)
    }

    if(length(omits) > 0){
      warning(paste("These species were excluded as they were not found:", paste(species[omits], collapse = ", ")))
    }
  } else specs <- species


  diagnost <- which(specs > speciesnumber | specs < 0)
  if(length(diagnost) == length(specs))
    return(NA)
  if(length(diagnost) > 0){
    warning(paste("numbers", paste(diagnost, sep = ", "), "are too high and did not match the community matrix"))
    specs <- specs[-diagnost]
  }

  if(as.name)
    specs <- specieslist[specs]
  
  specs
}

identify_sites <- function(site, distrib_data, as.name = FALSE)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type  distrib_data, nodiv_data or nodiv_result")
  
  if(is.character(site))
  {
    sitenames <- site
    site <- match(site, sites(distrib_data))
    omits <- which(is.na(site))
    site <- site[!is.na(site)]
    if(length(site) == 0){
      warning("Site names did not match") 
      return(NA)
    }
    
    if(length(omits) > 0){
      warning(paste("These sites were excluded as they were not found:", paste(sitenames[omits], collapse = "\t")))
    }
  } else site <- site
  
  
  diagnost <- which(site > Nsites(distrib_data) | site < 0)
  if(length(diagnost) == length(site))
    return(NA)
  if(length(diagnost) > 0){
    warning(paste("numbers", paste(diagnost, sep = ", "), "are too high and did not match the site number"))
    site <- site[-diagnost]
  }
  
  if(as.name)
    site <- sites(distrib_data)[site]
  
  site
}

