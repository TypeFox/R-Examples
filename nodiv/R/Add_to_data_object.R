# I made some major changes here and to infer_sites, that also need to be done for species!
add_sitestat <- function(distrib_data, site_stat, site = NULL){
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")

  if(is.matrix(site_stat))
    site_stat <- as.data.frame(site_stat, stringsAsFactors = FALSE)

  if(is.vector(site_stat)){
    nam <- deparse(substitute(site_stat))
    if(is.null(names(site_stat)) & length(site_stat) == Nsites(distrib_data)){
      if(nam %in% names(distrib_data$coords@data))
        warning(paste("Overwriting the contents of", nam))
      distrib_data$coords@data[[nam]] <- site_stat
      return(distrib_data)
    }
    rownam <- names(site_stat)
    site_stat <- as.data.frame(site_stat, stringsAsFactors = FALSE)
    names(site_stat) <- nam
    rownames(site_stat) <- rownam
  }

  if(is.null(site))
  {
    temp <- infer_sites(distrib_data, site_stat)
    site <- temp$site
    site_stat <- temp$site_stat
  }

  site <- identify_sites(site, distrib_data)

 # if(length(site) < Nsites(distrib_data))
  #  cat(paste(Nsites(distrib_data)- length(site), "sites from site_stat were not found in", deparse(substitute(distrib_data)), "\n"))

  mergeframe <- as.data.frame(lapply(site_stat, function(column) {
    ret <- vector(mode = typeof(column), length = Nsites(distrib_data))
    ret[] <- NA
    if(is.factor(column))
      ret <- factor(ret, levels = levels(column))
    ret[site] <- column
    ret
  }), stringsAsFactors = FALSE)

  if(sum(names(mergeframe) %in% names(distrib_data$coords@data)) > 0){
    matches <- which(names(distrib_data$coords@data) %in% names(mergeframe))
    deleted <- names(distrib_data$coords@data)[matches]
    distrib_data$coords <- distrib_data$coords[, -matches]
    warning(paste("Some data in the original distrib_data overwritten:\n"), paste(deleted, collapse = "\t"))
  }


  distrib_data$coords@data <- cbind(distrib_data$coords@data, mergeframe)
  distrib_data
}

add_species_stat <- function(distrib_data, species_stat, specs = NULL){
  #restorepoint::restore.point("add_species_stat")


  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  if(is.null(distrib_data$species_stats))
    stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")

  if(is.matrix(species_stat))
    species_stat <- as.data.frame(species_stat, stringsAsFactors = FALSE)

  if(is.vector(species_stat)){
    nam <- deparse(substitute(species_stat))
    species_stat <- as.data.frame(species_stat, stringsAsFactors = FALSE)
    names(species_stat) <- nam
  }

  num <- nrow(species_stat)

  if(is.null(specs))
  {
    temp <- infer_species(distrib_data, species_stat)
    specs <- temp$species
    species_stat <- temp$species_stat
  }

  specs <- identify_species(specs, distrib_data)

  if(length(specs) < Nspecies(distrib_data))
    cat(paste(num - length(specs), "species were not found in", deparse(substitute(distrib_data)), "\n"))

  mergeframe <- as.data.frame(lapply(species_stat, function(column) {
    ret <- vector(mode = typeof(column), length = Nspecies(distrib_data))
    ret[] <- NA
    if(is.factor(column)) ret <- factor(ret, levels = levels(column))
    ret[specs] <- column
    ret
  }), stringsAsFactors = FALSE)

  if(sum(names(mergeframe) %in% names(distrib_data$species_stats)) > 0){
    matches <- which(names(distrib_data$species_stats) %in% names(mergeframe))
    deleted <- names(distrib_data$species_stats)[matches]
    distrib_data$species_stats[,matches] <- NULL
    warning(paste("Some data in the original distrib_data overwritten:\n"), paste(deleted, sep = "\t"))
  }


  distrib_data$species_stats <- cbind(distrib_data$species_stats, mergeframe)
  distrib_data
}



infer_sites_intern <- function(sites, site_stat) # a non-exported convenience function
{
#   suppressWarnings(numsites <- as.numeric(sites)) #I removed this as it caused trouble
#   if(sum(is.na(numsites)) < 0.2*length(sites)){
#     if(all.equal(numsites, floor(numsites))) {      # if site names are just integers, matching is not attempted
#       if(nrow(site_stat) == length(sites))
#         return(list(site = sites, site_stat = site_stat)) else
#           warning("Site matching was done based on name matching, which is tricky when site names are integer values")
#     }
#   }
#
  continue <- FALSE
  temp <- 0
  
  potnams <- c("sites", "site","Sites", "Site", "plot", "Plot", "cell", "Cell", "ID", "id", "Centroid", "centroid") 
  
  potid <- which(names(site_stat) %in% potnams)
  
  
  if(is.null(rownames(site_stat)))
    continue <- TRUE else {
    if(length(potid) == 0){
      potentials <- as.list(as.data.frame(rownames(site_stat)))
      names(potentials) <- "rownames"
    } else {
      potentials <- as.list(as.data.frame(site_stat[[potid]]))
      names(potentials) <- names(site_stat)[potid]    
      potentials[["rownames"]] <- rownames(site_stat)
    }

    temp <- sapply(1:length(potentials), function(index){
      matches <- sum(potentials[[index]] %in% sites)
      return(matches/nrow(site_stat))
    })
  
    res <- which(temp == max(temp))
    if(length(res) > 1){
      res <- res[1]
      warning(paste(length(res), "variables had an equally good correspondence to the sitenames:", max(temp), ". Using the first of these,", names(potentials)[res], "to align"))
    }
    name <- names(potentials)[res]
    site <- potentials[[res]]
    
    if(temp[res] < 0.8) continue <- TRUE
  }
    
  if (continue) {
    temp2 <- sapply(1:length(site_stat), function(index){
      matches <- sum(site_stat[[index]] %in% sites)
      return(matches/nrow(site_stat))
    })
    
    if(max(temp2) > max(temp)){
      temp <- temp2
      res <- which(temp == max(temp))[1]
      name <- names(site_stat)[res]
      site <- site_stat[[res]]      
    }
  }


  ##### We need a matching function here to do the actual matching!

  site_stat <- subrow_data.frame(site_stat, which(!is.na(site)))
  site <- site[!is.na(site)]

  site_stat_ret <- subrow_data.frame(site_stat,match(sites, as.character(site)))

  cat(paste("Matching sites by", name, "\n"))
  hits = sum(site %in% sites)
  if(hits == length(site) & hits == length(site))
    cat("All sites matched\n") else
    cat(paste(hits, " sites were matched:\n\t", floor(hits/length(site) * 100), "% of ", length(site), " sites in site_stat\n\t", floor(hits/length(sites) * 100), "% of ", length(sites), " sites in distrib_data\n", sep = ""))

  if(name == "rownames"){
    site_stat_ret$sites <- as.character(site)[match(sites, as.character(site))]
    name <- "sites"
  }


  index <- which(sapply(site_stat_ret, function(x) identical(x, site_stat_ret[[name]])))[1]

  names(site_stat_ret)[index] <- "sites"

  if(sum(names(site_stat_ret) == "sites") > 1)
    stop(paste("Could not match on the variable called sites, as", name, "had a greater correspondence. Please rename"))

  site_stat_ret$sites <- as.character(site_stat_ret$sites)

  return(site_stat_ret)
}


infer_sites <- function(distrib_data, site_stat) # a non-exported convenience function
{
  ret <- infer_sites_intern(sites(distrib_data), site_stat)

  suppressWarnings(sitenames <- identify_sites(ret$sites, distrib_data, as.name = TRUE))
  matchsite <- match(ret$sites, sitenames)

  site_stat <- subrow_data.frame(ret, which(!is.na(matchsite)))
  site <- site_stat$sites
  site_stat$sites <- NULL

  return(list(site = site, site_stat = site_stat))
}


match_speciesnames <- function(reference_name, new_name, do_not_match = FALSE){
  chars <- c(" ", "[.]", "_")
  ref_ll <- sapply(chars, function(x) length(grep(x, reference_name)))
  ref_char <- chars[which(ref_ll == max(ref_ll))[1]]
  new_ll <- sapply(chars, function(x) length(grep(x, new_name)))
  new_char <- chars[which(new_ll == max(new_ll))[1]]
  new_name <- gsub(new_char, ref_char, new_name)
  if(do_not_match)
    return(new_name)

  ret <- match(new_name, reference_name)
  ret
}

infer_species <- function(distrib_data, species_stat) # a non-exported convenience function
{
    species_stat$rownames <- rownames(species_stat)
    temp <- sapply(1:length(species_stat), function(index){
      matches <- match_speciesnames(species(distrib_data), species_stat[[index]])
      return(sum(!is.na(matches))/nrow(species_stat))
    })
    res <- which(temp == max(temp))[1]
    name <- names(species_stat)[res]
    spec <- species(distrib_data)[match_speciesnames(species(distrib_data), species_stat[[res]])]

  if(temp[res] < 0.5 & sum(spec %in% species(distrib_data)) < 0.5 * Nspecies(distrib_data))
    stop("Species could not be matched automatically, please supply the species argument explicitly")

  if(!name == "rownames")
    species_stat[[name]] <- NULL

  species_stat$rownames <- NULL


  species_stat <- subrow_data.frame(species_stat, which(!is.na(spec)))
  spec <- spec[!is.na(spec)]


  suppressWarnings(specnames <- identify_species(spec, distrib_data, as.name = TRUE))
  matchspec <- match(spec, specnames)

  species_stat <- subrow_data.frame(species_stat, which(!is.na(matchspec)))
  spec <- spec[!is.na(matchspec)]

  cat(paste("Matching species by", name, "\n"))

  return(list(species = spec, species_stat = species_stat))
}

