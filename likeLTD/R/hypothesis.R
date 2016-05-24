# Reads crime scene profile from file.
# Documentation is in man directory.
read.csp.profile = function(path) {
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  fileType = file_ext(path)
  if((fileType=="xlsx")|(fileType=="xls"))
  	{
  	profile = read.xls(path, header=T, colClasses='character', sep=',')
  	} else {
  	profile = read.table(path, header=T, colClasses='character', sep=',')
  	}
  result = sapply(profile, function (n) strsplit(n, "(\\s*,|\\s)\\s*"))
  result = result[result[, "X"] != "Uncertain", 5:ncol(result), drop=FALSE]
  allNA = apply(result, 2, function(n) !all(is.na(n)))
  result[, allNA, drop=FALSE]
}

# Reads uncertain alleles in crime scene profile from file.
# Documentation is in man directory.
read.unc.profile = function(path) {
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  fileType = file_ext(path)
  if((fileType=="xlsx")|(fileType=="xls"))
  	{
  	profile = read.xls(path, header=T, colClasses='character', sep=',')
  	} else {
  	profile = read.table(path, header=T, colClasses='character', sep=',')
  	}
  result = sapply(profile, function (n) strsplit(n, "(\\s*,|\\s)\\s*"))
  result = result[result[, "X"] == "Uncertain", 5:ncol(result), drop=FALSE]
  allNA = apply(result, 2, function(n) !all(is.na(n)))
  result[, allNA, drop=FALSE]
}

# Reads known profiles from file.
# Documentation is in man directory.
read.known.profiles = function(path) {
  if(!file.exists(path)) stop(paste(path, "does not exist."))
  fileType = file_ext(path)
  if((fileType=="xlsx")|(fileType=="xls"))
  	{
  	profiles = read.xls(path, header=T, colClasses='character', row.names=1,sep=',', quote = "\"")
  	} else {
 	profiles = read.table(path, header=T, colClasses='character', row.names=1,sep=',', quote = "\"") 
  	}
  result = sapply(profiles, function (n) strsplit(n, "(\\s*,|\\s)\\s*"))
  # undoing R damage
  if(nrow(profiles) == 1) {
    result <- matrix(result, nrow=1)
    colnames(result) <- colnames(profiles)
  }
  row.names(result) = row.names(profiles)
  colnames(result)[1] = "queried"
  result[, "queried"] = result[, "queried"] == "queried"
  result
}

# Reads linkage information from file.
# Documentation is in man directory
load.linkage.info <- function(path=NULL) {
	# Loads allele database
	# Documentation is in man directory.
	if(is.null(path)) { # Load default database
	dummyEnv <- new.env()
	data('linkage', package='likeLTD', envir=dummyEnv)
	return(dummyEnv[['linkage']])
	}
	if(!file.exists(path)) stop(paste(path, "does not exist."))
	read.table(path, sep="\t", header=TRUE,as.is=TRUE)
	}

# Should check that hypothesis is syntactically correct
# This means matrices as opposed to lists, so on and so forth.
# It does not mean it the input should make any sort of sense.
sanity.check = function(hypothesis) {
  errors = c()
  if(!is.matrix(hypothesis$cspProfile)) 
     errors = c(errors, "cspProfile should be a matrix 'replicates' x 'loci'.")
  if(!is.matrix(hypothesis$dropoutProfs))
     errors = c( errors, 
                 "dropoutProfs should be a matrix 'profiles' x 'loci'.")
  if(!is.matrix(hypothesis$uncProf))
     errors = c(errors, "uncProf should be a matrix 'replicates' x 'loci'.")
  if(!is.matrix(hypothesis$queriedProfile))
     errors = c(errors, "queriedProfile should be a matrix 'profiles' x 'loci'.")
  if(length(errors) == 0) {
    if(nrow(hypothesis$cspProfile) != nrow(hypothesis$uncProf))
      errors = c( errors,
                  "Number of rows of uncProf and cspProfile do not match." )
    if(ncol(hypothesis$cspProfile) != ncol(hypothesis$uncProf))
      errors = c( errors,
                  "Number of loci of uncProf and cspProfile do not match." )
    if(ncol(hypothesis$cspProfile) != ncol(hypothesis$dropoutProfs))
      errors = c( errors,
                  "Number of loci of uncProf and dropoutProfs do not match." )
    if(ncol(hypothesis$cspProfile) != ncol(hypothesis$dropoutProfs))
      errors = c( errors,
                  "Number of loci of uncProf and dropoutProfs do not match." )
    if(ncol(hypothesis$cspProfile) != ncol(hypothesis$queriedProfile))
      errors = c( errors,
                  "Number of loci of uncProf and queriedProfs do not match." )
  }
  if(length(errors) != 0) {
    cat("There seems to be an error with the input hypothesis.\n")
    for(error in errors) cat(paste(c(error, "\n")))
    stop()
  }
}

# Adds dropout column to known profiles.
# Documentation is in man directory.
determine.dropout = function(knownProfiles, cspProfile) {
  # Alleles present in any replicate of the csp per locus.
  if(nrow(knownProfiles) == 0) return(vector(mode="logical"))
  per.indiv = function(indiv) {
    per.locus = function(locus) {
      all(sapply(cspProfile[, locus], function(n) indiv[[locus]] %in% n))
    }
    all(sapply(colnames(cspProfile), per.locus))
  }
  !apply(knownProfiles, 1, per.indiv)
}

masking.and.uncertain.profile = function(uncProfile, maskingProfiles) {
  # Adds masking profiles to each replicate of uncertain alleles.
  # 
  # Masking and uncertain alleles are somehow treated the same.
  if(length(maskingProfiles) == 0) return(uncProfile)
  for(locus in intersect(colnames(uncProfile), colnames(maskingProfiles))) {
    alleles = Reduce(union, maskingProfiles[, locus, drop=FALSE])
    for(i in 1:nrow(uncProfile))
      uncProfile[[i, locus]] = union(uncProfile[[i, locus]], alleles)
  }
  uncProfile
}

masking.profile = function(cspProfile, maskingProfiles) {
  # Remove masking alleles from CSP.
  if(length(maskingProfiles) == 0) return(cspProfile)
  for(locus in intersect(colnames(cspProfile), colnames(maskingProfiles))) {
    alleles = Reduce(union, maskingProfiles[, locus, drop=FALSE])
    for(i in 1:nrow(cspProfile))
      cspProfile[[i, locus]] = setdiff(cspProfile[[i, locus]], alleles)
  }
  cspProfile
}

missing.alleles = function(alleleDb, cspProfile, queriedProfile, dropoutProfiles) {
  # Adds missing alleles to database
  #
  # There may be alleles in the crime-scene profile, in the queried individual,
  # or in the individuals subject to dropout, which are not present in the
  # database. Theses alleles are added into the database with count 1 and
  # fragment length 0.
  if(!is.matrix(cspProfile)) stop("input should be a matrix.")
  if(is.null(colnames(cspProfile)))
    stop("input matrix does not have column names.")
  if(!is.matrix(dropoutProfiles)) stop("input should be a matrix.")
  if(!is.matrix(queriedProfile)) stop("input should be a matrix.")
  for(locus in colnames(cspProfile)) {
    dbAlleles = rownames(alleleDb[[locus]])
    cspAlleles = unique(unlist(cspProfile[, locus]))
    qAlleles = unique(unlist(queriedProfile[, locus]))
    dropAlleles = unique(unlist(dropoutProfiles[, locus]))
    allAlleles = unique(c(cspAlleles,qAlleles, dropAlleles))
    missingAlleles = allAlleles[!allAlleles%in%dbAlleles]
    # At this point, some of the missingAlleles might just be saying move
    # along, nothing to see.
    missingAlleles = setdiff(missingAlleles, c("", NA, "NA"))
    # Now add new rows to database. 
    if(length(missingAlleles)) {
      newrows = matrix(c(1, 0), nrow=length(missingAlleles), ncol=2,
                       byrow=TRUE)
      rownames(newrows) = missingAlleles
      alleleDb[[locus]] = rbind(alleleDb[[locus]], newrows)
    }
  }
  alleleDb
}

transform.to.locus.centric = function(hypothesis) {
  # Transform hypothesis to locus centric modes.
  # 
  # This means we reorganise the data to be in lists of the loci.
  result = list()
  for(locus in colnames(hypothesis$cspProfile)) {
    # Value of the resulting list for a given locus 
    locusValue = list()
    # Loop over all items in original list.
    for(key in names(hypothesis)) {
      value = hypothesis[[key]]
      # If a matrix and locus is either in rows or columns, then add only locus
      # part of the matrix. 
      # If a list, then add only the locus specific part of the list.
      if(is.matrix(value)) {
        if(locus %in% colnames(value)) 
          locusValue[[key]] = value[, locus, drop=FALSE]
        else if(locus %in% rownames(value))
          locusValue[[key]] = value[locus, , drop=FALSE]
      } else if(is.list(value) && (locus %in% names(value))) 
        locusValue[[key]] = value[[locus]]
      # If has not been added yet, then not locus dependent, and add it as a
      # whole. 
      if(!key %in% names(locusValue)) locusValue[[key]] = value 
    }
    # Only locus to result if not empty.
    if(length(locusValue) > 0) result[[locus]] = locusValue
  }
  result
}

# combine rare alleles into one joined allele for a single locus
combine.rares.locus = function(alleleDb,cspProfile,uncProfile,knownProfiles,queriedProfile,rareThreshold=0.05)
    {
    # index of alleles not in csp, unc, knowns or queried, and also probability < rareThreshold
    index = which(!rownames(alleleDb)%in%c(unlist(uncProfile),unlist(cspProfile),unlist(knownProfiles),unlist(queriedProfile))&alleleDb[,1]<rareThreshold)
    if(length(index)>0)
        {
        # remove indexed alleles, new allele has sum of probabilities, and mean of BP
        alleleDb = rbind(alleleDb[-index,],c(sum(alleleDb[index,1]),mean(alleleDb[index,2]))) 
        rownames(alleleDb)[nrow(alleleDb)] = "rares"
        }
    return(alleleDb)   
    }

# combine rare alleles into one joined allele for all loci
combine.rares = function(alleleDb, cspProfile, uncProfile, knownProfiles, queriedProfile,rareThreshold=0.05)
    {
    loci = colnames(cspProfile)
    sapply(loci,FUN=function(x) combine.rares.locus(alleleDb[[x]],cspProfile[,x],uncProfile[,x],knownProfiles[,x],queriedProfile[,x],rareThreshold=rareThreshold),simplify=FALSE)
    }


agnostic.hypothesis <- function(cspProfile, uncProfile, knownProfiles,
                                queriedProfile, alleleDb, ethnic='EA1',
                                adj=1e0, fst=0.02, combineRare=FALSE, 
                                rareThreshold=0.05) {
  # Helper function to figure out the input of most hypothesis.
  #
  # Basically, this groups operations that are done the same by defence and
  # prosection. 

 if(!ethnic%in%colnames(alleleDb)) stop("Chosen race code not included in database")

  # Read database and filter it down to requisite ethnicity and locus. 
  alleleDb = ethnic.database(ethnic, colnames(cspProfile), alleleDb)
 
  # Figure out which profiles show dropout.
  dropout = determine.dropout(knownProfiles, cspProfile)

  # Creates profiles with and without dropout.
  dropoutProfs = knownProfiles[dropout, colnames(cspProfile), drop=FALSE]
  noDropoutProfs = knownProfiles[!dropout, colnames(cspProfile), drop=FALSE]

  # Adjust uncertain profile to also contain masking (no-dropout) profiles.
  uncProf = masking.and.uncertain.profile(uncProfile, noDropoutProfs)

  # Remove masked alleles from  CSP.
  cspProf = masking.profile(cspProfile, noDropoutProfs)

  # Adjust database to contain all requisite alleles
  alleleDb = missing.alleles(alleleDb, cspProfile, queriedProfile, dropoutProfs)
  alleleDb = adjust.frequencies( alleleDb, 
                                 queriedProfile[1, colnames(cspProfile), 
                                                drop=FALSE],
                                 adj=adj, fst=fst )
  # combine rare alleles into a single allele
  if(combineRare) alleleDb = combine.rares(alleleDb, cspProf, uncProf, knownProfiles, queriedProfile, rareThreshold)

    
  # Construct all profiles as arrays of  
  list(cspProfile=cspProf, uncProf=uncProf, dropoutProfs=dropoutProfs,
       alleleDb=alleleDb,
       queriedProfile=queriedProfile[1, colnames(cspProfile), drop=FALSE]) 
}

# Creates prosecution hypothesis.
# Documentation is in man directory.
prosecution.hypothesis <- function(cspFile, refFile, ethnic='NDU1',
                                   nUnknowns=0, adj=1e0, fst=0.02,
                                   databaseFile=NULL, linkageFile=NULL,
                                   doDropin=FALSE, combineRare=TRUE,
				   rareThreshold=0.05, kit=NULL, relationship=0,...) {
  if(relationship>7|relationship<0) stop("Relationship must be specified between 0 and 7")
  if(is.null(databaseFile)&is.null(kit)) kit = "DNA17"
  relatedness = convertRelationship(relationship)
  alleleDb = load.allele.database(databaseFile,kit)
  cspProfile = read.csp.profile(cspFile)
  uncProfile = read.unc.profile(cspFile)
  if(!relationship%in%c(0,1))
	{
  	linkageInfo = load.linkage.info(linkageFile)
	rownames(linkageInfo) = linkageInfo[,1]
	linkageInfo = linkageInfo[,-1]
	linkIndex = which(colnames(linkageInfo)%in%colnames(cspProfile))
	linkageInfo = linkageInfo[linkIndex,linkIndex]
	}
  knownProfiles = read.known.profiles(refFile)
  if(sum(unlist(knownProfiles[, "queried"])) != 1)
    stop("Expect one queried profile on input.")
  qIndices = which(unlist(knownProfiles[, "queried"]))
  uIndices = which(!unlist(knownProfiles[, "queried"]))
  queriedProfile = knownProfiles[qIndices, , drop=FALSE]  
  # Puts queried profile at the end.
  knownProfiles = knownProfiles[c(uIndices, qIndices), , drop=FALSE] 

  result = agnostic.hypothesis(cspProfile, uncProfile, knownProfiles,
                               queriedProfile, alleleDb, ethnic=ethnic,
                               adj=adj, fst=fst, combineRare=combineRare,
				rareThreshold=rareThreshold)

  # If queried profile is subject to dropout, then it should be the reference
  # individual. Otherwise, the first individual subject to dropout will be set
  # as the reference individual.
  result[["refIndiv"]] = 1
  if(any(determine.dropout(queriedProfile, cspProfile)))
    {
  	# refIndiv is in relation to dropoutProfs rather than knownProfiles
  	nameRef = rownames(result$queriedProfile)
  	result[["refIndiv"]] = which(rownames(result$dropoutProfs)==nameRef)
	  }

  result = append(result, list(...))
  result[["nUnknowns"]] = nUnknowns
  result[["hypothesis"]] = "prosecution"
  result[["ethnic"]] = ethnic  
  result[["adj"]] = adj
  result[["fst"]] = fst
  result[["relatedness"]] = relatedness
  result[["relationship"]] = relationship 
  result[["doDropin"]] = doDropin
  result[["cspFile"]] = cspFile
  result[["refFile"]] = refFile
  result[["databaseFile"]] = databaseFile
  result[["kit"]] = kit

  if(!relationship%in%c(0,1))
	{
    result[["linkageInfo"]] = linkageInfo
    }

  sanity.check(result) # makes sure hypothesis has right type.
  result
}

# Creates defence hypothesis
# Documentation is in man directory.
defence.hypothesis <- function(cspFile, refFile, ethnic='NDU1',  nUnknowns=0,
                               adj=1e0, fst=0.02, databaseFile=NULL, linkageFile=NULL, 
                               doDropin=FALSE, combineRare=TRUE, 
			       rareThreshold=0.05, kit=NULL, relationship=0,...) {
  if(relationship>7|relationship<0) stop("Relationship must be specified between 0 and 7")
  if(is.null(databaseFile)&is.null(kit)) kit = "DNA17"
  relatedness = convertRelationship(relationship)
  alleleDb = load.allele.database(databaseFile,kit)
  cspProfile = read.csp.profile(cspFile)
  uncProfile = read.unc.profile(cspFile)
  if(!relationship%in%c(0,1))
	{
  	linkageInfo = load.linkage.info(linkageFile)
	rownames(linkageInfo) = linkageInfo[,1]
	linkageInfo = linkageInfo[,-1]
	linkIndex = which(colnames(linkageInfo)%in%colnames(cspProfile))
	linkageInfo = linkageInfo[linkIndex,linkIndex]
	}
  knownProfiles = read.known.profiles(refFile)
  if(sum(unlist(knownProfiles[, "queried"])) != 1)
    stop("Expect one queried profile on input.")
  queriedProfile = knownProfiles[unlist(knownProfiles[, "queried"]), ,
                                 drop=FALSE]
  # Removing queried individual from knownProfiles. 
  knownProfiles = knownProfiles[!unlist(knownProfiles[, "queried"]), ,
                                drop=FALSE]

  result = agnostic.hypothesis(cspProfile, uncProfile, knownProfiles,
                               queriedProfile, alleleDb, ethnic=ethnic,
                               adj=adj, fst=fst, combineRare=combineRare,
				rareThreshold=rareThreshold) 


  result[["refIndiv"]] = nrow(result[["dropoutProfs"]]) + 1

  result = append(result, list(...))
  result[["nUnknowns"]] = nUnknowns + 1
  result[["hypothesis"]] = "defence"
  result[["ethnic"]] = ethnic 
  result[["adj"]] = adj
  result[["fst"]] = fst
  result[["relatedness"]] = relatedness
  result[["relationship"]] = relationship 
  result[["doDropin"]] = doDropin
  result[["cspFile"]] = cspFile
  result[["refFile"]] = refFile
  result[["databaseFile"]] = databaseFile
  result[["kit"]] = kit
  if(!relationship%in%c(0,1)) 
	{
	result[["linkageInfo"]] = linkageInfo
	}
  sanity.check(result) # makes sure hypothesis has right type.
  result
}


add.args.to.hypothesis <- function(hypothesis, ...) {
  # Updates hypothesis with input arguments.
  #
  # Convenience function to easily modify hypothesis when creating the objective
  # functions. 
  #
  # Parameters:
  #   hypothesis: The hypothesis to modify.
  #   ...: Named parameters to add/update the hypothesis
  # Returns: The modified hypothesis  

  args = list(...)
  argnames = intersect(names(hypothesis), names(args))
  if(length(argnames)) hypothesis[argnames] = args[argnames]
  argnames = setdiff(names(args), names(hypothesis))
  if(length(argnames)) hypothesis = append(hypothesis, args[argnames])

  hypothesis
}
