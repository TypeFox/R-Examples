# function to subset data
subsetData = function(data,search)
	{
	# subset columns by search parameter
	dataout = data[,grep(search,colnames(data)),drop=FALSE]
	# add first column (for replicate names etc)
	# can add this back in if necessary
	#dataout = cbind(data[,1],dataout)
	#colnames(dataout)[1] = "info"
	# set rownames to marker names
	rownames(dataout) = data$Marker
	# remove columns that have no data
	index = which(apply(dataout,MARGIN=2,FUN=function(x) all(is.na(x))))
	if(length(index!=0)) dataout = dataout[,-index,drop=FALSE]
	# remove AMEL data
	index = which(rownames(dataout)%in%c("AMEL","amel","Amel"))
	if(length(index)!=0) dataout = dataout[-index,,drop=FALSE]
	return(dataout)
	}

# function to check that the data make sense
checkData1 = function(alleles,heights)
	{
	# check if dimensions are all the same
	check1 = identical(dim(alleles),dim(heights))
	# check if all have same loci
	check2 = all(mapply(FUN=function(a,b) a==b,rownames(alleles),rownames(heights)))
	return(all(check1,check2))
	}

# function to check data makes sense between replicates
checkData2 = function(alleleReps,heightReps=NULL)
	{
	# loci names for each replicate
	lociReps = sapply(alleleReps,FUN=function(x) rownames(x),simplify=FALSE)
	foo = function(a,b) all(a==b)
	check1 = all(sapply(lociReps,FUN=function(x) identical(lociReps[[1]],x)))
	return(all(check1))
	}




# function to read in peak height data
read.peaks.profile = function(FILE)
	{
	if(!file.exists(FILE)) stop(paste(FILE, "does not exist."))
	# get data
	data = read.csv(file=FILE,as.is=TRUE,na.strings=c("NA",""))
	# get replicate names
	conditions = names(table(data[,1]))
	# get allele data
	alleles = sapply(conditions, FUN=function(x) subsetData(data[which(data[,1]==x),],"Allele"),simplify=FALSE)
	# get height data
	heights = sapply(conditions, FUN=function(x) subsetData(data[which(data[,1]==x),],"Height"),simplify=FALSE)
	# checks
	check1 = !mapply(FUN=checkData1,alleles,heights)
	if(any(check1)) stop(paste("Error in the CSP provided - replicates ",paste(which(check1),collapse=",")))
	check2 = !checkData2(alleles,heights)
	if(check2) stop(paste("Error in the CSP provided - replicates do not match"))
	# return data
	return(list(alleles=alleles,heights=heights))
	}

convert.to.binary = function(data)
	{
	# get rid of spaces
	data = as.matrix(data)
	# if no alleles, return empty matrix
	if(ncol(data)==0) 
		{
		allelic = sapply(1:nrow(data),FUN=function(x) NULL,simplify=TRUE)
		names(allelic) = rownames(data)
		return(allelic)
		}
	# remove spaces
	data = matrix(gsub(" ","",data),ncol=ncol(data),dimnames=list(rownames(data),colnames(data)))
	# replace .0 with nothing
	data = gsub("[.]0","",data)
	# allelic calls
	allelic = sapply(1:nrow(data),FUN=function(x) data[x,which(!is.na(data[x,]))],simplify=FALSE)
	names(allelic) = rownames(data)
	return(allelic)
	}

# combine rare alleles into one joined allele for a single locus
combine.rares.locus.peaks = function(alleleDb,cspProfile,knownProfiles,queriedProfile,rareThreshold=0.05,doDoubleStutter=FALSE,combineStutters=TRUE)
    {
    # not in any profiles or in overstutter positions of profiles
    inProfile = !rownames(alleleDb)%in%c(unlist(cspProfile),unlist(knownProfiles),unlist(queriedProfile))
    isRare = alleleDb[,1]<rareThreshold
    if(!combineStutters)
	{
    	inOverStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))+1,as.numeric(unlist(knownProfiles))+1,as.numeric(unlist(queriedProfile))+1)
    	inUnderStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))-1,as.numeric(unlist(knownProfiles))-1,as.numeric(unlist(queriedProfile))-1)
    	}
    # index of alleles not in csp, unc, knowns or queried, and also probability < rareThreshold
    if(doDoubleStutter&!combineStutters)
        {
        inDoubleOverStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))+1,as.numeric(unlist(knownProfiles))+1,as.numeric(unlist(queriedProfile))+1)
        inDoubleUnderStutter = !as.numeric(rownames(alleleDb))%in%c(as.numeric(unlist(cspProfile))-1,as.numeric(unlist(knownProfiles))-1,as.numeric(unlist(queriedProfile))-1)
        index = which(inProfile&inOverStutter&inUnderStutter&isRare&inDoubleOverStutter&inDoubleUnderStutter)
        } else {
	if(!combineStutters)
		{
        	index = which(inProfile&inOverStutter&inUnderStutter&isRare)
		} else {
		index = which(inProfile&isRare)
        	} 
	}
    if(length(index)>0)
        {
        # remove indexed alleles, new allele has sum of probabilities, and mean of BP
        #alleleDb = rbind(alleleDb[-index,],c(sum(alleleDb[index,1]),mean(alleleDb[index,2]))) 
	if(all(is.na(alleleDb[index,3]))) meanLUS = mean(alleleDb[,3],na.rm=TRUE) else meanLUS = mean(alleleDb[index,3],na.rm=TRUE)
	alleleDb = rbind(alleleDb[-index,,drop=FALSE],c(sum(alleleDb[index,1]),mean(alleleDb[index,2]),meanLUS))
	#alleleDb = rbind(alleleDb[-index,],c(sum(alleleDb[index,1]),mean(alleleDb[index,2])))
        rownames(alleleDb)[nrow(alleleDb)] = "-1"
        }
    return(alleleDb)   
    }

# combine rare alleles into one joined allele for all loci
combine.rares.peaks = function(alleleDb, cspProfile, knownProfiles, queriedProfile,rareThreshold=0.05,doDoubleStutter=FALSE)
    {
    loci = colnames(cspProfile)
    sapply(loci,FUN=function(x) combine.rares.locus.peaks(alleleDb[[x]],cspProfile[,x],knownProfiles[,x],queriedProfile[,x],rareThreshold=rareThreshold,doDoubleStutter=doDoubleStutter),simplify=FALSE)
    }

add.stutter.index = function(alleleDb)
	{
 	# make i
  	index = which(!as.numeric(rownames(alleleDb))<0)
  	foo=function(x) 2/(max(x)-min(x))*(x-max(x))+1
  	outIndex = vector(length=nrow(alleleDb))
  	outIndex[index] = foo(as.numeric(rownames(alleleDb)[index]))
  	outIndex[-index] = 0
	return(cbind(alleleDb,outIndex))
	}

fill.unknown.LUS = function(alleleDb)
	{
 	# make i
  	index = which(is.na(alleleDb[,3]))
	toadd = sapply(index,FUN=function(x) extrapolateLUS(alleleDb,x))
	alleleDb[index,3] = toadd
	return(alleleDb)
	}

extrapolateLUS = function(alleleDb,toFill)
	{
	# get allele names
	alleleNames = as.numeric(rownames(alleleDb))
	# if a negative allele, give mean LUS
	if(alleleNames[toFill]<0) return(mean(alleleDb[,3],na.rm=TRUE))
	# get allele endings e.g. 9.3 ending is 0.3
	alleleEndings = alleleNames-floor(alleleNames)
	# if not partial add one (so not dividing by zero)
	if(alleleEndings[toFill]==0) alleleEndings = alleleEndings+1
	# get which db have same partial ending and arent NA
	sameEndingIndex = which(round(alleleEndings[toFill],1)%%round(alleleEndings,1)==0&!is.na(alleleDb[,3]))
	# find closest out of those
	closest = which.min(alleleNames[sameEndingIndex]-alleleNames[toFill])
	# if none have same partial ending, move to absolute closest e.g. closest to 9.3 will be 9
	if(length(closest)==0)
		{
		difference = abs(alleleNames-alleleNames[toFill])
		acceptableIndex = !is.na(alleleDb[,3])
		closestIndex = which.min(difference[acceptableIndex])
		closest = which(acceptableIndex)[closestIndex]
		out = alleleDb[closest,3]-(alleleNames[closest]-alleleNames[toFill])
		} else {
		out = alleleDb[sameEndingIndex[closest],3]-(alleleNames[sameEndingIndex[closest]]-alleleNames[toFill])
		}
	return(out)
	}

extrapolateBP = function(alleleDb, newAllele)
	{
	alleles = as.numeric(rownames(alleleDb))
	BPs = alleleDb[,2]
	newEnding = strsplit(newAllele,"[.]")[[1]][2]
	newAllele = as.numeric(newAllele)
	# if -10 give mean BP
	if(newAllele<0) 
		{
		return(mean(BPs[which(alleles>0)],na.rm=TRUE))
		}
	allEndings = sapply(rownames(alleleDb),FUN=function(x) strsplit(x,"[.]")[[1]][2])
	# get possible levels e.g. partial and whole repeats
	levels = names(table(allEndings,useNA="ifany"))
	# need two to find repeat length
	pairs = which(table(allEndings,useNA="ifany")>1)
	toCheck = levels[pairs[1]]
	if(is.na(toCheck))
		{
		matching = which(is.na(allEndings))
		} else {
		matching = which(allEndings==toCheck)
		}
	# get repeat length
	repeatLength = (BPs[matching[1]]-BPs[matching[2]])/(alleles[matching[1]]-alleles[matching[2]])
	# get alleles that have same ending as new allele
	if(is.na(newEnding))
		{
		matching = which(is.na(allEndings))
		} else {
		matching = which(allEndings==newEnding)
		}
	if(length(matching)>0)
		{
		# simple same ending
		newLength = BPs[matching[1]] + repeatLength*(newAllele-alleles[matching[1]])
		} else {
		# no matching ending
		if(any(is.na(allEndings)))
			{
			# try whole repeats first, simple
			naIndex = which(is.na(allEndings))
			newLength = BPs[naIndex[1]] + repeatLength*(floor(newAllele)-alleles[naIndex[1]])
			newLength = newLength + as.numeric(newEnding)
			} else {
			# only partials
			if(is.na(newEnding)) newEnding=0
			newLength = BPs[1] + repeatLength*(floor(newAllele)-floor(alleles[1]))
			newLength = newLength + (as.numeric(newEnding)-as.numeric(allEndings[1]))
			}
		}
	return(newLength)
	}

agnostic.hypothesis.peaks <- function(cspProfile, knownProfiles,
                                queriedProfile, alleleDb, ethnic='NDU1',
                                adj=1e0, fst=0.02, combineRare=FALSE, 
                                rareThreshold=0.05,doDoubleStutter=FALSE) {
  # Helper function to figure out the input of most hypothesis.
  #
  # Basically, this groups operations that are done the same by defence and
  # prosection. 

  if(!ethnic%in%colnames(alleleDb)) stop("Chosen race code not included in database")

  # Read database and filter it down to requisite ethnicity and locus. 
  alleleDb = ethnic.database.lus(ethnic, colnames(cspProfile), alleleDb)
  #alleleDb = likeLTD:::ethnic.database(ethnic, colnames(cspProfile), alleleDb)

  # Adjust database to contain all requisite alleles
  alleleDb = missing.alleles.peaks(alleleDb, cspProfile, queriedProfile, knownProfiles)
  # remove allele with 0 observations in database
  alleleDb = adjust.frequencies( alleleDb, 
                                 queriedProfile[1, colnames(cspProfile), 
                                                drop=FALSE],
                                 adj=adj, fst=fst )

  #alleleDb2 = sapply(alleleDb,FUN=fill.unknown.LUS,simplify=FALSE)

  # combine rare alleles into a single allele
  if(combineRare) alleleDb = combine.rares.peaks(alleleDb, cspProfile, knownProfiles, queriedProfile[1, colnames(cspProfile), drop=FALSE], rareThreshold,doDoubleStutter)

  # add index for stutter values
  #alleleDb = sapply(alleleDb,FUN=add.stutter.index,simplify=FALSE)


  # Construct all profiles as arrays of  
  list(binaryProfile=cspProfile,
       alleleDb=alleleDb,
       queriedProfile=queriedProfile[1, colnames(cspProfile), drop=FALSE],
       knownProfs=knownProfiles) 
}

# function to automatically make allele calls
make.allelic.calls = function(peaksProfile,stutterPercent=0.15)
    {
    out = peaksProfile$alleles
    for(i in 1:length(out))
        {
        tmpAlleles = out[[i]]
        tmpHeights = peaksProfile$heights[[i]]
        for(j in 1:nrow(tmpAlleles))
            {
            for(k in 1:ncol(tmpHeights))
               {
               if(is.na(tmpAlleles[j,k])) break;
               if(k==ncol(tmpHeights)) tmpAlleles[j,k] = 1;
               parentIndex = which(round(as.numeric(tmpAlleles[j,]),1)==round(as.numeric(tmpAlleles[j,k])+1,1))
               if(length(parentIndex)==0)
                    {
                    tmpHeights[j,k] = 1
                    next
                    }
               if(tmpHeights[j,k]<stutterPercent*tmpHeights[j,parentIndex])
                    {
                    tmpHeights[j,k] = 0
                    } else {
                    tmpHeights[j,k] = 1
                    }
               }
            
            }
        out[[i]] = tmpHeights

        }
    return(out)
    }

convertRelationship = function(relationship)
	{
	if(relationship==0) return(c(0,0))
	if(relationship==1) return(c(1,0))
	if(relationship==2) return(c(0.5,0.25))
	if(relationship==3) return(c(0.5,0))
	if(relationship==4) return(c(0.25,0))
	if(relationship==5) return(c(0.25,0))
	if(relationship==6) return(c(0.5,0))
	if(relationship==7) return(c(0.5,0))
	}


# Creates prosecution hypothesis.
# Documentation is in man directory.
prosecution.hypothesis.peaks <- function(peaksFile, refFile, ethnic='NDU1',
                                   nUnknowns=0, adj=1e0, fst=0.03, linkageFile=NULL,
                                   databaseFile=NULL, detectionThresh=20, 
                                   doDropin=FALSE, doDoubleStutter=TRUE,doOverStutter=TRUE, combineRare=TRUE, rareThreshold=1, kit=NULL, relationship=0,...) {
  if(relationship>7|relationship<0) stop("Relationship must be specified between 0 and 7")
  if(is.null(databaseFile)&is.null(kit)) kit = "DNA17"
  relatedness = convertRelationship(relationship)
  alleleDb = load.allele.database(databaseFile,kit)
  peaksProfile = read.peaks.profile(peaksFile)
  cspProfile = t(sapply(peaksProfile$alleles,FUN=convert.to.binary))
  #cspProfile = t(sapply(cspProfile,FUN=function(x) sapply(x,FUN=unlist,simplify=FALSE)))
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

  result = agnostic.hypothesis.peaks(cspProfile, knownProfiles,
                               queriedProfile, alleleDb, ethnic=ethnic,
                               adj=adj, fst=fst, combineRare=combineRare,
			       rareThreshold=rareThreshold, doDoubleStutter=doDoubleStutter)


  	# refIndiv is queried
  	nameRef = rownames(result$queriedProfile)
  	result[["refIndiv"]] = which(rownames(result$knownProfs)==nameRef)


  result = append(result, list(...))
  result[["peaksProfile"]] = peaksProfile$alleles
  result[["heightsProfile"]] = peaksProfile$heights
  result[["knownProfs"]] = knownProfiles
  result[["nUnknowns"]] = nUnknowns
  result[["hypothesis"]] = "prosecution"
  result[["ethnic"]] = ethnic  
  result[["adj"]] = adj
  result[["fst"]] = fst
  result[["relatedness"]] = relatedness
  result[["relationship"]] = relationship 
  result[["doDropin"]] = doDropin
  result[["doDoubleStutter"]] = doDoubleStutter
  result[["doOverStutter"]] = doOverStutter
  result[["peaksFile"]] = peaksFile
  result[["refFile"]] = refFile
  result[["databaseFile"]] = databaseFile
  result[["kit"]] = kit
  result[["detectionThresh"]] = detectionThresh
if(!relationship%in%c(0,1))
	{
    result[["linkageInfo"]] = linkageInfo
    }

  sanity.check.peaks(result) # makes sure hypothesis has right type.
  result
}


# Creates defence hypothesis
# Documentation is in man directory.
defence.hypothesis.peaks <- function(peaksFile, refFile, ethnic='NDU1',  nUnknowns=0,
                               adj=1e0, fst=0.03, databaseFile=NULL, linkageFile=NULL,
                               detectionThresh=20, doDropin=FALSE, doDoubleStutter=TRUE,
				doOverStutter=TRUE, combineRare=TRUE, 
			       rareThreshold=1, kit=NULL, relationship=0,...) {
  if(relationship>7|relationship<0) stop("Relationship must be specified between 0 and 7")
  if(is.null(databaseFile)&is.null(kit)) kit = "DNA17"
  relatedness = convertRelationship(relationship)
  alleleDb = load.allele.database(databaseFile,kit)
  peaksProfile = read.peaks.profile(peaksFile)
#  cspProfile = mapply(convert.to.binary,data=peaksProfile$alleles,allelicCalls=allelicCalls,SIMPLIFY=FALSE)
  cspProfile = t(sapply(peaksProfile$alleles,FUN=convert.to.binary))
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

  result = agnostic.hypothesis.peaks(cspProfile, knownProfiles,
                               queriedProfile, alleleDb, ethnic=ethnic,
                               adj=adj, fst=fst, combineRare=combineRare,
			       rareThreshold=rareThreshold, doDoubleStutter=doDoubleStutter) 
   
  result[["refIndiv"]] = 1

  result = append(result, list(...))
  result[["peaksProfile"]] = peaksProfile$alleles
  result[["heightsProfile"]] = peaksProfile$heights
  result[["knownProfs"]] = knownProfiles
  result[["nUnknowns"]] = nUnknowns + 1
  result[["hypothesis"]] = "defence"
  result[["ethnic"]] = ethnic 
  result[["adj"]] = adj
  result[["fst"]] = fst
  result[["relatedness"]] = relatedness
  result[["relationship"]] = relationship 
  result[["doDropin"]] = doDropin
  result[["doDoubleStutter"]] = doDoubleStutter
  result[["doOverStutter"]] = doOverStutter
  result[["peaksFile"]] = peaksFile
  result[["refFile"]] = refFile
  result[["databaseFile"]] = databaseFile
  result[["kit"]] = kit
  result[["detectionThresh"]] = detectionThresh
if(!relationship%in%c(0,1)) 
	{
	result[["linkageInfo"]] = linkageInfo  	
	}

  sanity.check.peaks(result) # makes sure hypothesis has right type.
  result
}


missing.alleles.peaks = function(alleleDb, cspProfile, queriedProfile, knownProfiles) {
  # Adds missing alleles to database
  #
  # There may be alleles in the crime-scene profile, in the queried individual,
  # or in the individuals subject to dropout, which are not present in the
  # database. Theses alleles are added into the database with count 1 and
  # fragment length 0.
  if(!is.matrix(cspProfile)) stop("input should be a matrix.")
  if(is.null(colnames(cspProfile)))
    stop("input matrix does not have column names.")
  if(!is.matrix(knownProfiles)) stop("input should be a matrix.")
  if(!is.matrix(queriedProfile)) stop("input should be a matrix.")
  for(locus in colnames(cspProfile)) {
    dbAlleles = rownames(alleleDb[[locus]])
    cspAlleles = unique(unlist(cspProfile[, locus]))
    qAlleles = unique(unlist(queriedProfile[, locus]))
    knownAlleles = unique(unlist(knownProfiles[, locus]))
    allAlleles = unique(c(cspAlleles,qAlleles, knownAlleles)) 
    missingAlleles = allAlleles[!allAlleles%in%dbAlleles]
    # At this point, some of the missingAlleles might just be saying move
    # along, nothing to see.
    missingAlleles = setdiff(missingAlleles, c("", NA, "NA"))
    # Now add new rows to database. 
    if(length(missingAlleles)) {
      newrows = matrix(c(1, 0, NA), nrow=length(missingAlleles), ncol=3,
                       byrow=TRUE)
      #newrows = matrix(c(1, 0), nrow=length(missingAlleles), ncol=2,
      #                 byrow=TRUE)
      rownames(newrows) = missingAlleles
      beforeLength = nrow(alleleDb[[locus]])
      addedRows = (beforeLength+1):(beforeLength+nrow(newrows))
      alleleDb[[locus]] = rbind(alleleDb[[locus]], newrows)
      # get missing LUS values
      alleleDb[[locus]][addedRows,3] = sapply(addedRows, FUN=function(x) extrapolateLUS(alleleDb[[locus]],x))
      # get missing BP values
      alleleDb[[locus]][addedRows,2] = sapply(missingAlleles, FUN=function(x) extrapolateBP(alleleDb[[locus]],x))
    }
  }
  alleleDb
}


# Should check that hypothesis is syntactically correct
# This means matrices as opposed to lists, so on and so forth.
# It does not mean it the input should make any sort of sense.
sanity.check.peaks = function(hypothesis) {
  errors = c()
  if(!is.matrix(hypothesis$binaryProfile)) 
     errors = c(errors, "binaryProfile should be a matrix 'replicates' x 'loci'.")
  if(!is.matrix(hypothesis$queriedProfile))
     errors = c(errors, "queriedProfile should be a matrix 'profiles' x 'loci'.")
  if(ncol(hypothesis$binaryProfile) != ncol(hypothesis$queriedProfile))
      errors = c( errors,
                  "Number of loci of binaryProf and queriedProf do not match." )
  if(any(unlist(hypothesis$detectionThresh)<0))
        errors = c( errors,
            "Detection threshold should not be set lower than 0." )
if(length(hypothesis$detectionThresh)==1)
	{
  if(any(as.numeric(unlist(hypothesis$heightsProfile))<hypothesis$detectionThresh,na.rm=TRUE))
	errors = c(errors,
    "Observed peak height below detection threshold.")
	} else {
	loci = colnames(hypothesis$binaryProfile)
	tmp = sapply(loci,FUN=function(x) 
		sapply(1:length(hypothesis$heightsProfile), FUN=function(y) 
			if(any(na.omit(unlist(hypothesis$heightsProfile[[y]][x,]))<
				hypothesis$detectionThresh[[x]])) 
    					paste0("Observed peak height below detection threshold: ", 
					x, " replicate ", y)
			)
		)
	errors = c(errors, unlist(tmp))
	}
  if(any(unlist(hypothesis$detectionThresh)>100))
	errors = c(errors,
    "Detection threshold set unusually high.")
  if(!length(hypothesis$detectionThresh)%in%c(1,ncol(hypothesis$binaryProfile)))
	errors = c(errors,
    "Detection threshold should be a vector with length 1 or a list with length ncol.")
  if(any(hypothesis$relatedness<0))
	errors = c(errors,
	"IBD probability cannot be < 0.")
  if(any(hypothesis$relatedness>1))
	errors = c(errors,
	"IBD probability cannot be > 1.")
  if(sum(hypothesis$relatedness)>1)
	errors = c(errors,
	"Combined IBD probability cannot be > 1.")
  if(hypothesis$fst<0)
	errors = c(errors,
	"Fst cannot be < 0.")
  if(hypothesis$fst>1)
	errors = c(errors,
	"Fst cannot be > 1.")
  if(hypothesis$nUnknowns<0)
	errors = c(errors,
	"nUnknowns cannot be negative.")
  if(hypothesis$nUnknowns>3)
	errors = c(errors,
	"likeLTD does not support more than 3 unknown contributors.")
  if(length(hypothesis$peaksProfile)!=length(hypothesis$heightsProfile))
	errors = c(errors,
	"Different number of replicates in allele profile and heights profile.")
  if(!identical(sapply(hypothesis$peaksProfile,nrow),sapply(hypothesis$heightsProfile,nrow)))
	errors = c(errors,
	"Different number of loci between allele profile and heights profile.")	
  cspInfo = do.call(Map,c(rbind,mapply(hypothesis$peaksProfile,hypothesis$heightsProfile,FUN=function(a,b) 
    list(peaksObs = sapply(1:nrow(a), FUN = function(x) which(!is.na(a[x,])),simplify=FALSE),
      peaksNa = sapply(1:nrow(a), FUN = function(x) which(is.na(a[x,])),simplify=FALSE),
      heightsObs = sapply(1:nrow(b), FUN = function(x) which(!is.na(b[x,])),simplify=FALSE),
      heightsNa = sapply(1:nrow(b), FUN = function(x) which(is.na(b[x,])),simplify=FALSE))
	,SIMPLIFY=FALSE)))
  if(!identical(cspInfo$peaksObs,cspInfo$heightsObs))
	errors = c(errors,
	"Mismatch of data points between alleles and heights.")
  if(!identical(cspInfo$peaksNa,cspInfo$heightsNa))
	errors = c(errors,
	"Mismatch of blank spaces between alleles and heights.")	
  if(any(sapply(1:nrow(cspInfo$peaksObs), FUN=function(a) sapply(1:ncol(cspInfo$peaksObs), 
	FUN=function(b) ifelse(length(unlist(cspInfo$peaksNa[a,b])==0)|length(unlist(cspInfo$peaksObs[a,b])==0),
		FALSE, any(sapply(unlist(cspInfo$peaksNa[a,b]),
		FUN=function(x) x<unlist(cspInfo$peaksObs[a,b]))))))))
	errors = c(errors,
	"Blank space in alleles data.")	
 if(any(sapply(1:nrow(cspInfo$heightsObs), FUN=function(a) sapply(1:ncol(cspInfo$heightsObs), 
	FUN=function(b) ifelse(length(unlist(cspInfo$heightsNa[a,b])==0)|length(unlist(cspInfo$heightsObs[a,b])==0),
		FALSE, any(sapply(unlist(cspInfo$heightsNa[a,b]),
		FUN=function(x) x<unlist(cspInfo$heightsObs[a,b]))))))))
	errors = c(errors,
	"Blank space in heights data.")	
  if(any(is.na(as.numeric(na.omit(unlist(hypothesis$peaksProfile))))))
	errors = c(errors,
	"Non-numeric values in alleles data.")	
  if(any(is.na(as.numeric(na.omit(unlist(hypothesis$heightsProfile))))))
	errors = c(errors,
	"Non-numeric values in heights data.")		
  if(length(errors) != 0) {
    cat("There seems to be an error with the input hypothesis.\n")
    for(error in errors) cat(paste(c(error, "\n")))
    stop()
  }
}


transform.to.locus.centric.peaks = function(hypothesis) {
  # Transform hypothesis to locus centric modes.
  # 
  # This means we reorganise the data to be in lists of the loci.
  loci =  colnames(hypothesis$binaryProfile)
  if(length(hypothesis$detectionThresh)==1) 
	{
	hypothesis$detectionThresh = rep(list(hypothesis$detectionThresh),times=length(loci))
	names(hypothesis$detectionThresh) = loci
	}
  result = list()
  for(locus in loci) {
    # Value of the resulting list for a given locus 
    locusValue = list()
    # Loop over all items in original list.
    for(key in names(hypothesis)) {
      value = hypothesis[[key]]
      # special treatment for peaks profiles (peaks, heights)
      if(key%in%c("peaksProfile","heightsProfile"))
        {
        locusValue[[key]] = sapply(value,FUN=function(x) x[locus, , drop=FALSE],simplify=FALSE)
	    for(i in 1:length(locusValue[[key]]))
    		{
    		tmp = as.numeric(locusValue[[key]][[i]][!is.na(locusValue[[key]][[i]])])
    		if(length(tmp)>0) names(tmp) = paste("Allele",1:length(tmp),sep=".")
    		locusValue[[key]][[i]] = tmp
    		}
        } else {
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
        }
      # If has not been added yet, then not locus dependent, and add it as a
      # whole. 
      if(!key %in% names(locusValue)) locusValue[[key]] = value 
    }
    # Only locus to result if not empty.
    if(length(locusValue) > 0) result[[locus]] = locusValue
  }
  result
}

