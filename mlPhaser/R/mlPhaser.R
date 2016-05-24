#' mlPhaser
#'
#' \tabular{ll}{
#' Package: \tab mlPhaser\cr
#' Type: \tab Package\cr
#' Version: \tab 0.01\cr
#' Date: \tab 2012-08-28\cr
#' Author: \tab Dave T. Gerrard <david.gerrard@@manchester.ac.uk>\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' mlPhaser
#'
#' Some text
#'
#' The main functions are: 
#' \code{\link{phaseReport}}, \code{\link{getValidHaploGroups}},  
#' 
#' ...
#' 
#' @docType package
#' @name mlPhaser-package
NULL

## Simulates genotypes from a table of haplotypes. Will use frequencies.
#' Simulate genotypes
#' 
#' Simulates genotypes from a table of haplotypes.
#'
#' Simulates n genotypes from a table of haplotypes. Genotypes will include one allele per 
#' ploidy level. 
#' 
#' @param haploTable The list of haplotypes in table format
#' @param haploFreqs A named vector of haplotype frequencies.
#' @param n	How many genotypes to simulate. 
#' @param  ploidy How many alleles per locus. Default = 2. 
#' 
#' @return A data.frame of genotypes. Each locus will have multiple colums as per the ploidy level.
#'
#' @export
#' @docType methods
#' @examples
#' # create a data frame to store alleles of haplotypes. Columns are loci.
#' haplotypes <- data.frame(	A= c("a","b","c","a","b","c","b"),
#' 					B= c("a","b","c","b","c","a","a"),
#' 					C= c("a","b","c","b","c","a","a") )
#' # give the haplotypes sensible names as rownames. 
#' rownames(haplotypes) <- apply(haplotypes, 1, paste,sep="" , collapse="")
#' # Create a named vector of haplotype frequencies.
#' haploFreqs <- c(0.4, 0.3, 0.15, 0.07,0.05, 0.02, 0.01)
#' names(haploFreqs) <- rownames(haplotypes)
#'
#' # simulate a set of genotypes
#' my.genotypes <- simGenoFromHaplo(haploTable=haplotypes, haploFreqs=haploFreqs , n=20, ploidy=2) 
simGenoFromHaplo <- function(haploTable, haploFreqs, n=1, ploidy=2)  {
	
	## put in checks that the haplotypes in haploTable have freqs in haploFreqs


	#haplos <- rep(names(haploFreqs), haploFreqs)

	sampleHaplos <- sample(names(haploFreqs), n*ploidy, replace=T, prob=haploFreqs )

	genoTable <- data.frame()

	for(i in 1:ploidy)  {
		newCols <- haploTable[sampleHaplos[((i-1)*n)+1:n],]
		names(newCols) <- paste(names(newCols),i,sep=".")
		if(i==1) {
			genoTable <- newCols
		} else {
			genoTable <- cbind(genoTable, newCols)
		}
	}
		
	rownames(genoTable) <- 1:n
	#return(sampleHaplos)
	genoTable <- genoTable[,order(colnames(genoTable))]		# should columns be re-ordered?
	return(genoTable)
}




#####  this is based only on using known haplotypes. 
##  Each valid haplotype is found by intersecting across each locus.
### What about if there are no valid haplotypes?  NA?
### Account for novel haplotypes?
## Could have optional call to list all possible haplotypes. (and assign very low probabilities to novel forms).
## This function relies on grep of loci names. Not good idea when loci names are "A1", "A10" etc. 
#' List valid haplotypes
#' 
#' Select valid haplotypes (from a list) for a genotype
#'
#' Tests a set of haplotypes against a genotype to see if they are contained within.
#' 
#' @param genotype The genotype. Can be in data.frame or list of list format.
#' @param haplotypes The haplotypes. Can be in data.frame or list of lists format.
#' @param ploidy How many alleles per locus. Default=2.
#' 
#' @return A list of valid haplotypes (each as a list).
#'
#' @export
#' @docType methods
listValidHaplotypes <- function(genotype,haplotypes, ploidy=2)  {
	# could generate haplotype list from fresh. Might be very long.	
	
	# TEST IF genotype and haplotypeList are list() or data.frame().

	if(class(haplotypes) == "data.frame")  {
		haploList <- tableHaploToList(haplotypes)
		haploTable <- haplotypes
	} else if(class(haplotypes) == "list") {
		haploList <- haplotypes
		haploTable <- listHaploToTable(haplotypes)
	} else {
		stop("bad format for haplotypes")
	}

	#lociNames <- c("A", "B", "C")	# TODO find these from genotype
	lociNames <- names(haploList[[1]])	# currently finds lociNames from first entry in haplotypeList

	if(class(genotype) == "data.frame")  {
		genoTable <- genotype
		genoList <- tableGenoToList(genotype, lociNames)
	} else if (class(genotype) == "list") {
		genoList <- genotype
		genoTable  <- listGenoToTable(genotype)
	} else {
		stop("bad format for genotype")
	}




	# NEW
	validHaplos <- character()
	for(i in 1:length(lociNames))  {
		thisLocus <- lociNames[i]
		genotypes <- genoList[[thisLocus]]
		theseHaplos <- character()
		for(thisGenotype in genotypes) {
			theseHaplos <- c(theseHaplos, rownames(haploTable)[haploTable[,thisLocus] == thisGenotype])
		}
		if(i == 1)  {
			validHaplos <- theseHaplos
		} else  {
			validHaplos <- intersect(validHaplos,theseHaplos)
		}
		
	}

	validHaploList  <- haploList[validHaplos]

	return(validHaploList)
}





# Tries to extract a single haplotype from a compound genotype and return, amongst other things, the remainder genotype.
# this version to test and return if remaining alleles.
#' Extract haplotype from genotype
#' 
#' Attempts to extract a haplotype from a genotype
#'
#' Tries to extract a single haplotype from a compound genotype and return, amongst other things,
#' the remainder genotype.
#' 
#' @param haplotype The haplotype to be removed
#' @param genotypeList The genotype in list of lists format.
#' 
#' @return A list giving the original haplotype extracted (haplotype), 
#' a table of TRUE/FALSE for each locus with TRUE if the allele was successfully extracted (passTable),
#' and a list giving the genotype remaining after extraction (remList). 
#'
#' @export
#' @docType methods
remGeno <- function(haplotype, genotypeList)  {
	# haplotype should be a list, preferably named.
	# and entries as alleles at each locus.
	# genotype should be formatted as a list with loci as entries and present alleles under each entry.
	# can use tableGenoToList() to format.
	
	passTable <- logical(length=length(haplotype[[1]]))
	names(passTable) <- names(haplotype[[1]])
	remList <- list()
	for(thisLocus in names(haplotype[[1]])) {
		matchIndex <-  match(haplotype[[1]][[thisLocus]] , genotypeList[[thisLocus]])
		if(is.na(matchIndex))  {
			passTable[thisLocus] <- FALSE
		} else {
			passTable[thisLocus] <- TRUE
			remList[[thisLocus]] <- genotypeList[[thisLocus]][-which(genotypeList[[thisLocus]] == haplotype[[1]][[thisLocus]])[1]]
		}
	}
	if(sum(unlist(lapply(remList,length))) == 0) {
		remList <- NULL
		#remList <- list()	# better as empty list?
	}

	return(list(haplotype=haplotype, passTable=passTable,remList = remList))
}

## test if a single haplotype is consistent with a genotype.
#' Test genotype for presence of haplotype
#' 
#' Test if a genotype contains a haplotype.
#'
#' An early implementation to test if a genotype contained a haplotype.
#'	N.B. I don't think this is used anymore.
#' 
#' @param haplotype The haplotype as a one line data.frame
#' @param genotypeList The genotype as a list of lists.
#' 
#' @return TRUE/FALSE if haplotype is present in the genotype
#'
#' @export
#' @docType methods
testHaploInGeno <- function(haplotype, genotypeList)  {
	# haplotype should be one line dataframe with columns as loci names
	# and entries as alleles at each locus.
	# genotype should be formatted as a list with loci as entries and present alleles under each entry.
	# can use tableGenoToList() to format.
	haploPresent <- FALSE
	testList <- logical()
	for(thisLocus in names(genotypeList)) {
		testList <- c(testList,which(genotypeList[[thisLocus]] == haplotype[,thisLocus]) > 0)
	}

	haploPresent <- all(testList)

	return(haploPresent)
}




# perhaps should be called extractHaploFromGeno
# return genoList with certain haplo taken out, error if haplo cannot be taken out.



## recursion
## requires access to a global list to store results:
### NOT THIS ANYMORE  groupStorageEnv <- new.env(parent = baseenv())
#assign("validHaploGroups", list(), globalenv())
## then call function
## then retrieve validGroups.
#get("validHaploGroups", globalenv())
#' Recurse to get haplo groups
#' 
#' Recursive function to extract valid groups of haplotypes explaining a genotype
#'
#' This recursive function subtracts haplotypes from a genotypes to find 'groups' 
#' of haplotypes that can fully explain a genotype. To make the function general and
#' cope with ploidy > 2, I made it recursive. It will keep going until it has run out 
#' of genotype and/or it has run out of valid haplotypes.
#' This should probably stay as an internal function because of its recursive nature.
#' N.B.  requires access to a globally accessible storage variable: validHaploGroups
#' 
#' @param validHaplotypes A list of haplotypes to choose from 
#' @param remGenotype The remaining part of the genotype
#' @param group The valid group to this point.
#' 
#' @return NULL. N.B.  requires access to a globally accessible storage variable: validHaploGroups
#'
#' @docType methods
recurseHaplos <- function(validHaplotypes, remGenotype, group) {
	#print("New recurse")
		#print("Rem geno, try more haplotypes")
		for(i in 1:length(validHaplotypes))  {
			thisHaplotype <- validHaplotypes[i]
			#print(names(validHaplotypes)[i] )
			remGeneList <- remGeno(thisHaplotype , remGenotype)
			#print(remGeneList[['passTable']])
			if(all(remGeneList[['passTable']])) {
				#print("Extracted haplotype")
				#myGroup <- c(group,remGeneList[['haplotype']])
				#myGroup <- c(group,list(thisHaplotype))
				myGroup <- c(group,thisHaplotype)
				if(is.null(remGeneList[['remList']]) )  {
					#print("Valid pair")
					#print(myGroup)
					#groupStore <- list()
					#for(g in 1:length(myGroup))  {
					#	groupStore[[g]] <- myGroup[[g]]
					#}
					#groupStoreList <- list(groupStore)
					groupStoreList <- list(myGroup)
					# store the valid group within a previously set up globally accesible list.	List of list of list.
					#assign("validHaploGroups", c(get("validHaploGroups", groupStorageEnv),groupStoreList), groupStorageEnv)
					assign("validHaploGroups", c(get("validHaploGroups", globalenv()),groupStoreList), globalenv())
					#print("breaking")
					break()	
				} else  {  # further genotypes to extract
					recurseHaplos(validHaplotypes,remGeneList[['remList']],myGroup)
				}
			}
		}
		#print("tried all haplos or broke")
}

##wrapper function to set up and control the recursive search for groups of haplotyps, each of which are consistent with the genotype in question.
#' Get haplo groups for a genotype
#' 
#' Get all valid groups of haplotypes that fully explain a genotype.
#'
#' Wrapper function to set up and control the recursive search for groups of haplotyps, 
#' each of which are consistent with the genotype in question. Makes use of recursion via
#' the function  \code{\link{recurseHaplos}}
#' 
#' @param genotype The genotype in question. Can be data.frame or list of lists format
#' @param haplotypes The set of candidate haplotypes.
#' 
#' @return A list of valid haplotype groups (each itself a list of haplotypes).
#'
#' @export
#' @docType methods
#' @seealso \code{\link{phaseReport}}
#' @examples
#' # create a data frame to store alleles of haplotypes. Columns are loci.
#' haplotypes <- data.frame(	A= c("a","b","c","a","b","c","b"),
#'					B= c("a","b","c","b","c","a","a"),
#'					C= c("a","b","c","b","c","a","a") )
#' # give the haplotypes sensible names as rownames. 
#' rownames(haplotypes) <- apply(haplotypes, 1, paste,sep="" , collapse="")
#' # load a genotype as a table
#' thisGenotype <- data.frame(A.1="a", A.2="b", B.1="a", B.2="b",C.1="a", C.2="b")
#' # find groups of haplotypes as a list of lists
#' my.valid.groups <- getValidHaploGroups(thisGenotype,haplotypes)
#' # look at the list structure of the valid groups list
#' str(my.valid.groups)
#' # see phaseReport() for more friendly function
getValidHaploGroups <- function(genotype, haplotypes)  {

	if(class(haplotypes) == "data.frame")  {
		haploList <- tableHaploToList(haplotypes)
		haploTable <- haplotypes
	} else if(class(haplotypes) == "list") {
		haploList <- haplotypes
		haploTable <- listHaploToTable(haplotypes)
	} else {
		stop("bad format for haplotypes")
	}

	if(class(genotype) == "data.frame")  {
		genoTable <- genotype
		genoList <- tableGenoToList(genotype, names(haploTable))
	} else if (class(genotype) == "list") {
		genoList <- genotype
		genoTable  <- listGenoToTable(genotype)
	} else {
		stop("bad format for genotype")
	}


	
	startGroup=list()
	#startRemGeno <- list('remList'=genoList )
	startRemGeno <- list()
	startRemGeno[['remList']] <- genoList

	haploList.valid <- listValidHaplotypes(genoTable, haploTable )

	if(length(haploList.valid) < 1)  {
		validHaploGroups.nonRedund <- list()
	} else {
		assign("validHaploGroups", list(), globalenv() )
		# then call function
		# then retrieve validGroups.
		recurseHaplos(validHaplotypes=haploList.valid,startRemGeno[['remList']] ,
					group=startGroup )
		validHaploGroups <- get("validHaploGroups", envir=globalenv())
		#rm(groupStorageEnv)
		validHaploGroups.nonRedund <- reduceRedundantList(validHaploGroups )
		rm(validHaploGroups, envir=globalenv())
	
		if(length(validHaploGroups.nonRedund) < 1) {
			validHaploGroups.nonRedund <- list()
		}
	}
	return(validHaploGroups.nonRedund)
}



## first attempt at assigning probabilities/likelihoods to competing haplotype groups.
## prints and summed log-likelihood and reconstituted (exp()) likelihood.
#' Print haplo group probabilities
#' 
#' Print out haplotype groups and their relative probabilities
#'
#' This was a first attempt to order competing haplotype groups. It only prints and does not
#'  return a useable object.
#' 
#' @param namedHaploGroups A list of haplotype groups (each a list) that explain a genotype.
#' @param haploFrequencies A named numeric vector giving the probabilities of haplotypes. 
#'  Names store the haplotype names.
#' 
#' @return Nothing. Prints only
#'
#' @export
#' @docType methods
printHaploProbs <- function(namedHaploGroups, haploFrequencies) {
	for(i in 1:length(namedHaploGroups)) {	
		thisGroup <- namedHaploGroups[[i]]
			sumLikelihood <- 0
			for(thisHaplo in thisGroup)  {
				thisProb <- haploFrequencies[thisHaplo]
				sumLikelihood <- sumLikelihood + log(thisProb)
				cat(thisHaplo)
				cat("/")
			}
			cat("\t")
			cat(sumLikelihood)
			cat("\t")
			cat(exp(sumLikelihood))
			cat("\n")
	}
}

#printHaploProbs(namedGroups,haploFreqs )


# simple function to print out list of list haplotypes, one per line
#' Print haplotype groups
#' 
#' Prints haplotype groups from a list.
#'
#' Works through a list of haplotype groups (quite a complex list structure) 
#' and prints out one group per line. For presentation only, results cannot be re-used.
#' 
#' @param haploListOfLists A list of haplotype groups.
#' 
#' @return Nothing. Just prints.
#'
#' @export
#' @docType methods
printHaploGroups <- function(haploListOfLists)  {
	# assume this is a list of list of lists
	for(i in 1:length(haploListOfLists)) {
		thisGroup <- haploListOfLists[[i]]
		for(j in 1:length(thisGroup)) {
			thisHaplo <- thisGroup[[j]]
				for(k in 1:length(thisHaplo))  {
					cat(thisHaplo[[k]])
					cat(".")
				}
			cat("/")
		}
		cat("\n")

	}

}

#printHaploGroups(result )



# and some way of scoring/ranking them based on relative probabilities.




#' Remove redundant haplotype groups
#' 
#' Removes redundant groups of haplotypes from a common list. 
#'
#' The recursive method \code{\link{recurseHaplos}} of finding groups of consistent haplotypes does not differentiate,
#' re-arranged versions of the same set. e.g. keeps aaa/bbb AND bbb/aaa.
#' This function removes that redundancy from the results.
#' uses the length of intersect to determine if two lists contain all the same elements. 
#' 
#' @param startList A list of haplotype groups (each is a list of haplotypes).
#' 
#' @return  A list of haplotype groups but each group is unique.
#'
#' @export
#' @docType methods
reduceRedundantList <- function(startList)  {	
	listLength <- length(startList)
	if(listLength < 2) {
		return(startList)

	} else {
		removeIndex <- integer()

		for(i in 1:(listLength-1)) {
			if(!is.na(match(i,removeIndex))) {
				next()
			} else {
				for(j in (i+1):listLength)  {
					if(!is.na(match(j,removeIndex))) {	
						next()
					} else {
						if(length(intersect(startList[[i]],startList[[j]])) == length(startList[[j]])) {
							removeIndex <- c(removeIndex, j)
						}
					}
				}
			}		
		}	
		return(startList[-removeIndex])
	}
}



#### TODO this function uses grep to find columns for each locus. Potentially very bad if
#### 		one locus name is a substring of another. e.g. HLA_1, HLA_10
# converts a table genotype with multi columns per locus to a list with one item per locus, each listing the alleles present.
#' Convert genotype table to list of lists
#' 
#' Converts a table of genotypes to a list of lists, one sub-list per genotype.
#'
#' Converts a table of genotypes to a list of lists, one sub-list per genotype.
#' 
#' @param genoTable A data.frame containg genotypes. One row per genotype. Multiple columns per locus
#'   as per the ploidy. 
#' @param locusNames A character vector of locus names.
#' 
#' @return genotypes as a list of lists
#'
#' @export
#' @docType methods
tableGenoToList <- function(genoTable, locusNames)  {
	# need to check if only one row in genoTable
	
	genoList <- list()
	for(thisLocus in locusNames)  {
		locusColumns <- grep(thisLocus, names(genoTable))
		genoList[[thisLocus]] <- as.character(unlist(genoTable[,locusColumns]))
	}
	return(genoList)

}

# converts a table (with rownames) of haplotypes to a list of haplotypes.
#' Haplotype table to haplotype list
#' 
#' Converts a data.frame of haplotypes to a list of lists with haplotypes
#' at the top level and list of loci (with their alleles) beneath.
#'
#' Summary paragraph outlining method
#' 
#' @param haploTable A data.frame of alleles making up the haplotypes. One columm per locus, 
#' 	one row per haplotype. 	The rownames should contain the haplotype ids.
#' @param locusNames A character vector giving the names of the loci which should match
#'  the column names of the haploTable
#' 
#' @return Haplotypes as a list of lists.
#'
#' @export
#' @docType methods
tableHaploToList <- function(haploTable, locusNames=colnames(haploTable))  {
	haploList <- list()
	for(thisRow in 1:nrow(haploTable)) {
		thisList <- list()
		thisHaplo <- rownames(haploTable)[thisRow]
		for(thisLocus in locusNames)  {
			#locusColumns <- grep(thisLocus, names(haploTable))
			thisList[[thisLocus]] <- as.character(unlist(haploTable[thisRow,thisLocus]))
		}
		haploList[[thisHaplo]] <- thisList
	}
	return(haploList )
	
}

# tableHaploToList(haplotypes[listValidHaplotypes(thisGenotype, haplotypes ),])

#' Genotype list to genotype table
#' 
#' Converts a list of genotypes to a table with several columns per locus.
#'
#' The multiple columns per locus are differentiated with a numerical suffic. e.g. locA.1, locA.2
#' 
#' @param genoList A list of lists. One entry per genotype (sample) each containing a list of loci 
#' to store the alleles. 
#' 
#' @return A data.frame a row for each genotype and n colums for each locus (where n is ploidy of locus)
#'
#' @export
#' @docType methods
listGenoToTable <- function(genoList) {
	exportTable <- data.frame()
	for(i in 1:length(genoList)) {
		thisGeno <- genoList[i]
		genoTable <- as.data.frame(thisGeno)
		thisRow <- data.frame(dummy=NA)		# require dummy column to use cbind with first row.
		for(j in 1:nrow(genoTable)) {
			partRow <- genoTable[j,]
			colnames(partRow) <- paste(colnames(genoTable),j,sep=".")
			#print(partRow)
			thisRow <- cbind(thisRow,partRow)
			#print(thisRow)
			
		}
		thisRow <- subset(thisRow, select=-dummy)	# need to remove dummy column again.
		exportTable <- rbind(exportTable,thisRow)
	}
	exportTable <- exportTable[,order(names(exportTable))]

	return(exportTable)
}

#listGenoToTable(list(tableGenoToList(thisGenotype, locusNames=c("A","B","C"))))



#' Haplotype list to Haplotype table.
#' 
#' Converts a set of haplotypes from list format to table format
#'
#' Each list element becomes a row. Each locus becomes a column.
#' 
#' @param haploList A list of lists. One top level element per haplotype. 
#' Each haplotype should be named and have the same set of loci as a sublist. 
#' 
#' @return A set of haplotypes in table format as a data.frame
#'
#' @export
#' @docType methods
listHaploToTable <- function(haploList)  {
	exportTable <- data.frame()
	for(i in 1:length(haploList)) {
		thisHaplo  <- haploList[i]
		thisRow <- data.frame(thisHaplo[[1]])
		rownames(thisRow) <- names(thisHaplo)
		exportTable <- rbind(exportTable , thisRow)
	}
	return(exportTable)
}

#listHaploToTable(fullHaploList)

#' Get haplotype group probability
#' 
#' Combine probabilities across a haplotype group
#'
#' Each haplotype might have a probability of being found. e.g. the population frequency.
#' This function combines probabilities of a group of haplotypes
#' 
#' @param haploGroup A list of haplotypes. Only the names attribute is important.
#' @param haploFreqs The frequencies of haplotypes as a named vector.
#' @param method Only one method currently ("basic").
#' @param returnLog Whether to output the combined probability as a natural log. Default=FALSE.
#' 
#' @return A numeric likelihood representing the combined probability across a group.
#'
#' @export
#' @docType methods
getHaploGroupProb <- function(haploGroup, haploFreqs, method="basic", returnLog=FALSE)  {
	score <- sum(log(haploFreqs[names(haploGroup)]))
	if(!returnLog) {
		score <- exp(score)
	}
	return(score)
}

#' Best/all hapltype groups for a genotype
#' 
#' Attempts to find best/all haplotype groups that fully explain observed
#' multi-locus genotypes.
#'
#' This wrapper function takes a set of genotypes, a set of haplotypes 
#' and a set of haplotype frequencies and attempts to report either all
#' groups or just the single most likely group of known haplotypes that
#' fully explains each observed genotype.
#' 
#' @param genotypes The table/list of genotypes
#' @param haplotypes The table/list of candidate haplotypes
#' @param haploFreqs The frequencies of haplotypes as a named vector.
#' @param outFormat Whether to output all valid haplotype groups or just the best (based on joint probability).
#' 
#' @return A data.frame of results...
#'
#' @export
#' @docType methods
#' @seealso \code{\link{getValidHaploGroups}}
#' @examples
#' # create a data frame to store alleles of haplotypes. Columns are loci.
#' haplotypes <- data.frame(	A= c("a","b","c","a","b","c","b"),
#' 					B= c("a","b","c","b","c","a","a"),
#' 					C= c("a","b","c","b","c","a","a") )
#' # give the haplotypes sensible names as rownames. 
#' rownames(haplotypes) <- apply(haplotypes, 1, paste,sep="" , collapse="")
#' # Create a named vector of haplotype frequencies.
#' haploFreqs <- c(0.4, 0.3, 0.15, 0.07,0.05, 0.02, 0.01)
#' names(haploFreqs) <- rownames(haplotypes)
#' # load a genotype as a table
#' thisGenotype <- data.frame(A.1="a", A.2="b", B.1="a", B.2="b",C.1="a", C.2="b")
#' phaseReport(thisGenotype,haplotypes)
#' # use haplotype frequencies to rank candidate haplotype groups.
#' phaseReport(thisGenotype,haplotypes, haploFreqs)
#' # return only the best haplotype group for each genotype.
#' phaseReport(thisGenotype,haplotypes, haploFreqs, outFormat="top")
#'
#' # simulate a set of genootypes
#' my.genotypes <- simGenoFromHaplo(haploTable=haplotypes, haploFreqs=haploFreqs , n=20, ploidy=2) 
#' # get phase report on all genotypes
#' phaseReport(my.genotypes,haplotypes, haploFreqs, outFormat="all")	# outFormat="all" is the default
#' phaseReport(my.genotypes,haplotypes, haploFreqs, outFormat="top")
phaseReport <- function(genotypes,haplotypes,haploFreqs,outFormat="all")   {
	## do conversion for genotypes and haplotypes.
	genoTable <- genotypes
	haploTable <- haplotypes	

	phaseResults <- data.frame()
	for(i in 1:nrow(genoTable)) {
		thisGenotype <- genoTable[i,]
		genoName <- rownames(thisGenotype)
		result <- getValidHaploGroups(thisGenotype , haplotypes=haploTable)
		if(length(result) < 1)  {
			thisRow <- data.frame(id=genoName,n.validGroups = length(result),haploGroup=NA, likelihood=NA)
			phaseResults <- rbind(phaseResults, thisRow)
		} else {
			for(j in 1:length(result))  {
				if(missing(haploFreqs))  {
					thisLH=NA	
				} else {
					thisLH <- getHaploGroupProb(result[[j]], haploFreqs=haploFreqs)
				}
				thisRow <- data.frame(id=genoName,n.validGroups = length(result),haploGroup=paste(sort(names(result[[j]])), collapse="/"), likelihood=thisLH)
				phaseResults <- rbind(phaseResults, thisRow)
			}
		}
			
	}
	
	phaseResults <- phaseResults[order(phaseResults$id,phaseResults$likelihood, decreasing=T),]
	if(outFormat=="top")  {
		phaseResults <- phaseResults[match(unique(phaseResults$id),phaseResults$id),]
	}
	phaseResults <- phaseResults[order(phaseResults$id),]	# need second reorder to retrieve sample order.
	return(phaseResults )
}

#test <- phaseReport(genotypes=my.genotypes,haplotypes)
#test <- phaseReport(genotypes=my.genotypes,haplotypes, outFormat="top")