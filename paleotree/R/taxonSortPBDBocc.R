#' Sorting Unique Taxa of a Given Rank from Paleobiology Database Occurrence Data
#'
#' Functions for sorting out unique taxa from Paleobiology Database occurrence downloads,
#' which should accept several different formats resulting from different versions of the
#' PBDB API and different vocabularies available from the API.

#' @details
#' Data input for \code{taxonSortPBDBocc} are expected to be from version 1.2 API
#' with the 'pbdb' vocabulary. However, datasets are passed to internal function \code{translatePBDBocc},
#' which attempts to correct any necessary field names and field contents used by
#' \code{taxonSortPBDBocc}.
#'
#' This function can pull either \emph{just} the 'formally' identified and synonymized taxa in a given table of occurrence
#' data or pull \emph{in addition} occurrences listed under informal taxa of the sought taxonomic rank. Only formal taxa
#' are sorted by default; this is controlled by argument \code{onlyFormal}. Pulling the informally-listed taxonomic
#' occurrences is often necessary in some groups that have received little focused taxonomic effort, such that many
#' species are linked to their generic taxon ID and never received a species-level taxonomic ID in the PBDB.
#' Pulling both formal and informally listed taxonomic occurrences is a hierarchical process and performed in
#' stages: formal taxa are identified first, informal taxa are identified from the occurrences that are
#' 'leftover', and informal occurrences with name labels that match a previously sorted formally listed
#' taxon are concatenated to the 'formal' occurrences for that same taxon, rather than
#' being listed under separate elements of the list as if they were separate taxa.

#' This function is simpler than similar functions that inspired it by using the input"rank" to both filter occurrences
#' and directly reference a taxon's accepted taxonomic placement, rather than a series of specific \code{if()} checks.
#'
#' Unlike some similar functions in other packages, such as version 0.3 \code{paleobioDB}'s \code{pbdb_temp_range},
#' \code{taxonSortPBDBocc} does not check if sorted taxa have a single 'taxon_no' ID number. This makes the blanket 
#' assumption that if a taxon's listed name in relevant fields is identical, the taxon is identical, with the important
#' caveat that occurrences with accepted formal synonymies are sorted first based on their accepted names, followed by
#' taxa without formal taxon IDs. This should avoid mistakingly linking the same occurrences to multiple taxa or assigning
#' occurrences listed under separate formal taxa to the same taxon based on their 'identified' taxon name, as long as all
#' formal taxa have unique names (which is an untested assumption). In some cases, this procedure is helpful, such as when
#' taxa with identical generic and species names are listed under separate taxon ID numbers because of a difference in the
#' listed subgenus for some occurrences (example, "Pseudoclimacograptus (Metaclimacograptus) hughesi' and
#' 'Pseudoclimacograptus hughesi' in the PBDB as of 03/01/2015). Presumably any data that would be affected by differences
#' in this procedure is very minor.
#'
#' Occurrences with taxonomic uncertainty indicators in the listed identified taxon name are removed
#' by default, as controlled by argument \code{cleanUncertain}. This is done by removing any occurrences that
#' have an entry in \code{primary_reso} (was "\code{genus_reso}" in v1.1 API) when \code{rank} is a
#' supraspecific level, and \code{species_reso} when rank=species, if that entry is not found in
#' \code{cleanResoValues}. In some rare cases, when \code{onlyFormal=FALSE}, supraspecific taxon names may be
#' returned in the output that have various 'cruft' attached, like 'n.sp'.
#'
#' Empty values in the input data table ("") are converted to NAs, as they may be due to issues
#' with using read.csv to convert API-downloaded data.

#' @param data A table of occurrence data collected from the Paleobiology Database. 

#' @param rank The selected taxon rank; must be one of 'species', 'genus', 'family', 'order',
#' 'class' or 'phylum'.

#' @param onlyFormal If TRUE (the default) only taxa formally accepted by the Paleobiology
#' Database are returned. If FALSE, then the identified name fields are searched for any
#' additional 'informal' taxa with the proper taxon. If their taxon name happens to match
#' any formal taxa, their occurrences are merged onto the formal taxa. This argument generally
#' has any appreciable effect when rank=species.

#' @param cleanUncertain If TRUE (the default) any occurrences with an entry in the respective
#' 'resolution' field that is *not* found in the argument cleanResoValue will be removed from
#' the dataset. These are assumed to be values indicating taxonomic uncertainty, i.e. 'cf.' or '?'.

#' @param cleanResoValues The set of values that can be found in a 'resolution' field that do not
#' cause a taxon to be removed, as they do not seem to indicate taxonomic uncertainty.

#' @return
#' Returns a list where each element is different unique taxon obtained by the sorting function,
#' and named with that taxon name. Each element is composed of a table containing all the same
#' occurrence data fields as the input (potentially with some fields renamed and some field
#' contents change, due to vocabulary translation).
	
#' @seealso
#' \code{\link{occData2timeList}}, \code{\link{plotOccData}} and the
#' example graptolite dataset at \code{\link{graptPBDB}}

#' @author 
#' David W. Bapst, but partly inspired by Matthew Clapham's \code{cleanTaxon} 
#' (found at \href{https://github.com/mclapham/PBDB-R-scripts/blob/master/taxonClean.R}{this location} on github) and
#' R package paleobioDB's \code{pbdb_temp_range} function (found
#' at  \href{https://github.com/ropensci/paleobioDB/blob/master/R/pbdb_temporal_functions.R#L64-178 }{this location} 
#' on github.

#' @examples
#' #load example graptolite PBDB occ dataset
#' data(graptPBDB)
#' 
#' #get formal genera
#' occGenus<-taxonSortPBDBocc(graptOccPBDB, rank="genus")
#' length(occGenus)
#' 
#' #get formal species
#' occSpeciesFormal<-taxonSortPBDBocc(graptOccPBDB, rank="species")
#' length(occSpeciesFormal)
#' 
#' #yes, there are fewer 'formal' graptolite species in the PBDB then genera
#' 
#' #get formal and informal species
#' occSpeciesInformal<-taxonSortPBDBocc(graptOccPBDB, rank="species",
#' 	 onlyFormal=FALSE)
#' length(occSpeciesInformal)
#' 
#' #way more graptolite species are 'informal' in the PBDB
#' 
#' #get formal and informal species 
#' 	#including from occurrences with uncertain taxonomy
#' 	#basically everything and the kitchen sink
#' occSpeciesEverything<-taxonSortPBDBocc(graptOccPBDB, rank="species",
#' 		onlyFormal=FALSE, cleanUncertain=FALSE)
#' length(occSpeciesEverything)
#' 
#' \dontrun{
#'
#' # simple function for getting occurrence data from API v1.1 
#' easyGetPBDBocc<-function(taxa,show=c("ident","phylo")){
#'   #cleans PBDB occurrence downloads of warnings
#'   taxa<-paste(taxa,collapse=",")
#' 	taxa<-paste(unlist(strsplit(taxa,"_")),collapse="%20")
#' 	show<-paste(show,collapse=",")
#' 	command<-paste0("http://paleobiodb.org/data1.2/occs/list.txt?base_name=",
#' 		taxa,"&show=",show,"&limit=all",
#' 		collapse="")
#' 	command<-paste(unlist(strsplit(command,split=" ")),collapse="%20")
#' 	downData<-readLines(command)
#' 	if(length(grep("Warning",downData))!=0){
#' 		start<-grep("Records",downData)
#' 		warn<-downData[1:(start-1)]
#' 		warn<-sapply(warn, function(x) 
#' 			paste0(unlist(strsplit(unlist(strsplit(x,'"')),",")),collapse=""))
#' 		warn<-paste0(warn,collapse="\n")
#' 		names(warn)<-NULL
#' 		mat<-downData[-(1:start)]
#' 		mat<-read.csv(textConnection(mat))
#' 		message(warn)
#' 	}else{
#' 		mat<-downData
#' 		mat<-read.csv(textConnection(mat))
#' 		}
#' 	return(mat)
#' 	}
#' 
#' #try a PBDB API download with lots of synonymization
#' 	#this should have only 1 species
#' #old way:
#' #acoData<-read.csv(paste0("http://paleobiodb.org/data1.1/occs/list.txt?",
#' #	"base_name=Acosarina%20minuta&show=ident,phylo&limit=all"))
#' # with easyGetPBDBocc:
#' acoData<-easyGetPBDBocc("Acosarina minuta")
#' x<-taxonSortPBDBocc(acoData, rank="species", onlyFormal=FALSE)
#' names(x)
#'
#' #make sure works with API v1.1
#' dicelloData<-read.csv(paste0("http://paleobiodb.org",
#' 	"/data1.1/occs/list.txt?base_name=Dicellograptus",
#' 	"&show=ident,phylo&limit=all"))
#' dicelloOcc2<-taxonSortPBDBocc(dicelloData, rank="species", onlyFormal=FALSE)
#' names(dicelloOcc2)
#' 
#' #make sure works with compact vocab v1.1
#' dicelloData<-read.csv(paste0("http://paleobiodb.org",
#' 	"/data1.1/occs/list.txt?base_name=Dicellograptus",
#' 	"&show=ident,phylo&limit=all&vocab=com"))
#' dicelloOccCom1<-taxonSortPBDBocc(dicelloData, rank="species", onlyFormal=FALSE)
#' names(dicelloOccCom1)
#' head(dicelloOccCom1[[1]])[,1:7]
#'
#' #make sure works with compact vocab v1.2
#' dicelloData<-read.csv(paste0("http://paleobiodb.org",
#' 	"/data1.2/occs/list.txt?base_name=Dicellograptus",
#' 	"&show=ident,phylo&limit=all&vocab=com"))
#' dicelloOccCom1<-taxonSortPBDBocc(dicelloData, rank="species", onlyFormal=FALSE)
#' names(dicelloOccCom1)
#' head(dicelloOccCom1[[1]])[,1:7]
#'
#' }
#' 

#' @name taxonSortPBDBocc
#' @rdname taxonSortPBDBocc
#' @export
taxonSortPBDBocc<-function(data,rank, onlyFormal=TRUE, cleanUncertain=TRUE, 
								cleanResoValues=c(NA, '"', "", "n. sp.", "n. gen."," ","  ")){
	#this function inspired by Matt Clapham's taxonClean and paleobioDB's pbdb_temp_range
		#onlyFormal=FALSE;rank="species"
		#onlyFormal=FALSE;rank="genus"
	#pull occurrences out of a data table and sort by unique taxa into a list structure
	#translated vocabs!
	data<-translatePBDBocc(data)
	#Second, some warning checks to see if occurrences downloaded correctly: 
	# need phylo and ident data
		#the following is taken with minor modification from code in paleobioDB package
	if (!any("primary_name"==colnames(data)) | !any("genus"==colnames(data))){	
		stop("need to add 'show=c('phylo', 'ident')' to pbdb_occurrences query\n *or*\n  'show=ident,phylo' to PBDB API query")
		}
	#additional checks for rank
	if(length(rank)!=1){stop("length of rank must be 1")}
	if(!any(sapply(c("species","genus","family","order","class","phylum"),function(x) x==rank))){
		stop("rank must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")}
	#for inconsistencies between rank and onlyFormal - NOT TRUE, IGNORE THIS CHECK
	#if(!onlyFormal & (rank!="species" | rank!="genus")){
	#  stop("Informal taxon does not exist for above genus level, please use 'onlyFormal=TRUE'")}
	#need to replace any empty string values with NAs (due perhaps to use of read.csv with the API)
	data[data==""]<-NA
	#now,  pull taxa using the relevant variable from 'phylo' for formal ID
		#this matches what paleobioDB now does
	#sort occurrences by unique taxa in each level and then append valid names
	if(rank=="species"){
		if(cleanUncertain){ #remove uncertain species taxonomy
		data<-data[data[,"species_reso"] %in% cleanResoValues,,drop=FALSE]
		}
	#first formal taxa
    #get species names for taxa formally recognized as species level
    whichFormal<-(data[,"accepted_rank"]==rank)
    taxonVar<-as.character(data[,"accepted_name"])
    taxaNames<-as.character(unique(taxonVar[whichFormal]))  #get unique taxa
    #drop empty entries (occurrences listed at a higher taxonomic level, formally)
    taxaNames<-taxaNames[!is.na(taxaNames)]
    #take rows these taxa are in and drop them into elements of a list
		sortedOcc<-lapply(taxaNames,function(x) data[which(x==taxonVar),,drop=FALSE])
    names(sortedOcc)<-taxaNames
    if(!onlyFormal){
		#now use taxon_rank to identify useful informal occurrence
		taxonVar2<-data[,c("primary_name","species_name")]
		taxonVar2<-apply(taxonVar2,1,paste,collapse=" ")
		#which of these occs are useful as informal occs
		stillUseful<-which(!whichFormal & (data[,"identified_rank"]==rank))
		taxaNames2<-as.character(unique(taxonVar2[stillUseful]))
		#drop any weird empties
		taxaNames2<-taxaNames2[taxaNames2!=" " & !is.na(taxaNames2)]  
		sortedOcc2<-lapply(taxaNames2,function(x) data[which(x==taxonVar2),,drop=FALSE])
		names(sortedOcc2)<-taxaNames2
		# merge sortedocc2 with sortedOcc
		share1<-sapply(taxaNames,function(x) any(x==taxaNames2))
		share2<-sapply(taxaNames2,function(x) any(x==taxaNames))
		if(sum(!share1)>0){sortedOccU<-sortedOcc[!share1]}else{sortedOccU<-list()}
		if(sum(!share2)>0){sortedOcc2U<-sortedOcc2[!share2]}else{sortedOcc2U<-list()}
		#and for shared names
		if(sum(share1)>0){
			shared<-lapply(taxaNames[share1],function(x) 
				cbind(sortedOcc[[taxaNames==x]],sortedOcc2[[taxaNames2==x]]))
			names(shared)<-names(taxaNames[share1])
		}else{shared<-list()}
		sortedOcc<-c(sortedOccU,sortedOcc2U,shared)
		}
	}else{   #if not at the species rank
		if(cleanUncertain){  #removing uncertain taxonomy if appended to primary_name
			data<-data[data[,"primary_reso"] %in% cleanResoValues,,drop=FALSE]
			}
		taxonVar<-data[,rank] #then our taxonomic variable of interest is just so!
		taxaNames<-as.character(unique(taxonVar))  #get unique taxa
		taxaNames<-taxaNames[!is.na(taxaNames)]
		#take rows these taxa are in and drop them into elements of a list
		sortedOcc<-lapply(taxaNames,function(x) data[which(x==taxonVar),,drop=FALSE])
		names(sortedOcc)<-taxaNames
		if(!onlyFormal){
			#now use taxon_rank to identify useful informal occurrence
			taxonVar2<-data[,"primary_name"]
			#which of these occs are useful as informal occs
			stillUseful<-which(is.na(taxonVar) & (data[,"identified_rank"]==rank))
			taxaNames2<-as.character(unique(taxonVar2[stillUseful]))
			taxaNames2<-taxaNames2[!is.na(taxaNames2)]
			sortedOcc2<-lapply(taxaNames2,function(x) data[which(x==taxonVar2),,drop=FALSE])
			names(sortedOcc2)<-taxaNames2
			# merge sortedocc2 with sortedOcc
			share1<-sapply(taxaNames,function(x) any(x==taxaNames2))
			share2<-sapply(taxaNames2,function(x) any(x==taxaNames))
			if(sum(!share1)>0){sortedOccU<-sortedOcc[!share1]}else{sortedOccU<-list()}
			if(sum(!share2)>0){sortedOcc2U<-sortedOcc2[!share2]}else{sortedOcc2U<-list()}
			#and for shared names
			if(sum(share1)>0){
				shared<-lapply(taxaNames[share1],function(x)
					cbind(sortedOcc[[taxaNames==x]],sortedOcc2[[taxaNames2==x]]))
				names(shared)<-names(taxaNames[share1])
			}else{shared<-list()}
			sortedOcc<-c(sortedOccU,sortedOcc2U,shared)
			}
		}
	# sort occurrences by taxon name
	sortedOcc<-sortedOcc[order(names(sortedOcc))]
	return(sortedOcc)
	}

translatePBDBocc<-function(data){
	#translate PBDB occ data	
	if(any(colnames(data)=="taxon_name")){
		#from PBDB API version 1.1 to 1.2
		colnames(data)[colnames(data)=="taxon_name"]<-"identified_name"	
		colnames(data)[colnames(data)=="taxon_no"]<-"identified_no"
		colnames(data)[colnames(data)=="taxon_rank"]<-"identified_rank"	
		colnames(data)[colnames(data)=="matched_name"]<-"accepted_name"
		colnames(data)[colnames(data)=="matched_no"]<-"accepted_no"
		colnames(data)[colnames(data)=="matched_rank"]<-"accepted_rank"
		colnames(data)[colnames(data)=="genus_name"]<-"primary_name"
		colnames(data)[colnames(data)=="genus_reso"]<-"primary_reso"
		}
	if(any("tna"==colnames(data))){
		#need to translate data from 1 vocab to another
			#do on by-function basis - this is for taxonSort
		#for this just translate colname of relevant taxon variables
			#also will need to translate 'taxon_rank' and 'matched_rank' contents
		if(any("mna"==colnames(data))){
			#then compact vocab v1.1
			if(!all(c("tna","rnk","tid","mna","mra","mid","idt","ids",
				"gnl","fml","odl","cll","phl","rst","rss") %in% colnames(data))){
					stop("Not all fields founds for compact vocab v1.1")}	
			colnames(data)[colnames(data)=="tna"]<-"identified_name"	
			colnames(data)[colnames(data)=="tid"]<-"identified_no"
			colnames(data)[colnames(data)=="rnk"]<-"identified_rank"	
			colnames(data)[colnames(data)=="mna"]<-"accepted_name"
			colnames(data)[colnames(data)=="mid"]<-"accepted_no"
			colnames(data)[colnames(data)=="mra"]<-"accepted_rank"
			colnames(data)[colnames(data)=="idt"]<-"primary_name"
			colnames(data)[colnames(data)=="ids"]<-"species_name"
			colnames(data)[colnames(data)=="rss"]<-"primary_reso"
			colnames(data)[colnames(data)=="rst"]<-"species_reso"
			colnames(data)[colnames(data)=="gnl"]<-"genus"
			colnames(data)[colnames(data)=="fml"]<-"family"
			colnames(data)[colnames(data)=="odl"]<-"order"
			colnames(data)[colnames(data)=="cll"]<-"class"	
			colnames(data)[colnames(data)=="phl"]<-"phylum"			
			}
		if(any("idn"==colnames(data))){		
			#then compact vocab v1.2
			if(!all(c(	"idn","idr","iid","tna","rnk","tid","idg","ids",
				"gnl","fml","odl","cll","phl","rsg","rss") %in% colnames(data))){
					stop("Not all fields founds for compact vocab v1.1")}	
			colnames(data)[colnames(data)=="idn"]<-"identified_name"	
			colnames(data)[colnames(data)=="iid"]<-"identified_no"
			colnames(data)[colnames(data)=="idr"]<-"identified_rank"	
			colnames(data)[colnames(data)=="tna"]<-"accepted_name"
			colnames(data)[colnames(data)=="tid"]<-"accepted_no"
			colnames(data)[colnames(data)=="rnk"]<-"accepted_rank"
			colnames(data)[colnames(data)=="idg"]<-"primary_name"
			colnames(data)[colnames(data)=="ids"]<-"species_name"
			colnames(data)[colnames(data)=="rss"]<-"primary_reso"
			colnames(data)[colnames(data)=="rsg"]<-"species_reso"
			colnames(data)[colnames(data)=="gnl"]<-"genus"
			colnames(data)[colnames(data)=="fml"]<-"family"
			colnames(data)[colnames(data)=="odl"]<-"order"
			colnames(data)[colnames(data)=="cll"]<-"class"	
			colnames(data)[colnames(data)=="phl"]<-"phylum"			
			#stop("need to add 'vocab='pbdbd'' to pbdb_occurrences query\n  *or*\n 'vocab=pbdb' to your PBDB API query")
			}
		# taxon rank translation vectors for compact vocab
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")
		taxRankCOM<-2:26
		#change contents of "identified_rank" and "accepted_rank"
		data$identified_rank<-sapply(data$identified_rank,function(x) taxRankPBDB[x==taxRankCOM])
		data$accepted_rank<-sapply(data$accepted_rank,function(x) taxRankPBDB[x==taxRankCOM])
		}
	return(data)
	}
	
	