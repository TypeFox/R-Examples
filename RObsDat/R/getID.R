getID <- function(table, value,#, allowNoValue=FALSE
		  #allowNoValue does not make sense as it is implemented
		  remove.special.character=TRUE
		  ){
	  #return directly if no value is passed
	if(length(value)==0) return(c())
	
	#convert table to string if passed as factor
	table <- as.character(table)
# Generate Table where to search for information - this could be done once only
	lookup <- list(SpatialReference = c("ID", "SRSName", "SRSID"), 
		       Site=c('ID', 'Name','Code'),
		       Method=c('ID','Description'),
		       Qualifier=c('ID','Description','Code'),
		       QualityControlLevel=c('ID','Definition','Explanation', 'Code'),
		       Sample=c('ID','LabSampleCode'),
		       Source=c('ID','Organization','Description','Citation'),
		       Variable=c('ID','Code','Name'),
		       OffsetType=c('ID', 'Description'),
		       Units=c('ID', 'Name', 'Abbreviation'),
		       ISOMetadata=c('ID', 'Title','Abstract')
		       )
	old.names <- names(lookup)
	cv.tables <- CVtables()
	for(i in seq(along=cv.tables)){
		lookup[[i+length(old.names)]] <- c('Term','Definition')
	}
	names(lookup) <- c(old.names, cv.tables)

	#select info for correct table


	if(length(lookup[[table]])==0) stop(paste("No information for table",table,"in function getID"))

	#make sure, value is a vector
	value.orig <- value
	value <- as.vector(value)


	#avoid special characters
	#and remove trailing spaces
	if(remove.special.character){
		#we need to be able to turn this off after lookup with
		#levenstein algorithm, if special characters are present
		value <- sub(' +$', '', iconv(value, "","ascii",sub=" "))
	}

	#remove trailing spaces
	value <- sub('^ ', '', value)

	#get ids for unique values
	uvalue <- unique(value)
	uvalueID <- c()
	for(i in 1:length(uvalue)){
		#make sure numeric values are treated as such
		theUvalue <- uvalue[i]
		phraseA <- gsub("'", "`", theUvalue)
		numUvalue <- suppressWarnings(as.numeric(phraseA))
		if(!is.na(numUvalue)){
			if(numUvalue==phraseA){
				phraseA <- numUvalue
			}
		}
		#aaa
		entry <- c()
		for(exact in c(TRUE,FALSE)){
			#avoid running loop again if we found a value
			if(NROW(entry)==0){
				for(field in lookup[[table]]){
					#Skip ID field if we have characters (causes trouble with postgres
					if(field == "ID" & is.character(phraseA)) next
					command <- paste('entry <- unique(rbind(entry, Iget',table,'(options("odm.handler")[[1]], ',field,'="',phraseA,'", exact=',exact,')))', sep='')
					eval(parse(text=command))
					#stop if we found too many results
					if(NROW(entry)>1){ 
						if(exact){
							#aaa
							synonymID <- IgetSynonymID(getOption("odm.handler"), table=table, phrase=phraseA)
							if(length(synonymID)==1){
								uvalueID[i] <- synonymID
								break
							}
						}
						print(entry)
						cat("Term '",theUvalue,"' returns more than 1 values in table ",table , "\n", sep="")
						cat("\n\n Please select matching row or hit 0 to stop\n")
						if(!interactive()) stop("Error: no unique value in non-interactive session")
						choice <- "impossible"
						inval <- 0
						while(choice == "impossible"){
							next.step <- readline("What is your choice? ")
							choice <- as.numeric(next.step)
							if(is.na(choice)){ choice <- "impossible"}
							if(choice < 0) { choice <- "impossible"}
							if(choice > NROW(entry)) { choice <- "impossible"}
							inval <- inval+1
							if(inval==10) choice = 0
						} 
						if(choice==0){
							stop(paste("No value found related to", theUvalue, 'in table', table))
						} 
						entry <- entry[choice,]

					}
					#store value if we are doing ok
					if(NROW(entry)==1){ 
						if(table %in% CVtables()){
							#necessary because postgres does not support capital letters
							coln <- which(names(entry) %in% c("term","Term"))
							uvalueID[i]  <- entry[,coln]
						} else {
							#necessary because postgres does not support capital letters
							coln <- which(names(entry) %in% c("ID","id"))
							uvalueID[i]  <- entry[,coln]
						}
						break
					}
				}
			}
		}
		if(NROW(entry)==0){
			synonymID <- IgetSynonymID(getOption("odm.handler"), table=table, phrase=phraseA)
			if(length(synonymID)==1){
				uvalueID[i] <- synonymID
			} else {
				all.table <- NULL # to avoid warnings during check
				command <- paste('all.table <- Iget',table,'(options("odm.handler")[[1]])', sep='')
				eval(parse(text=command))
				if(NROW(all.table)==0){
					if(uvalue[i]=='No' | uvalue[i]=="Unknown"){
						uvalueID[i] <- IgetNo(options("odm.handler")[[1]], table)
						next
					} else {
						stop("Table ", table, " has no entries. Please enter a record for '", value, "'")
					}
				}

				if(!all(tolower(lookup[[table]]) %in% tolower(names(all.table)))){
					fields <- paste(lookup[[table]][!lookup[[table]] %in% names(all.table)], collapse="; ")
					existing <- paste(names(all.table), collapse="; ")
					stop(paste("getID error: lookup not defined correctly, missing fields for table ", table, ":", fields, ". Existing fields are:", existing ))
				}


				words <- uvalue[i]
				for(splitstr in c("(")){
					the.split0 <- strsplit(uvalue[i], splitstr, fixed=TRUE)
					if(length(the.split0[[1]])>1){
						words <- c(words, the.split0[[1]][1])
					}
				}
				the.split <- c(strsplit(uvalue[i], "-"), strsplit(uvalue[i], ","), strsplit(uvalue[i], "'"))
				# go through all splitting possibilities
				for(mutate in the.split[sapply(the.split, length) > 1]){
					perm <- permutations(length(mutate))
					for(n in 2:NROW(perm)){
						words <- c(words, paste(mutate[perm[n,]], collapse=" "))
					}
				}
				#calculate levenshtein distance for all possible combinations
				the.dist <- c()
				for(word in words){
					for(field in lookup[[table]]){
						column <- match(tolower(field), tolower(names(all.table)))

						lev.dist <- levenshtein.distance(xsource=tolower(word), targets=tolower(all.table[,column]))
						names(lev.dist) <- all.table[[field]]
						the.dist <- c(the.dist, lev.dist)
					}
				}
				
				#return best matches
				if(length(unique(names(sort(the.dist))))<10){
					possible.answer <- unique(names(sort(the.dist)))
				} else {
					m <- 10
					while(length(unique(names(sort(the.dist)[1:m])))<10){
						m <- m+1
					}
					possible.answer <- unique(names(sort(the.dist)[1:m]))[1:10]
				}
				cat("\nNo match found for ",uvalue[i],". Fuzzy search in table ",table, " returns:\n")
				print(possible.answer)
				cat("\n\n Please select matching record or hit 0 to stop or -1 to enter a search phrase\n")
				choice <- "impossible"
				if(!interactive()) stop("Error: no unique value in non-interactive session")
				while(choice == "impossible"){
					next.step <- readline("What is your choice? ")
					choice <- as.numeric(next.step)
					if(is.na(choice)){ choice <- "impossible"}
					if(choice < -1) { choice <- "impossible"}
					if(choice > length(possible.answer)) { choice <- "impossible"}
				} 
				if(choice==0){
					#if(allowNoValue){
					#	entry <- list(ID=IgetNo(options('odm.handler')[[1]], table))
					#} else {
						stop(paste("No value found related to", uvalue[i], 'in table', table))
					#}
				} else if(choice==-1){
					search.term <- readline("Please enter serach term: ")
					uvalueID[i]  <- addSynonym(table=table, uvalue[i], search.term)
				} else {
					uvalueID[i]  <- getID(table, possible.answer[choice], remove.special.character=FALSE)
					addSynonym(table, uvalue[i], uvalueID[i])
					warning(paste("Stored", uvalue[i], "as synonym of",  possible.answer[choice]))
				}
			}
		}
		
	}
	lookup <- data.frame(value=uvalue, ID=uvalueID, stringsAsFactors=FALSE)
	lookedUp <- merge(data.frame(value=value, rownr=1:length(value), stringsAsFactors=FALSE),lookup, sort=FALSE)
	lookedUp$value[is.na(lookedUp$value)] <- "NA"
	value[is.na(value)] <- "NA"
	lu.ordered <- lookedUp[order(lookedUp$rownr),]
	stopifnot(lu.ordered$value==value)
	return(lu.ordered$ID)
}
