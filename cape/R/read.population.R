read.population <-
function(filename = NULL, pheno.col = NULL, geno.col = NULL, delim = ",", na.strings = "-", check.chr.order = TRUE) {

		if(is.null(filename)){
			filename <- file.choose()
			}
			
		cross.data <- read.table(filename, na.strings = na.strings, stringsAsFactors = FALSE, sep = delim, header = TRUE)
	
	
	    #determine where phenotypes end and genotypes begin by blanks in row 1
		beginGeno = match(FALSE, is.na(suppressWarnings(as.numeric(cross.data[1,]))))
		# beginGeno = min(which(!is.na(cross.data[1,])))
	
		#add a check for non-numeric phenotypes
		pheno.classes <- NULL
		for(i in 1:(beginGeno-1)){
			pheno.classes <- c(pheno.classes, class(cross.data[,i]))
			}
			
		char.pheno <- which(pheno.classes == "character")
	
		if(length(char.pheno) > 0){
			message("All phenotypes must be numeric.")
			cat("The following phenotype columns have character values:", colnames(cross.data)[char.pheno], sep = "\n")
			cat("This error can occur if NA's are coded with multiple characters. Make sure NA coding is consistent throughout the data set.")
			return(NULL)
			}
	
		chr <- as.vector(as.matrix(cross.data[1,beginGeno:dim(cross.data)[2]]))
		
		if(check.chr.order){
			x.locale <- grep("x", chr, ignore.case = TRUE)
			y.locale <- grep("y", chr, ignore.case = TRUE)
			just.num.chr <- setdiff(1:length(chr), c(x.locale, y.locale))
			consec.chr <- consec.pairs(as.numeric(chr[just.num.chr]))
			order.check <- apply(consec.chr, 1, function(x) x[2] - x[1])
			if(length(which(order.check < 0)) > 0){
				warning("The chromosomes appear to be out of order.\nIt is best to sort the chromosomes before beginning the analysis.")
				}
			}
			
		marker.loc <- as.numeric(cross.data[2,beginGeno:dim(cross.data)[2]])
	
		#take out the genotype matrix
		#It begins in the 3rd row after chromosome numbers and marker locations
		geno <- as.matrix(cross.data[3:dim(cross.data)[1],beginGeno:dim(cross.data)[2]])
	
		#if no phenotypes are specified, just take all phenotypes
		if(is.null(pheno.col)){
			pheno.col <- 1:(beginGeno-1)
			}
			
		#if phenotypes are specified as characters, find their
		#locations
		pheno.columns <- get.col.num(cross.data, pheno.col)
			
		#take out the phenotype matrix
		#It begins in the third row and includes the columns
		#specified by the user
		pheno <- as.matrix(cross.data[3:dim(cross.data)[1],pheno.columns])
	
		#convert pheno into a numeric matrix, so we can do 
		#matrix algebra on it later. If there are any non-numeric
		#phenotypes, convert them to numeric 
		pheno <- matrix(apply(pheno, 2, as.numeric), ncol = dim(pheno)[2], byrow = FALSE)
		colnames(pheno) <- colnames(cross.data)[pheno.columns]
		rownames(pheno) <- 1:dim(pheno)[1]
	   
	   if(is.null(geno.col)){
		   	geno.col <- 1:dim(geno)[2]
	   		}
	   	geno.columns <- get.col.num(geno, geno.col)
	  	geno <- geno[,geno.columns] 
	  	chr <- chr[geno.columns]
		marker.loc <- marker.loc[geno.columns]
	   
	    #run a check to see how the genotypes are stored.
		#genotypes can be stored as (0,1,2), ("A","H","B")
		#or as probabilities between 0 and 1
		#if the genotypes are stored as (0,1,2) or 
		#letters, we need to convert them to probabilities
	    
	    found.genotype.mode <- 0 #create a flag to determine whether we have figured out how the genotypes are encoded
	    
	    genotype.class <- class(cross.data[,beginGeno])
	    if(genotype.class != "character"){
		    all.genotypes <- sort(unique(na.omit(as.numeric(as.matrix(geno))))) #get a vector of all genotypes used
		    }else{
			    all.genotypes <- sort(unique(na.omit(as.vector(as.matrix(geno)))))
			    if(length(all.genotypes) > 3){
			    	#look for empty genotypes
			    	cat("I have detected", length(all.genotypes), "genotypes:\n")
			    	print(all.genotypes)
					stop("Please check for missing genotype values or other errors in the genotype data.")
			    	}
		    	}
		
		#check to see if the genotypes are encoded as letters
	
	  	if(genotype.class == "character"){
	  		het.present <- grep("H", all.genotypes)
	  		if(length(all.genotypes) > 2 && length(het.present) == 0){
	  			stop("Heterozygotes must be coded by H")
	  			}
	  			
	 		
	  		found.genotype.mode <- 1 #indicate that we have found the genotype mode for the file
			#assign baseGeno and notBaseGeno
		    #the baseGeno is assigned a numeric value
		    #of 0. The notBaseGeno is assigned 1
		    #by default we make the first letter
		    #alphabetically the base genotype
		    if(length(all.genotypes) == 3){
				baseGeno <- all.genotypes[all.genotypes != "H"][1]
				notBaseGeno <- all.genotypes[all.genotypes != "H"][2]
				cat("The genotypes are encoded as ", baseGeno, ", H, ", notBaseGeno, "\nConverting to 0, 0.5, 1.\n", sep = "")
				}else{
				# baseGeno <- all.genotypes[all.genotypes != "H"][1]
				# notBaseGeno <- "H"
				baseGeno <- sort(all.genotypes)[1]
				notBaseGeno <- sort(all.genotypes)[2]
				cat("The genotypes are encoded as ", baseGeno, " and ", notBaseGeno, "\nConverting to 0 and 1.\n", sep = "")
		 		}
		
	
		
		    #turn baseGeno, H, and notBaseGeno to 0, .5, and 1 respectively
			#This function takes in a vector and converts the letters to
			#the appropriate numbers
		
			convert.geno.letter <- function(genotypes){
				genotypes[which(as.character(genotypes) == baseGeno)] <- 0
		    	genotypes[which(as.character(genotypes) == notBaseGeno)] <- 1
		    	if(length(all.genotypes) == 3){
		    		genotypes[which(as.character(genotypes) == "H")] <- 0.5
		    		}else{
		    		genotypes[which(as.character(genotypes) == "H")] <- 1
		    		}
		    	return(as.numeric(genotypes))
		    	}
		 	geno <- apply(geno, 2, convert.geno.letter) 
	  		}
	    
		
	 	#check to see if the genotypes are encoded as (0, 1, 2)
	 	numeric.test <- which(all.genotypes == 2) #check for 2, since 2 is unique to this encoding
	 	if(length(numeric.test) > 0){
	 		outside.upper.bound <- which(all.genotypes > 2)
	 		outside.lower.bound <- which(all.genotypes < 0)
	 		if(length(outside.upper.bound) > 0 || length(outside.lower.bound) > 0){
	 			stop("Assuming (0,1,2) coding, but I detected genotypes greater than 2 or less than 0.")
	 			}
	 		cat("The genotypes are encoded as 0, 1, 2.\nConverting to 0, 0.5, 1.\n")
	 		found.genotype.mode <- 1 #set the flag indicating we've figured out the encoding
	 		#turn 0, 1, 2 into 0, 0.5 and 1 respectively
			convert.geno.number <- function(genotypes){
		        genotypes[which(as.numeric(genotypes) == 1)] <- 0.5
		        genotypes[which(as.numeric(genotypes) == 2)] <- 1
		        return(as.numeric(genotypes))
				}
		
		 	geno <- apply(geno, 2, convert.geno.number) 
	
		 	}
	 	
	 	#if we still haven't found the genotype mode yet
	 	#check to see if the genotypes are encoded as probabilities
	 	if(found.genotype.mode == 0){
	 		min.geno <- min(all.genotypes); max.geno <- max(all.genotypes) #find the max and min values for genotypes
	 		if(min.geno >= 0 && max.geno <= 1){ #and make sure they are bounded as probabilities are
	 			found.genotype.mode <- 1 #set the flag to indicate we have found the genotype mode
	 			cat("The genotypes are encoded as probabilities.\nNo conversion needed.\n")
	 			
				#we still need to conver the data frame to a numeric matrix for later compatibility
	 			convert.geno.prob <- function(genotypes){
			        return(as.numeric(genotypes))
					}
		
			 	geno <- apply(geno, 2, convert.geno.prob) 
	
	 			}
	 		}
	 	
	 	
	 	#If after all this, we haven't found the genotype encoding, stop and warn the user
	 	if(found.genotype.mode == 0){
	 		stop("\nGenotypes must be encoded as (0, 1, 2), (A,H,B), or probabilities.\n")
	 		}
	
	
	
		#take out the sex chromosomes and invariant markers
		x.locale <- grep("X", chr, ignore.case = TRUE)
		if(length(x.locale) > 0){
			message("\nRemoving markers on the X chromosome")
			geno <- geno[,-x.locale]
			chr <- chr[-x.locale]
			marker.loc <- marker.loc[-x.locale]
			}
			
		y.locale <- grep("Y", chr, ignore.case = TRUE)
		if(length(y.locale) > 0){
			message("\nRemoving markers on the Y chromosome")
			geno <- geno[,-y.locale]
			chr <- chr[-y.locale]
			marker.loc <- marker.loc[-y.locale]
			}
			
		#take out markers with only 1 allele
		num.allele <- apply(geno, 2, function(x) length(unique(x)))
		mono.allele <- which(num.allele == 1)
		if(length(mono.allele) > 0){
			message("\nRemoving invariant markers.\n")
			geno <- geno[,-mono.allele]
			chr <- chr[-mono.allele]
			marker.loc <- marker.loc[-mono.allele]
			}
		
	
	
		na.locale <- which(is.na(geno))
		if(length(na.locale) > 0){
			message("I have detected missing values in the genotype matrix.\n\tIf you are planning to use the kinship correction, please use impute.geno() to impute the genotype data.\n")
			}
	
		#put in code here to distribute the genotypes between -1 and 1 so we get symmetric m12/m21 null distributions
		#construct the data object
		marker.names <- colnames(geno)
		colnames(geno) <- 1:dim(geno)[2]
		rownames(geno) <- 1:dim(geno)[1]
	
		#scale the genotypes to lie between 0 and 1
		#even if it's a backcross
		geno <- geno/max(geno, na.rm = TRUE)
		marker.num <- 1:dim(geno)[2]
	
		final.data <- list(pheno, geno, chr, marker.names, marker.num, marker.loc)
		names(final.data) <- c("pheno", "geno", "chromosome", "marker.names", "marker.num", "marker.location")
		     
		cat("Read in the following data:\n")
		cat("\t-", dim(pheno)[1], "individuals -\n")
		cat("\t-", dim(geno)[2], "markers -\n")
		cat("\t-", dim(pheno)[2], "phenotypes -\n")
	
	
	    return(final.data) 
	   
	    
}
