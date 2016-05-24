read.pheno <-
function(file.format = c("cape", "csv"), filename = NULL, pheno.col = NULL, id.col = 1, delim = ",", na.strings = "-") {

	if(length(file.format) > 1){
		stop("A file format must be specified.")
		}

	if(file.format == "cape"){
		if(is.null(filename)){
			filename <- file.choose()
			}
			
		pop.data <- read.table(filename, na.strings = na.strings, stringsAsFactors = FALSE, sep = delim, header = TRUE)
	
		
		ids <- pop.data[,id.col]
	
		beginGeno = match(FALSE, is.na(suppressWarnings(as.numeric(pop.data[1,]))))
		endPheno <- beginGeno-1
	
		not.na.rows <- which(apply(pop.data[,1:endPheno], 1, function(x) !all(is.na(x))))
		
		#add a check for non-numeric phenotypes
		pheno.classes <- NULL
		for(i in 2:endPheno){
			pheno.classes <- c(pheno.classes, class(pop.data[,i]))
			}
			
		char.pheno <- which(pheno.classes == "character")
	
		if(length(char.pheno) > 0){
			message("All phenotypes must be numeric.\nMake sure NA strings are specified correctly.")
			cat("The following phenotype columns have character values:", colnames(pop.data)[char.pheno], sep = "\n")
			return(NULL)
			}
	
	
		#if no phenotypes are specified, just take all phenotypes
		if(is.null(pheno.col)){
			pheno.col <- 1:endPheno
			pheno.col <- pheno.col[-id.col]
			}
			
		#if phenotypes are specified as characters, find their
		#locations
		pheno.columns <- get.col.num(pop.data, pheno.col)
			
		#take out the phenotype matrix
		#It begins in the third row and includes the columns
		#specified by the user
		pheno <- as.matrix(pop.data[not.na.rows,pheno.columns])
	
		#convert pheno into a numeric matrix, so we can do 
		#matrix algebra on it later. If there are any non-numeric
		#phenotypes, convert them to numeric 
		pheno <- matrix(apply(pheno, 2, as.numeric), ncol = dim(pheno)[2], byrow = FALSE)
		colnames(pheno) <- colnames(pop.data)[pheno.columns]
		rownames(pheno) <- ids[not.na.rows]
		
	    return(pheno)
		}	    
    

	if(file.format == "csv"){
		if(is.null(filename)){
			filename <- file.choose()
			}
			
		pop.data <- read.table(filename, na.strings = na.strings, stringsAsFactors = FALSE, sep = delim, header = TRUE)
	
		
		ids <- pop.data[,id.col]
		
		not.na.rows <- which(apply(pop.data, 1, function(x) !all(is.na(x))))
		
		#add a check for non-numeric phenotypes
		pheno.classes <- NULL
		for(i in 1:dim(pop.data)[2]){
			pheno.classes <- c(pheno.classes, class(pop.data[,i]))
			}
			
		char.pheno <- which(pheno.classes == "character")
	
		if(length(char.pheno) > 0){
			message("All phenotypes must be numeric.\nMake sure NA strings are specified correctly.")
			cat("The following phenotype columns have character values:", colnames(pop.data)[char.pheno], sep = "\n")
			return(NULL)
			}
	
	
		#if no phenotypes are specified, just take all phenotypes
		if(is.null(pheno.col)){
			pheno.col <- 1:dim(pop.data)[2]
			pheno.col <- pheno.col[-id.col]
			}
			
		#if phenotypes are specified as characters, find their
		#locations
		pheno.columns <- get.col.num(pop.data, pheno.col)
			
		#take out the phenotype matrix
		#It begins in the third row and includes the columns
		#specified by the user
		pheno <- as.matrix(pop.data[not.na.rows,pheno.columns])
	
		#convert pheno into a numeric matrix, so we can do 
		#matrix algebra on it later. If there are any non-numeric
		#phenotypes, convert them to numeric 
		pheno <- matrix(apply(pheno, 2, as.numeric), ncol = dim(pheno)[2], byrow = FALSE)
		colnames(pheno) <- colnames(pop.data)[pheno.columns]
		rownames(pheno) <- ids[not.na.rows]
		
	    return(pheno)
		}	    


}
