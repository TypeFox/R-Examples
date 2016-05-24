apexFeatures <- function(x, ...) UseMethod("apexFeatures")

apexFeatures.default <- function(x, ...) {
	if (class(x)!="data.frame" || dim(x)[1] < 1 || dim(x)[2] != 2 || sort(names(x)) != c("apex","peptide_sequence")) {stop("The input for apexFeatures is a mandatory data frame containing the variables in the model. The data frame requires the columns \"peptide_sequence\", \"apex\". The data may contain training data (with boolean \"apex\" and test data (with \"apex\"=NA))", call. = FALSE)}

	x <- unique(x)
	
	object<-as.data.frame(rbindlist(apply(x,1,pcfeatures.apexFeatures)))

	object$apex<-as.factor(object$apex)
	
	class(object) <- "apexFeatures"
	
	return(object)
}

pcfeatures.apexFeatures <- function(x, ...) {
	aaspecies <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

	featurelist <- c("FASG760101","CHOP780201","CHOP780202","CHOP780203","WERD780101","ZIMJ680104","KLEP840101","EISD860102","FAUJ880111","VINM940101","FAUJ880103","GUYH850105","NOZY710101")
	
	features <- data.frame(peptide_sequence=x[1][[1]], apex=x[2][[1]], stringsAsFactors=FALSE)

	aasequence <- sapply(strsplit(features[["peptide_sequence"]],"")[[1]],function(X){if(X %in% aaspecies){return(X)}else{return(NA)}})
	aasequence <- aasequence[!is.na(aasequence)]

	features[["length"]] <- length(aasequence)

	a <- 1
	while (a <= length(aaspecies)) {
		features[[paste(aaspecies[a],"_sum",sep="")]] <- 0
		a <- a + 1
	}
	
	b <- 1
	while (b <= features[["length"]]) {
		features[[paste(aasequence[b],"_sum",sep="")]] <- features[[paste(aasequence[b],"_sum",sep="")]] + 1
		b <- b + 1
	}
	
	c <- 1
	while (c <= length(aaspecies)) {
		features[[paste(aaspecies[c],"_avg",sep="")]] <- features[[paste(aaspecies[c],"_sum",sep="")]] / features[["length"]]
		c <- c + 1
	}
	
	i <- 1
	while (i <= length(featurelist)) {
		feature_sum <- 0
	
		j <- 1
		while (j <= features[["length"]]) {
			if (!is.na(bio3d::aa.index[[featurelist[i]]]$I[[aasequence[j]]])) {
				feature_sum <- feature_sum + bio3d::aa.index[[featurelist[i]]]$I[[aasequence[j]]]
			}
			j <- j + 1
		}
		features[[paste(featurelist[i],"_sum",sep="")]] <- feature_sum
		features[[paste(featurelist[i],"_avg",sep="")]] <- (feature_sum / features[["length"]])
		i <- i + 1
	}

	row.names(features) <- NULL
	return(features)
}

print.apexFeatures <- function(x, ...) {
    if (!inherits(x, "apexFeatures")) stop("Is not a apexFeatures object")		
	class(x) = "data.frame"
	cat("apexFeatures\n")
	cat("Features: ")
	cat(dim(x)[2])
	cat("\n")
	cat("Dataset size: ")
	cat(dim(x)[1])
	cat("\n")
}

