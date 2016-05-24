setClass("BiclustSet",
		representation(GenesMembership="matrix",ColumnMembership="matrix",
				Number="numeric")
)
#constructor for the class BiclustSet
BiclustSet <-function(x,fabia.thresZ=0.5,fabia.thresL=NULL){myBiclustSet = new("BiclustSet",GenesMembership=x@RowxNumber, ColumnMembership=x@NumberxCol, 
			Number=x@Number)}
		
setGeneric("BiclustSet")
setMethod("BiclustSet",signature(x="Biclust"),
		function(x){myBiclustSet = new("BiclustSet",GenesMembership=x@RowxNumber, ColumnMembership=x@NumberxCol, 
				Number=x@Number)})
setMethod("BiclustSet",signature("Factorization"),function(x,fabia.thresZ=0.5,fabia.thresL=NULL){
	#require(fabia)
	N <-  x@p1
	tmp <- extractBic(x,thresZ=fabia.thresZ,thresL=fabia.thresL)
	X <- x@X
	FabiaBicl <- tmp$numn
	DataFabiaColB <- c()
	DataFabiaRowB <- c()
	for (j in 1:N){
		colMem <- vector(length= ncol(X), mode="logical")
		rowMem <- vector (length = nrow(X), mode = "logical")
		
		if(length(unlist(FabiaBicl[j,2])) > 0){
			for (i in unlist(FabiaBicl[j,2])) {colMem[i] <- TRUE} 
			for (i in unlist(FabiaBicl[j,1])) {rowMem[i] <- TRUE}
			DataFabiaColB <- rbind(DataFabiaColB, colMem)
			DataFabiaRowB <- cbind(DataFabiaRowB,  rowMem)
		}else next
	}
	myBiclustSet = new("BiclustSet",GenesMembership=DataFabiaRowB , ColumnMembership=DataFabiaColB, 
			Number=N)
} )

setMethod("BiclustSet",signature("list"),
		function(x){
			biclust1Number<- ncol(x$rows)
			bicArows <- matrix(as.logical(x$rows),ncol=biclust1Number,byrow=T)
			bicAcols <- matrix(as.logical(x$columns),nrow=biclust1Number)
			myBiclustSet = new("BiclustSet",GenesMembership=bicArows, ColumnMembership=bicAcols, 
					Number=biclust1Number)
		})
###**show and summary the same as for Biclust class******************************

setMethod("show", "BiclustSet",
		function(object)
		{
			cat("\nAn object of class",class(object),"\n\n")			
			n<-object@Number
			n<-min(c(n,5))
			if(n>1)
			{
				cat("\nNumber of Clusters found: ",object@Number, "\n")
				cat("\nFirst ",n," Cluster sizes:\n")
				
				rowcolsizes<-rbind(colSums(object@GenesMembership[,1:n]),rowSums(object@ColumnMembership[1:n,]))
				rownames(rowcolsizes)<-c("Number of Rows:","Number of Columns:")
				colnames(rowcolsizes)<-paste("BC", 1:n)
				#print.default(format(rowcolsizes, print.gap = 2, quote = FALSE))
				print(rowcolsizes)
			}
			else
			{
				if(n==1) cat("\nThere was one cluster found with\n ",sum(object@GenesMembership[,1]), "Rows and ", sum(object@ColumnMembership), "columns")
				if(n==0) cat("\nThere was no cluster found")
			}
			cat("\n\n")
		})

setGeneric("summary")
setMethod("summary", "BiclustSet",
		function(object)
		{
			cat("\nAn object of class",class(object),"\n\n")
			n<-object@Number
			
			if(n>1)
			{
				cat("\nNumber of Clusters found: ",object@Number, "\n")
				cat("\nCluster sizes:\n")
				
				rowcolsizes<-rbind(colSums(object@GenesMembership[,1:n]),rowSums(object@ColumnMembership[1:n,]))
				rownames(rowcolsizes)<-c("Number of Rows:","Number of Columns:")
				colnames(rowcolsizes)<-paste("BC", 1:n)
				#print.default(format(rowcolsizes, print.gap = 2, quote = FALSE))
				print(rowcolsizes)
			}
			else
			{
				if(n==1) cat("\nThere was one cluster found with\n ",sum(object@GenesMembership[,1]), "Rows and ", sum(object@ColumnMembership), "columns")
				if(n==0) cat("\nThere was no cluster found")
			}
			cat("\n\n")
			
			
		})

