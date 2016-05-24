# package created using:
# http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
#' @import rJava arules
#' @importFrom utils write.table

.onLoad <- function(libname, pkgname ){
	.jinit()
	if(J("java.lang.System")$getProperty("java.version") < "1.8.0") {
		stop("rCBA requires Java >= 1.8 ", call. = FALSE)
	} 
}

init <- function(){
	# initialize rJava
	.jinit()
	if(J("java.lang.System")$getProperty("java.version") < "1.8.0") {
		stop("rCBA requires Java >= 1.8 ", call. = FALSE)
	}
	# add java implementation to classpath
	.jaddClassPath(dir(paste(path.package("rCBA"), "/java/", sep=""), full.names=TRUE))
	# add jar archives to classpath
	jars <- list.files(paste(path.package("rCBA"), "/java/lib/", sep=""))
	for(jar in jars){
		if(!grepl("validation-api",jar)){
			.jaddClassPath(paste(path.package("rCBA"), "/java/lib/", jar, sep=""))
		}
	}
	# add configuration files to classpath
	confs <- list.files(paste(path.package("rCBA"), "/java/conf/", sep=""))
	for(conf in confs){
		.jaddClassPath(paste(path.package("rCBA"), "/java/conf/", conf, sep=""))
	}
}

#' A Pruning function
#'
#' @param train data.frame with training data
#' @param rules data.frame with rules
#' @param method pruning method m2cba(default)|m1cba|dcbrcba
#' @return data.frame with pruned rules
#' @export
#' @examples
#' library("arules")
#' library("rCBA")
#' data("iris")
#' 
#' train <- sapply(iris,as.factor)
#' train <- data.frame(train, check.names=FALSE)
#' txns <- as(train,"transactions")
#' 
#' rules = apriori(txns, parameter=list(support=0.03, confidence=0.03, minlen=2), 
#'	appearance = list(rhs=c("Species=setosa", "Species=versicolor", "Species=virginica"),default="lhs"))
#' rulesFrame <- as(rules,"data.frame")
#' 
#' print(nrow(rulesFrame))
#' prunedRulesFrame <- pruning(train, rulesFrame, method="m2cba")
#' print(nrow(prunedRulesFrame))
pruning <- function(train, rules, method="m2cba"){
	# init java
	init()
	print(paste(Sys.time()," rCBA: initialized",sep=""))
	# init interface
	jPruning <- J("cz.jkuchar.rcba.r.RSpring")$initializePruning()
	# set column names
	.jcall(jPruning, , "setColumns", .jarray(colnames(train)))
	# add train items
	
	# sequence
	# trainConverted <- data.frame(lapply(train, as.character), stringsAsFactors=FALSE)
	# for(i in 1:nrow(trainConverted)){
	# 	.jcall(jPruning, , "addItem", as.character(unname(unlist(trainConverted[i,]))))
	# }

	# all
	# trainConverted <- data.frame(lapply(train, as.character), stringsAsFactors=FALSE)
	# .jcall(jPruning, , "addAll", as.vector(t(trainConverted)))
	
	# file
	f<-tempfile()
	write.table(train, file=f, sep=",", row.names=FALSE, col.names=FALSE)
	.jcall(jPruning, , "loadItemsFromFile", as.character(f))
	print(paste(Sys.time()," rCBA: dataframe ",nrow(train),"x",ncol(train),sep=""))
	# add rules
	for (i in 1:nrow(rules)){
		.jcall(jPruning, , "addRule", as.character(rules[i,]$rules), as.numeric(rules[i,]$confidence), as.numeric(rules[i,]$support), as.numeric(rules[i,]$lift))
	}
	print(paste(Sys.time()," rCBA: rules ",nrow(rules),"x",ncol(rules),sep=""))
	# perform pruning
	jPruned <- .jcall(jPruning, "[Lcz/jkuchar/rcba/rules/Rule;", "prune", as.character(method))
	print(paste(Sys.time()," rCBA: pruning completed",sep=""))
	# build dataframe
	pruned <- data.frame(rules=rep("", 0), support=rep(0.0, 0), confidence=rep(0.0, 0), lift=rep(0.0, 0), stringsAsFactors=FALSE) 
	for(jRule in jPruned){
		pruned[nrow(pruned) + 1,] <- c(.jcall(jRule, "S", "getText"), .jcall(jRule, "D", "getSupport"), .jcall(jRule, "D", "getConfidence"), .jcall(jRule, "D", "getLift"))
	}
	print(paste(Sys.time()," rCBA: pruned rules ",nrow(pruned),"x",ncol(pruned),sep=""))
	pruned$support <- as.double(pruned$support)
	pruned$confidence <- as.double(pruned$confidence)
	pruned$lift <- as.double(pruned$lift)
	pruned
}