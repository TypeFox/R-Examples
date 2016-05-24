########### optCluster Methods ###########

########### Create Accessor Functions ##########
         
## clValid object accessor
setGeneric("getDataset", function(object, ...) standardGeneric("getDataset"))
setMethod("getDataset",signature(object="optCluster"),
          function(object) return(object@inputData))
          
setGeneric("getClValid", function(object, ...) standardGeneric("getClValid"))
setMethod("getClValid",signature(object="optCluster"),
          function(object) return(object@clVal))
          
## cluster algorithm rankings accessor
setGeneric("methodRanks", function(object, ...) standardGeneric("methodRanks"))
setMethod("methodRanks",signature(object="optCluster"),
          function(object) return(object@ranksWeights$ranks))

## validation score rankings accessor
setGeneric("scoreRanks", function(object, ...) standardGeneric("scoreRanks"))
setMethod("scoreRanks",signature(object="optCluster"),
          function(object) return(object@ranksWeights$weights))
          
## raggr object accessor
setGeneric("getRankAggreg", function(object, ...) standardGeneric("getRankAggreg"))
setMethod("getRankAggreg",signature(object="optCluster"),
          function(object) return(object@rankAgg))

## top method accessor
setGeneric("topMethod", function(object, ...) standardGeneric("topMethod"))
setMethod("topMethod",signature(object="optCluster"),
          function(object) return(object@rankAgg$top.list[1]))

## validation measure names accessor
setGeneric("measureNames", function(object, ...) standardGeneric("measureNames"))
setMethod("measureNames",signature(object="optCluster"),
          function(object) return(getClValid(object)@measNames))

## clustering algorithm name accessor
setGeneric("methodNames", function(object, ...) standardGeneric("methodNames"))
setMethod("methodNames",signature(object="optCluster"),
          function(object) return(getClValid(object)@clMethods))           
            
########## Additional Methods ##########
         
## clustering results 
setGeneric("clusterResults", function(object, ...) standardGeneric("clusterResults"))
setMethod("clusterResults",signature(object="optCluster"),
          function(object,method=methodNames(object)) {
            method <- match.arg(method,methodNames(object))
            clValObj <- getClValid(object)
            return(clValObj@clusterObjs[[method]])})

## validation measures 
setGeneric("valScores", function(object, ...) standardGeneric("valScores"))
setMethod("valScores",signature(object="optCluster"),
          function(object,measures=measureNames(object)) {
            measures <- match.arg(measures,measureNames(object),several.ok=TRUE)
            clValObj <- getClValid(object)         
            return(clValObj@measures[measures,,,drop=FALSE])})          

## optimal validation scores          
setGeneric("optimalScores",function(object, ...) standardGeneric("optimalScores"))
setMethod("optimalScores","optCluster",
          function(object,digits = max(3,getOption("digits")-3)) {
            ## Find best scores
            best <- scoreRanks(object)[,1]
			bestMethods <- methodRanks(object)[,1]
			methodNums <-as.numeric(gsub("\\D", "", bestMethods))
			methodNames <- as.character(gsub("\\d", "", bestMethods))
			methodNames <- gsub("-", "", methodNames)            
            cat("Optimal Scores:\n\n") 
            return(data.frame("Score"=round(best,digits),"Method"=methodNames,"Clusters"=methodNums))
          })
          
########### Print, Show, and Summary Methods ##########

setMethod("print","optCluster",
          function(x) {
            cat("\nThe overall optimal clustering method and number of clusters is: \n\t   ", topMethod(x), "\n\n")
            print(getRankAggreg(x))
            cat("  Iterations: ", getRankAggreg(x)$num.iter, "\n")
          })

setMethod("show","optCluster",
          function(object) {
            cat("\nThe overall optimal clustering method and number of clusters is: \n\t   ", topMethod(object), "\n\n")
            print(getRankAggreg(object))
            cat("  Iterations: ", getRankAggreg(object)$num.iter, "\n")
            })

setMethod("summary","optCluster",
          function(object,digits = max(3,getOption("digits")-3)) {
            cat("\nClustering Methods:\n",methodNames(object),"\n\n")
            cat("Cluster sizes:\n",nClusters(getClValid(object)),"\n\n")
            cat("Validation Measures:\n")
            print(ftable(round(valScores(object),digits),row.vars=c(3,1)))
            cat("\n")
            print(optimalScores(object))            
            cat("\nThe overall optimal clustering method and number of clusters is: \n\t   ", topMethod(object), "\n\n")
            print(getRankAggreg(object))
            cat("  Iterations: ", getRankAggreg(object)$num.iter, "\n")
            })
