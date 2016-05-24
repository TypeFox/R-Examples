# summary.r: various functions to print summarized information on an object

summary.RecLinkData <- function(object,...)
{
    if (!("RecLinkData" %in% class(object)) &&
        !("RecLinkResult" %in% class(object)))
        stop(sprintf("Wrong type for object: %s!", class(object)))
    if (object$type=="linkage")
    {
        cat("\nLinkage Data Set\n\n")
        cat(sprintf("%d records in data set 1",nrow(object$data1)),"\n")
        cat(sprintf("%d records in data set 2",nrow(object$data2)),"\n")
    }
    else
    {
        cat("\nDeduplication Data Set\n\n")
        cat(sprintf("%d records",nrow(object$data)),"\n")
    }

    cat(sprintf("%d record pairs",nrow(object$pairs)),"\n")
    cat("\n")
    # the expression "length(which(..." is needed to eliminate NAs
	cat(sprintf("%d matches\n",
        length(which(object$pairs$is_match==TRUE))))
    cat(sprintf("%d non-matches\n",
        length(which(object$pairs$is_match==FALSE))))
	cat(sprintf("%d pairs with unknown status\n",
        length(which(is.na(object$pairs$is_match)))))
    cat("\n")

	if (!is.null(object$Wdata))
	{
		cat("\n")
		cat("Weight distribution:\n\n")
		h=hist(object$Wdata,plot=FALSE)
		c=h$counts
		# nehme Gewichtsintervalle als Indizes, um Histogrammansicht zu erhalten
    names(c)=sapply(1:(length(h$breaks)-1),
      function(x) sprintf("(%g,%g]",h$breaks[x],h$breaks[x+1]))
    # erstes Intervall ist auch links geschlossen
    names(c)[1]=sprintf("[%g,%g]", h$breaks[1], h$breaks[2])
		print(c)
	}
	return(invisible(NULL))
}



summary.RecLinkResult <- function (object, ...)
{
    if (!("RecLinkResult" %in% class(object)))
        stop(sprintf("Wrong type for object: %s!", class(object)))

    summary.RecLinkData(object,...)
    crossTable <- table(as.logical(object$pairs$is_match),object$prediction,
          dnn=list("true status","classification"),useNA="ifany")

    cat("\n")
    cat(sprintf("%d links detected", sum(crossTable[,"L"])),"\n")
    cat(sprintf("%d possible links detected", sum(crossTable[,"P"])),"\n")
    cat(sprintf("%d non-links detected", sum(crossTable[,"N"])),"\n")

    cat("\n")

    # print the following summary only if matching status is known at leas
    # for some pairs
    if (nrow(crossTable) >= 2)
    {
      TP=crossTable["TRUE", "L"] # true positive
      FP=crossTable["FALSE", "L"] # false positive
      TN=crossTable["FALSE", "N"] # true negative
      FN=crossTable["TRUE", "N"] # false negative

      alpha=FN/(TP+FN)
      beta=FP/(TN+FP)
      accuracy=(TP+TN)/(TP+TN+FP+FN)
      cat(sprintf("alpha error: %f\n",alpha))
      cat(sprintf("beta error: %f\n",beta))
      cat(sprintf("accuracy: %f\n",accuracy))
      cat("\n\n")
    }
    cat("Classification table:\n\n")
    print(crossTable)
  	return(invisible(NULL))
}

texSummary <- function (object)
{
    TP=length(which(object$pairs$is_match & object$prediction=="L")) # true positive
    FP=length(which(!object$pairs$is_match & object$prediction=="L")) # false positive
    TN=length(which(!object$pairs$is_match & object$prediction=="N")) # true negative
    FN=length(which(object$pairs$is_match & object$prediction=="N")) # false negative
    
    alpha=FN/(TP+FN)
    beta=FP/(TN+FP)
    accuracy=(TP+TN)/(TP+TN+FP+FN)

    cat("\\begin{description}\n")
    cat(sprintf("\\item[alpha error] %f\n",alpha))
    cat(sprintf("\\item[beta error] %f\n",beta))
    cat(sprintf("\\item[accuracy] %f\n",accuracy))
    cat("\\end{description}\n")
    cat("\n")
#    cat("Classification table:\n\n")
    print(xtable(table(as.logical(object$pairs$is_match),object$prediction,
          dnn=list("true status","classification"),useNA="ifany")),
          floating=FALSE, latex.environments=NULL)
}


setMethod(
  f = "show",
  signature = "RLBigData",
  definition = function(object)
  {
    # shortcut to database connection
    if (class(object)=="RLBigDataDedup")
    {
      cat("\nDeduplication Data Set for large number of data\n\n")
      cat(sprintf("%d records",nrow(object@data)),"\n")
    }

    if (class(object)=="RLBigDataLinkage")
    {
      cat("\nLinkage Data Set for large number of data\n\n")
      cat(sprintf("%d records in first data set",nrow(object@data1)),"\n")
      cat(sprintf("%d records in second data set",nrow(object@data2)),"\n")
    }
    cat(sprintf("%d record pairs\n", nrow(object@pairs)))
  }
)

setMethod(
  f = "show",
  signature = "RLResult",
  definition = function(object)
  {
    cat("\nClassification result for large data set\n\n")
    cat(sprintf("%d record pairs\n", nrow(object@pairs)))
  }
)


#setMethod(
#  f = "summary",
#  signature = "RLBigDataDedup",
#  definition = function(object)
summary.RLBigDataDedup <- function(object, ...)  {
    val <- list()
    val[["nData"]] <- nrow(object@data)
    val[["attributes"]] <- if (length(object@excludeFld) == 0)
        names(object@data) else names(object@data)[-object@excludeFld]
    val[["blockFld"]] <- lapply(object@blockFld, function(x) names(object@data)[x])
    val$nPairs <- nrow(object@pairs)
    val$nMatches <- getMatchCount(object)
    val$nNonMatches <- getNonMatchCount(object)
    val$nUnknown <- getNACount(object)

    if (hasWeights(object))
    {
  		h=hist(as.ram(object@Wdata), plot=FALSE)
  		c=h$counts
  		# nehme Gewichtsintervalle als Indizes, um Histogrammansicht zu erhalten
      names(c)=sapply(1:(length(h$breaks)-1),
        function(x) sprintf("(%g,%g]",h$breaks[x],h$breaks[x+1]))
      # erstes Intervall ist auch links geschlossen
      names(c)[1]=sprintf("[%g,%g]", h$breaks[1], h$breaks[2])
      val$weightHist <- c
    }
    class(val) <- c("summaryRLBigDataDedup", "list")
    val
  }
#)

print.summaryRLBigDataDedup <- function(x, ...)
{
  cat("RLBigDataDedup, Deduplication object\n")
  cat("\n")
  cat("Number of records:", x$nData, "\n")
  cat("Attributes:", paste(x$attributes, collapse=", "), "\n")
  cat("Blocking definition:", paste(sapply(x$blockFld,
    function(x) paste("[", paste(x, collapse=", "), "]", sep="")), collapse = ", "), "\n")
  cat("Number of record pairs:", x$nPairs, "\n")
  cat("Number of matches:", x$nMatches, "\n")
  cat("Number of non-matches:", x$nNonMatches, "\n")
  if (x$nUnknown > 0) cat("Number of pairs with unknown status:", x$nUnknown, "\n")
  if(!is.null(x$weightHist))
  {
    cat("Weight histogram:\n")
    print(x$weightHist)
  }
}


#setMethod(
#  f = "summary",
#  signature = "RLBigDataLinkage",
#  definition = function(object)
summary.RLBigDataLinkage <- function(object, ...)  {
    val <- list()
    val[["nData1"]] <- as.vector(dbGetQuery(object@con, "select count(*) from data1")$count)
    val[["nData2"]] <- as.vector(dbGetQuery(object@con, "select count(*) from data2")$count)
    val[["attributes"]] <- if (length(object@excludeFld) == 0)
        names(object@data1) else names(object@data)[-object@excludeFld]
    val[["blockFld"]] <- lapply(object@blockFld, function(x) names(object@data1)[x])
    val$expectedSize <- getExpectedSize(object)
    val$nMatches <- getMatchCount(object)
    val$nNonMatches <- getNonMatchCount(object)
    val$nUnknown <- getNACount(object)
    if (hasWeights(object))
    {
  		h=hist(as.ram(object@Wdata), plot=FALSE)
  		c=h$counts
  		# nehme Gewichtsintervalle als Indizes, um Histogrammansicht zu erhalten
      names(c)=sapply(1:(length(h$breaks)-1),
        function(x) sprintf("(%g,%g]",h$breaks[x],h$breaks[x+1]))
      # erstes Intervall ist auch links geschlossen
      names(c)[1]=sprintf("[%g,%g]", h$breaks[1], h$breaks[2])
      val$weightHist <- c
    }

    class(val) <- c("summaryRLBigDataLinkage", "list")
    val
  }
#)

print.summaryRLBigDataLinkage <- function(x, ...)
{
  cat("RLBigDataLinkage, Linkage object\n")
  cat("\n")
  cat("Number of records in dataset 1:", x$nData1, "\n")
  cat("Number of records in dataset 2:", x$nData2, "\n")
  cat("Attributes:", paste(x$attributes, collapse=", "), "\n")
  cat("Blocking definition:", paste(sapply(x$blockFld,
    function(x) paste("[", paste(x, collapse=", "), "]", sep="")), collapse = ", "), "\n")
  cat("Number of record pairs:", x$nPairs, "\n")
  cat("Number of matches:", x$nMatches, "\n")
  cat("Number of non-matches:", x$nNonMatches, "\n")
  if (x$nUnknown > 0) cat("Number of pairs with unknown status:", x$nUnknown, "\n")
  if(!is.null(x$weightHist))
  {
    cat("Weight histogram:\n")
    print(x$weightHist)
  }
}


#setMethod(
#  f = "summary",
#  signature = "RLResult",
#  definition = function(object)
summary.RLResult <- function(object, ...)
  {
    val <- summary(object@data)
    val[["nLinks"]] <- sum(chunkify(function(x) sum(x=="L"))(object@prediction))
    val[["nNonLinks"]] <- sum(chunkify(function(x) sum(x=="N"))(object@prediction))
    val[["nPossibleLinks"]] <- sum(chunkify(function(x) sum(x=="P"))(object@prediction))
    class(val) <- c("summaryRLResult", class(val))
    val
  }
#)

print.summaryRLResult <- function(x, ...)
{
  cat("RLResult object\n")

  cat("\n")
  NextMethod("print", object="x"  )
  cat("\n")
  cat("Number of detected links:", x$nLinks, "\n")
  cat("Number of detected possible links:", x$nPossibleLinks, "\n")
  cat("Number of detected non-links:", x$nNonLinks, "\n")
}


# get accuracy, alpha-error, beta-error etc. from result object
setGeneric(
  name = "getErrorMeasures",
  def = function(object, ...) standardGeneric("getErrorMeasures")
)

# method for 'old' S3 class
setMethod(
  f = "getErrorMeasures",
  signature = "RecLinkResult",
  definition = function(object, ...)
  {
    TP=length(which(object$pairs$is_match & object$prediction=="L")) # true positive
    FP=length(which(!object$pairs$is_match & object$prediction=="L")) # false positive
    TN=length(which(!object$pairs$is_match & object$prediction=="N")) # true negative
    FN=length(which(object$pairs$is_match & object$prediction=="N")) # false negative

    return(list(
      alpha=FN/(TP+FN),
      beta=FP/(TN+FP),
      accuracy=(TP+TN)/(TP+TN+FP+FN),
      precision=TP/(TP+FP),
      sensitivity=TP/(TP+FN),
      specificity=TN/(TN+FP),
      ppv=TP/(TP+FP),
      npv=TN/(TN+FN)
    ))
  }
)

# method for 'new' S4 class (for big data objects)
setMethod(
  f = "getErrorMeasures",
  signature = "RLResult",
  definition = function(object, ...)
  {

    # TP: true positive, FP: false positive, TN: true negative,
    # FN: false negative, PM: possible links that are matches,
    # PN: possible links that are non-matches
    nUnknown <- nMatch <- TP <- FP <- PM <- PN <- numeric(1)
    
    resultTable <- ffdf(is_match=object@data@pairs$is_match, class=object@prediction)
    ffrowapply(
    {
      is_match <- resultTable[i1:i2,"is_match"]==1
      is_link <-  resultTable[i1:i2,"class"]=="L"
      is_possible <- resultTable[i1:i2,"class"]=="P"
      TP <- TP + sum(is_match & is_link, na.rm=TRUE)
      FP <- FP + sum(!is_match & is_link, na.rm=TRUE)
      PM <- PM + sum(is_match & is_possible, na.rm=TRUE)
      PN <- PN + sum(!is_match & is_possible, na.rm=TRUE)
      nMatch <- nMatch + sum(is_match, na.rm=TRUE)
      nUnknown <- nUnknown + sum(is.na(is_match))
    }, X=resultTable)


    FN <- nMatch - TP - PM
    TN <- nrow(object@data@pairs) - TP - FN - FP - PM - PN - nUnknown
    return(list(
      alpha=FN/(TP+FN),
      beta=FP/(TN+FP),
      accuracy=(TP+TN)/(TP+TN+FP+FN),
      precision=TP/(TP+FP),
      sensitivity=TP/(TP+FN),
      specificity=TN/(TN+FP),
      ppv=TP/(TP+FP),
      npv=TN/(TN+FN)
    ))
  }
)

# wrapper for new S4 method for backward compatibility
errorMeasures <- function(result)
{
  if (!("RecLinkResult" %in% class(result)))
      stop(sprintf("Wrong type for result: %s!", class(result)))
  getErrorMeasures(result)
}


setGeneric(
  name = "getTable",
  def = function(object, ...) standardGeneric("getTable")
)

# constructs a contengency table of matches
setMethod(
  f = "getTable",
  signature = "RecLinkResult",
  definition = function(object, ...)
  {
    TP=length(which(object$pairs$is_match & object$prediction=="L")) # true positive
    FP=length(which(!object$pairs$is_match & object$prediction=="L")) # false positive
    TN=length(which(!object$pairs$is_match & object$prediction=="N")) # true negative
    FN=length(which(object$pairs$is_match & object$prediction=="N")) # false negative

    tab <- table(as.logical(object$pairs$is_match),object$prediction,
            dnn=list("true status","classification"),useNA="ifany")
    # if "NA" row appears in the table (for pairs with unknown true status),
    # put it in the middle
    if (nrow(tab) == 3)
      tab[c(1,3,2),]
    else
      tab
  }
)

setMethod(
  f = "getTable",
  signature = "RLResult",
  definition = function(object, ...)
  {
    tab <- table(object@data@pairs$is_match, object@prediction,
      useNA = "ifany")
    names(dimnames(tab)) <- c("true status", "classification")
    dimnames(tab)[dimnames(tab)=="1"] <- "TRUE"
    dimnames(tab)[dimnames(tab)=="0"] <- "FALSE"

    # if "NA" row appears in the table (for pairs with unknown true status),
    # put it in the middle
    if (nrow(tab) == 3)
      tab[c(1,3,2),]
    else
      tab

  }
)

