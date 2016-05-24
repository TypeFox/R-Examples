#' Mining Association rules and Frequent Itemsets
#'
#' Mine frequent itemsets or association rules using MarkLogic Server built in Range Index functions.
#' The function require that there is a Range Index on the underlying field of itemField, a range indexe can
#' be created with the \link{ml.add.index} function. It will return a object that is of class rules or itemsets
#' as defined in the arules package. It will need the arules package installed.
#'
#' The frequent itemset and association rules extraction method is using the same method as the Apriori
#' algorithm by first identify all 1-n itemsets that satisfy the support threshold and based on these
#' extract rules that satisfy the confidence threshold.
#'
#' It is depended on that there are a Range Index on the underlying field for the itemField.
#' Information about the name of the field can be shown by \code{mlDataFrame$itemField}, where mlDataFrame
#' is a ml.data.frame object and itemField is the name of the field.
#'
#' @param data an \link{ml.data.frame} object
#' @param itemField a ml.data.frame field which is the field that the itemsets will be created of. The underlying field needs to have a Range Index defined.
#' @param support a numeric value for the minimal support of an item set (default: 0.5)
#' @param confidence a numeric value for the minimal confidence of rules/association hyperedges (default: 0.8)
#' @param maxlen an integer value for the maximal number of items per item set (default: 5)
#' @param target a character string indicating the type of association mined. One of "frequent itemsets" or "rules", default is "rules"
#' @return Returns an object of class rules or itemsets.

#' @export
ml.arules <- function(data, itemField, support = 0.5, confidence = 0.8, maxlen = 5, target = "rules") {
#"frequent itemsets" "maximally frequent itemsets" "closed frequent itemsets" "rules" (only available for Apriori) "hyperedgesets"

  # need to check for the arules package ...
  if (!requireNamespace("arules", quietly = TRUE)) {
    stop("arules is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  conn <- data@.conn
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  queryComArgs <- data@.queryArgs

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.arules", sep="")

  queryArgs <- c(queryComArgs, 'rs:supp'=support, 'rs:conf'=confidence, 'rs:maxlen'=maxlen, 'rs:target'=target)

  if (!inherits(itemField, "ml.col.def")) {
    stop("itemField parameter must be a valid ml.data.frame field.")
  }
  fields <- "{"
  fields <- paste(fields, '"',itemField@.name , '":{"fieldDef":"',itemField@.expr ,'","orgField":"', itemField@.org_name, '","orgFormat":"', itemField@.format , '"}',sep='')
  fields <- paste(fields, '}', sep='')
  queryArgs <- c(queryArgs, 'rs:fields'=fields)

  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())

  rContent <- content(response) #, as = "text""
  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  if (target == "frequent itemsets") {
    # check what we should return ..
    # create a arules itemset...
    # first generate a list with all itemsets ...
    rItemsets <- rContent$itemsets
    extrItemsets <- list()
    supportDf <- data.frame(support=numeric())
    if (length(rItemsets) == 0) {
      warning("No itemsets returned. You might need to decrease the support or confidence thresholds.")
      return()
    }
    sets <- 1
    for (i in 1:length(rItemsets)) {
      for (j in 1:length(rItemsets[[i]])) {
        extrItemsets[[sets]] <- as.vector(unlist(rItemsets[[i]][[j]]$'itemSet'))
        supportDf <- rbind(supportDf,data.frame(support=rItemsets[[i]][[j]]$'support'))
        sets <- sets + 1
      }
    }
    # we need arules loaded here!
    result <- new("itemsets")
    # seems to work but not fully...
    # itemsetInfo not right...
    result@items <- as(extrItemsets, "itemMatrix")
    result@quality <- supportDf
    validObject(result@items@data)
  } else if (target == "rules") {
    # get rules ...
    extrLhs <- list()
    extrRhs <- list()
    qualityRule <- data.frame(support=numeric(),confidence=numeric(),lift=numeric())
    rRules <- rContent$rules
    if (length(rRules) == 0) {
      warning("No rules returned. You might need to decrease the support or confidence thresholds.")
      return()
    }
    for (i in 1:length(rRules)) {
      if (length(rRules[[i]]$'lhs') > 0) {
        extrLhs[[i]] <- as.vector(unlist(rRules[[i]]$'lhs')) #c(extrLhs, toString(rRules[[i]]$'lhs'))
      } else {
        extrLhs[[i]] <- ""
      }
      if (length(rRules[[i]]$'rhs') > 0) {
        extrRhs[[i]] <- as.vector(unlist(rRules[[i]]$'rhs')) #c(extrRhs, toString(rRules[[i]]$'rhs'))
      }else{
        extrRhs[[i]] <- ""
      }
       qualityRule <- rbind(qualityRule, data.frame(support = rRules[[i]]$'support',confidence= rRules[[i]]$'confidence',lift= rRules[[i]]$'lift'))
     }
     # we need arules loaded here!
    suppressWarnings(result <- new("rules", lhs=as(extrLhs, "itemMatrix"), rhs=as(extrRhs, "itemMatrix"), quality=qualityRule))
     validObject(result@lhs@data)
     validObject(result@rhs@data)
  }
  ## add some reflectance
  call <- match.call()
  result@info <- list(data = call$data,
                       ntransactions = data@.nrows,
                       support = support,
                       confidence = confidence)
  result
}
