################ Correlation ############################
#' Correlation
#'
#' Returns the Pearson correlation coefficient between two \link{ml.data.frame} fields.
#'
#' The function eliminates all pairs for which either the first element or the second
#' element is empty. After the elimination, if the length of the input is less than 2,
#' the function returns the empty sequence. After the elimination, if the standard
#' deviation of the first column or the standard deviation of the second column is 0,
#' the function returns the empty sequence.
#'
#' @param x a ml.data.frame field.
#' @param y a ml.data.frame field
#' @param use not used currently
#' @param method not used currently
#' @return The correlation coefficient
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # return the correlation
#'  cor(mlIris$Sepal.Length, mlIris$Petal.Length)
#' }
#' @export
setMethod(f="cor", signature=c(x="ml.col.def",y="ml.col.def"),

          function(x,y,use = NULL,method = NULL ) {

            # use
            if (!missing(use) && !is.null(use))
             stop(simpleError("use option is not implemented yet"))

             # method
            if (!missing(method) && !is.null(method))
               stop(simpleError("method option is not implemented yet"))

            if(x@.parent@.name!=y@.parent@.name) {
               stop("Cannot combine two columns from different ml.data.frames.")
            }
            if(x@.data_type!="number" || y@.data_type != "number") {
              stop("Can only use columns of number type")
            }
            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns , '"},"', y@.name, '":{"fieldDef":"',y@.expr, '","orgField":"', y@.org_name, '","orgFormat":"', y@.format, '","xmlns":"', y@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.correlation", "noindex": "math.correlation"}'
            return(.ml.stat.func(x@.parent, fields, func))
        }
)
################ Correlation Matrix ############################
#' Correlation Matrix
#'
#' Returns the Pearson correlation coefficient matrix of all numeric fields in a \link{ml.data.frame}
#'
#' The function eliminates all fields pairs for which either the first element or the second
#' element is empty. After the elimination, if the length of the input is less than 2,
#' the function returns the empty sequence. After the elimination, if the standard
#' deviation of the first column or the standard deviation of the second column is 0,
#' the function returns the empty sequence.
#'
#' @param x a \link{ml.data.frame}
#' @param y not used when doing a matrix
#' @param use not implemented
#' @param method not implemented
#' @return The correlation coefficient matrix
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # return the correlation matrix
#'  cor(mlIris)
#' }
#' @export
setMethod(f="cor", signature=c(x="ml.data.frame"),

          function(x,y = NULL,use = NULL,method = NULL ) {

            # use
            if (!missing(use) && !is.null(use))
              stop(simpleError("use option is not implemented yet"))

            # method
            if (!missing(method) && !is.null(method))
              stop(simpleError("method option is not implemented yet"))
            # get correlation matrix data
            corMatResult <- .ml.matrix(x, "correlation")
            # create the matrix
            corMat <- matrix(1:(length(corMatResult)),nrow=length(corMatResult),ncol=length(corMatResult),dimnames = list(names(corMatResult),names(corMatResult)),byrow=T)
            for(i in 1:length(corMatResult)) {
              for(j in 1:length(corMatResult[[i]])) {
                if (is.null(corMatResult[[i]][[j]])) {
                  corMat[i,j] <- NA
                } else {
                  corMat[i,j] <- corMatResult[[i]][[j]]
                }
              }
            }
            corMat
          }
)
################ Covariance ############################
#' Covariance
#'
#' Returns the sample covariance of two variables, \link{ml.data.frame} fields.
#'
#'The function eliminates all pairs for which either the first element or the second
#'element is empty. After the elimination, if the length of the input is less than 2,
#'the function returns the empty sequence.
#'
#' @param x a ml.data.frame field.
#' @param y a ml.data.frame field
#' @param use not implemented
#' @param method not implemented
#' @return The sample covariance
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn,collection = "iris")
#'  # return the Covariance
#'  cov(mlIris$Sepal.Length, mlIris$Petal.Length)
#' }
#' @export
setMethod(f="cov", signature=c(x="ml.col.def",y="ml.col.def"),

          function(x,y,use = NULL,method = NULL ) {

            # use
            if (!missing(use) && !is.null(use))
              stop(simpleError("use option is not implemented yet"))

            # method
            if (!missing(method) && !is.null(method))
              stop(simpleError("method option is not implemented yet"))

            if(x@.parent@.name!=y@.parent@.name) {
              stop("Cannot combine two columns from different ml.data.frames.")
            }
            if(x@.data_type!="number" || y@.data_type != "number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns , '"},"', y@.name, '":{"fieldDef":"',y@.expr, '","orgField":"', y@.org_name, '","orgFormat":"', y@.format, '","xmlns":"', y@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.covariance", "noindex": "math.covariance"}'
            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Population Covariance ############################
#' Population Covariance
#'
#' Returns the population covariance of two variables, \link{ml.data.frame} fields.
#'
#' The function eliminates all pairs for which either the first element or the
#' second element is empty. After the elimination, if the length of the input is 0,
#' the function returns the empty sequence.
#'
#' @param x a ml.data.frame field.
#' @param y a ml.data.frame field
#' @return The population covariance
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # return the population covariance
#'  cov.pop(mlIris$Sepal.Length, mlIris$Petal.Length)
#' }
#' @export
cov.pop <- function(x,y) {

  if(x@.parent@.name!=y@.parent@.name) {
    stop("Cannot combine two columns from different ml.data.frames.")
  }
  if(x@.data_type!="number" || y@.data_type != "number") {
    stop("Can only use columns of number type")
  }

  fields <- "{"
  fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns, '"},"', y@.name, '":{"fieldDef":"',y@.expr, '","orgField":"', y@.org_name, '","orgFormat":"', y@.format, '","xmlns":"', y@.xmlns ,'"}' ,sep='')
  fields <- paste(fields, '}', sep='')
  func <- '{"index":"cts.covarianceP", "noindex": "math.covarianceP"}'

  return(.ml.stat.func(x@.parent, fields, func))
}
################ Variance ############################
#' Variance
#'
#' Returns the sample variance of a \link{ml.data.frame} field.
#'
#' The function returns a empty value if the number of rows of the ml.data.frame
#' that x belongs to is less than 2.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not used currently
#' @return The sample variance
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # return the variance
#'  var(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="var", signature=c(x="ml.col.def"),

          function(x,na.rm = FALSE ) {

            # use
            if (na.rm )
              stop(simpleError("na.rm option is not implemented yet"))

            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.variance", "noindex": "math.variance"}'

            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Population Variance ############################
#' Population variance
#'
#' Returns the population variance of of a \link{ml.data.frame} field.
#'
#' The function returns a empty value if the number of rows of the ml.data.frame
#' that x belongs to is less than 2.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not used currently
#' @return The population variance
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # population variance
#'  var.pop(mlIris$Sepal.Length)
#' }
#' @export
var.pop <- function(x,na.rm = FALSE ) {

  # use
  if (na.rm )
    stop(simpleError("na.rm option is not implemented yet"))

  if(x@.data_type!="number") {
    stop("Can only use columns of number type")
  }

  fields <- "{"
  fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
  fields <- paste(fields, '}', sep='')
  func <- '{"index":"cts.varianceP", "noindex": "math.varianceP"}'
  return(.ml.stat.func(x@.parent, fields, func))
}
################ Standard Deviation ############################
#' Standard Deviation
#'
#' Returns the sample standard deviation of a \link{ml.data.frame} field.
#'
#' The function returns a empty value if the number of rows of the ml.data.frame
#' that x belongs to is less than 2.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not used currently
#' @return The sample standard deviation
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # standard deviation
#'  sd(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="sd", signature=c(x="ml.col.def"),
          function(x,na.rm=NULL) {

            # na.rm
            if (!missing(na.rm) && !is.null(na.rm))
              warning(simpleError("na.rm option is not implemented yet"))

            if(nrow(x@.parent)==0) {
              stop("No rows in ml.data.frame.")
            }
            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.stddev", "noindex": "math.stddev"}'

            return(.ml.stat.func(x@.parent, fields, func))

          }
)

################ Standard Deviation population ############################
#' Standard Deviation of a population
#'
#' Returns the sample standard deviation of a population.
#'
#' @param x a ml.data.frame field.
#' @return The sample standard deviation of a population.
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # standard deviation
#'  sd.pop(mlIris$Sepal.Length)
#' }
#' @export
sd.pop <- function(x) {

  if(nrow(x@.parent)==0) {
    stop("No rows in ml.data.frame.")
  }
  if(x@.data_type!="number") {
    stop("Can only use columns of number type")
  }

  fields <- "{"
  fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
  fields <- paste(fields, '}', sep='')
  func <- '{"index":"cts.stddevP", "noindex": "math.stddevP"}'
  return(.ml.stat.func(x@.parent, fields, func))

}
################ Median ############################
#' Median
#'
#' Returns the median of a \link{ml.data.frame} field.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not currently used.
#' @return The median
#' @examples
#' \dontrun{
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # median
#'  median(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="median", signature=c(x="ml.col.def"),

          function(x, na.rm = FALSE) {

            # use
            if (na.rm)
              warning(simpleError("na.rm option is not implemented yet"))

            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.median", "noindex": "math.median"}'
            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Mean ############################
#' Mean
#'
#' Returns the mean of a \link{ml.data.frame} field.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not currently used.
#' @return The mean
#' @examples
#' \dontrun{
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # mean
#'  mean(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="mean", signature=c(x="ml.col.def"),

          function(x, na.rm = FALSE) {

            # use
            if (na.rm)
              warning(simpleError("na.rm option is not implemented yet"))

            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.avgAggregate", "noindex": "fn.avg"}'
            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Sum ############################
#' Sum
#'
#' Returns the sum of a \link{ml.data.frame} field.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not currently used.
#' @return The sum
#' @examples
#' \dontrun{
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # sum
#'  sum(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="sum", signature=c(x="ml.col.def"),

          function(x, na.rm = FALSE) {

            # use
            if (na.rm)
              warning(simpleError("na.rm option is not implemented yet"))

            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.sumAggregate", "noindex": "fn.sum"}'

            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Max ############################
#' Max
#'
#' Returns the maximum value of a \link{ml.data.frame} field.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not currently used.
#' @return The maximum value
#' @examples
#' \dontrun{
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(locConn, collection = "iris")
#'  # max
#'  max(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="max", signature=c(x="ml.col.def"),

          function(x, na.rm = FALSE) {

            # use
            if (na.rm)
              warning(simpleError("na.rm option is not implemented yet"))

            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.max", "noindex": "fn.max"}'
            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Min ############################
#' Min
#'
#' Returns the minimum value of a ml.data.frame field.
#'
#' @param x a ml.data.frame field.
#' @param na.rm not currently used.
#' @return The minimum value
#' @examples
#' \dontrun{
#'  ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(collection = "iris")
#'  # min
#'  min(mlIris$Sepal.Length)
#' }
#' @export
setMethod(f="min", signature=c(x="ml.col.def"),

          function(x, na.rm = FALSE) {

            # use
            if (na.rm)
              warning(simpleError("na.rm option is not implemented yet"))

            if(x@.data_type!="number") {
              stop("Can only use columns of number type")
            }

            fields <- "{"
            fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
            fields <- paste(fields, '}', sep='')
            func <- '{"index":"cts.min", "noindex": "fn.min"}'
            return(.ml.stat.func(x@.parent, fields, func))
          }
)
################ Percentile ############################
# Currently no implemented since the range index version,cts.percentile, does not
# follow how the majority  of the cts functions works.
percentile <- function(x, p) {
 # use
  if(x@.data_type!="number") {
    stop("Can only use columns of number type")
  }

  fields <- "{"
  fields <- paste(fields, '"',x@.name , '":{"fieldDef":"',x@.expr ,'","orgField":"', x@.org_name, '","orgFormat":"', x@.format, '","xmlns":"', x@.xmlns ,'"}' ,sep='')
  fields <- paste(fields, '}', sep='')
  func <- '{"index":"cts.percentile", "noindex": "math.percentile"}'
  return(.ml.stat.func(x@.parent, fields, func))
}

################ Summary ############################
#' ml.data.frame Summaries
#'
#' @param object an ml.data.frame object
#' @param digits integer, used for number formatting
#' @param maxsum not used.
#' @param ... not used.
#' @aliases summary
#' @export
setMethod(f="summary", signature=c("ml.data.frame"),
          function (object,digits=max(3L, getOption("digits") -3L), maxsum = 7L, ...) {
            mlDf<-object

            options(scipen=999)

            summaryTbl <- list()
            labelNum<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
            #labelCat <- c("Length","Class","Mode")
            # Get statsitics per field
            sumResult <- .ml.matrix(mlDf, "summary")
            # since the result is a list with lists (one for each field) we need to transform it to
            # a list with one value, all the statistics, for each field.
             for (i in 1:length(sumResult)) {
              if (sumResult[[i]]$valType == 'NUMERIC') {
                values <- c(sumResult[[i]]$min, sumResult[[i]]$q1, sumResult[[i]]$median,sumResult[[i]]$mean, sumResult[[i]]$q3, sumResult[[i]]$max)
                values <- paste0(format(labelNum), ":", format(values,digits=digits,nsmall=3), "  ")
              } else {
                labelCat <- names(sumResult[[i]]$levels)
                #values <- c(sumResult[[i]]$length, "character", "character")
                values <- unlist(sumResult[[i]]$levels, use.names = FALSE)
                values <- paste0(format(labelCat), ":", format(values), "  ")
              }
              # make sure that the length is same for all values
              length(values) <- maxsum
              summaryTbl <- c(summaryTbl,list(c(values)))

            }
            # we need to add names on summaryTbl so we can do the intersect later
            names(summaryTbl)<-names(sumResult)
            # We are only intrested in keeping the fields that exists in our ml.data.frame
            summaryTbl <- summaryTbl[intersect(mlDf@.col.name, names(sumResult))]

            summaryTbl<-unlist(summaryTbl)

            # Since there is a possibility that our summary result has more fields than our
            # ml.data.frame we need to use the number of fields in mlDf
            dim(summaryTbl) <- c(maxsum, length(mlDf@.col.name))

            # adjusting the output, so field names are aligned with the statistic information
            blanks <- paste(character(max(10, na.rm = TRUE) + 2L),collapse = " ")
            pad <- floor(nchar(summaryTbl[1,])/2 - nchar(intersect(mlDf@.col.name,names(sumResult)))/2)
            names <- paste0(substring(blanks, 1, pad), intersect(mlDf@.col.name,names(sumResult)))
            dimnames(summaryTbl)<-list(rep.int("",maxsum),names)
            attr(summaryTbl, "class") <- c("table")
            summaryTbl
          }
)
