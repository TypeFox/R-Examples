##################################################################################
##                      METHODS FOR S3 GENERICS                                 ##
##################################################################################

##CAUTION NOTE:
##(1) Please DO NOT access to the S4 objects by using the slots this causes inconsistent
## behaviour, please use the correspong RLum-methods instead!
##
##(2) Especially, please DO NOT include S3-methods for which no S4-method is implemented! Especially
##for coercing.
##
##(3) Finally, what ever you want to implemnt, check whether a S4-method exists, it should
##be just passed to this methods, not the opposite, otherwise this will yield in undesired behaviour.
##
##TODO: For this S3 generics so far no proper documentation exists ... we should consider
##to provide an overview within a separat document, as it becomes otherwise rather
##complicated for beginners to work with the documentation.
##


## -------------------- INTRODUCED WITH 0.5.0 ----------------------- ##


#' methods_RLum
#'
#' Methods for S3-generics implemented for the package 'Luminescence'.
#' This document summarises all implemented S3-generics. The name of the function
#' is given before the first dot, after the dot the name of the object that is supported by this method
#' is given, e.g. \code{plot.RLum.Data.Curve} can be called by \code{plot(object, ...)}, where
#' \code{object} is the \code{RLum.Data.Curve} object.
#'
#' The term S3-generics sounds complicated, however, it just means that something has been implemented
#' in the package to increase the usability for users new in R and who are not familiar with the
#' underlying \code{RLum}-object structure of the package. The practical outcome is that
#' operations and functions presented in standard books on R can be used without knowing the specifica
#' of the R package 'Luminescence'. For examples see the example section.
#'
#' @param x \code{\linkS4class{RLum}} (\bold{required}): input opject
#'
#' @param object \code{\linkS4class{RLum}} (\bold{required}): input opject
#'
#' @param y \code{\link{integer}} (optional): the row index of the matrix, data.frame
#'
#' @param z \code{\link{integer}} (optional): the column index of the matrix, data.frame
#'
#' @param i \code{\link{character}} (optional): name of the wanted record type or data object
#'
#' @param drop \code{\link{logical}} (with default): keep object structure or drop it
#'
#' @param row.names \code{\link{logical}} (with default): enables or disables row names (\code{as.data.frame})
#'
#' @param recursive \code{\link{logical}} (with default): enables or disables further subsetting (\code{unlist})
#'
#' @param optional \code{\link{logical}} (with default): logical. If TRUE, setting row names and
#' converting column names (to syntactic names: see make.names) is optional (see \code{\link[base]{as.data.frame}})
#'
#' @param ... further arguments that can be passed to the method
#'
#' @note \code{methods_RLum} are not really new functions, everything given here are mostly just
#' surrogates for existing functions in the package.
#'
#' @examples
#'
#' ##load example data
#' data(ExampleData.RLum.Analysis, envir = environment())
#'
#' @name methods_RLum
NULL

####################################################################################################
# methods for generic: plot()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
plot.RLum.Results <- function(x, y, ...) plot_RLum(object = x, ...)

#' @rdname methods_RLum
#' @export
plot.RLum.Analysis <- function(x, y, ...) plot_RLum(object = x, ...)

#' @rdname methods_RLum
#' @export
plot.RLum.Data.Curve <- function(x, y, ...) plot_RLum(object = x, ...)

#' @rdname methods_RLum
#' @export
plot.RLum.Data.Spectrum <- function(x, y, ...) plot_RLum(object = x, ...)

#' @rdname methods_RLum
#' @export
plot.RLum.Data.Image <- function(x, y, ...) plot_RLum(object = x, ...)

#' @rdname methods_RLum
#' @export
plot.Risoe.BINfileData <- function(x, y, ...) plot_Risoe.BINfileData(BINfileData = x, ...)

####################################################################################################
# methods for generic: hist()
# ##################################################################################################

#' @rdname methods_RLum
#' @export
hist.RLum.Results <- function(x, ...) plot_Histogram(data = x, ...)

#' @rdname methods_RLum
#' @export
hist.RLum.Data.Image <- function(x, ...) hist(x =get_RLum(x)@data@values, ...)

#' @rdname methods_RLum
#' @export
hist.RLum.Data.Curve <- function(x, ...) hist(as(get_RLum(x),"matrix")[,2])

#' @rdname methods_RLum
#' @export
hist.RLum.Analysis <- function(x, ...) lapply(1:length_RLum(x), function(z){
  hist(as(get_RLum(x, record.id = z, ...),"matrix")[,2])})

####################################################################################################
# methods for generic: summary()
# ##################################################################################################
# methods for generic: summary()
#' @rdname methods_RLum
#' @export
summary.RLum.Results <- function(object, ...) get_RLum(object = object, ...)

#' @rdname methods_RLum
#' @export
summary.RLum.Analysis <- function(object, ...) lapply(object@records, function(x) summary(x@data))

#' @rdname methods_RLum
#' @export
summary.RLum.Data.Image <- function(object, ...) summary(object@data@data@values)

# summary.RLum.Data.Spectrum <- function(object, ...)

#' @rdname methods_RLum
#' @export
summary.RLum.Data.Curve <- function(object, ...) summary(object@data, ...)

####################################################################################################
# methods for generic: length()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
length.RLum.Results <- function(x, ...) length_RLum(x)

#' @rdname methods_RLum
#' @export
length.RLum.Analysis <- function(x, ...) length_RLum(x)

#' @rdname methods_RLum
#' @export
length.RLum.Data.Curve <- function(x, ...) length_RLum(x)

#' @rdname methods_RLum
#' @export
length.Risoe.BINfileData <- function(x, ...) length(x@METADATA$ID)

####################################################################################################
# methods for generic: dim()
# ##################################################################################################
# methods for generic: dim()
#' @rdname methods_RLum
#' @export
dim.RLum.Data.Curve <- function(x) dim(as(x, "matrix"))

#' @rdname methods_RLum
#' @export
dim.RLum.Data.Spectrum <- function(x) dim(as(x, "matrix"))

####################################################################################################
# methods for generic: rep()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
rep.RLum <- function(x, ...) replicate_RLum(x, ...)

####################################################################################################
# methods for generic: name()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
names.RLum.Data.Curve <- function(x, ...) names_RLum(x)

#' @rdname methods_RLum
#' @export
names.RLum.Data.Spectrum <- function(x, ...) names_RLum(x)

#' @rdname methods_RLum
#' @export
names.RLum.Data.Image <- function(x, ...) names_RLum(x)

#' @rdname methods_RLum
#' @export
names.RLum.Analysis <- function(x, ...) names_RLum(x)

#' @rdname methods_RLum
#' @export
names.RLum.Results <- function(x, ...) names_RLum(x)

#' @rdname methods_RLum
#' @export
names.Risoe.BINfileData <- function(x)  as.character(x@METADATA$LTYPE)

####################################################################################################
# methods for generic: row.name()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
row.names.RLum.Data.Spectrum <- function(x, ...) rownames(as(x, "matrix"))

####################################################################################################
# methods for generic: as.data.frame()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
as.data.frame.RLum.Data.Curve <- function(x, row.names = NULL, optional = FALSE, ...) as(x, "data.frame")

#' @rdname methods_RLum
#' @export
as.data.frame.RLum.Data.Spectrum <- function(x,  row.names = NULL, optional = FALSE, ...) as(x, "data.frame")
# for RLum.Results ... makes no sense and may yield in unpredictable behaviour

####################################################################################################
# methods for generic: as.list()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
as.list.RLum.Results <- function(x, ...) as(x, "list")

#' @rdname methods_RLum
#' @export
as.list.RLum.Data.Curve <- function(x, ...) as(x, "list")

#' @rdname methods_RLum
#' @export
as.list.RLum.Analysis <- function(x, ...) as(x, "list")

####################################################################################################
# methods for generic: as.matrix()
# ##################################################################################################
#' @rdname methods_RLum
#' @export
as.matrix.RLum.Data.Curve <- function(x, ...) as(x, "matrix")

#' @rdname methods_RLum
#' @export
as.matrix.RLum.Data.Spectrum <- function(x, ...) as(x, "matrix")
# for RLum.Results ... makes no sense and may yield in unpredictable behaviour

####################################################################################################
# methods for generic: merge()
####################################################################################################
#' @rdname methods_RLum
#' @export
merge.RLum <- function(x, y, ...) merge_RLum(append(list(...), values = c(x, y)))

####################################################################################################
# methods for generic: unlist()
####################################################################################################
#' @rdname methods_RLum
#' @method unlist RLum.Analysis
#' @export
unlist.RLum.Analysis <- function(x, recursive = TRUE, ...){

  temp <- get_RLum(object = x, recursive = recursive, ... )
  if(recursive){
    unlist(lapply(1:length(temp), function(x){
      get_RLum(temp)
    }), recursive = FALSE)

  }else{
    return(temp)

  }

}

####################################################################################################
# methods for generic: `+`
####################################################################################################
#' @rdname methods_RLum
#'
#' @examples
#'
#' ##combine curve is various ways
#' curve1 <- IRSAR.RF.Data[[1]]
#' curve2 <-  IRSAR.RF.Data[[1]]
#' curve1 + curve2
#' curve1 - curve2
#' curve1 / curve2
#' curve1 * curve2
#'
#' @export
`+.RLum.Data.Curve` <- function(x, y) merge_RLum(list(x, y), merge.method = "sum")

####################################################################################################
# methods for generic: `-`
####################################################################################################
#' @rdname methods_RLum
#' @export
`-.RLum.Data.Curve` <- function(x, y) merge_RLum(list(x, y), merge.method = "-")

####################################################################################################
# methods for generic: `*`
####################################################################################################
#' @rdname methods_RLum
#' @export
`*.RLum.Data.Curve` <- function(x, y) merge_RLum(list(x, y), merge.method = "*")

####################################################################################################
# methods for generic: `/`
####################################################################################################
#' @rdname methods_RLum
#' @export
`/.RLum.Data.Curve` <- function(x, y) merge_RLum(list(x, y), merge.method = "/")

####################################################################################################
# methods for generic: `[`
####################################################################################################
#' @rdname methods_RLum
#' @export
`[.RLum.Data.Curve` <- function(x,y,z, drop = TRUE) {as(x, "matrix")[y,z, drop = drop]}

#' @rdname methods_RLum
#' @export
`[.RLum.Data.Spectrum` <- function(x,y,z, drop = TRUE) {as(x, "matrix")[y,z, drop = drop]}

#' @rdname methods_RLum
#' @export
`[.RLum.Data.Image` <- function(x,y,z, drop = TRUE) {as(x, "matrix")[y,z, drop = drop]}

#' @rdname methods_RLum
#' @export
`[.RLum.Analysis` <- function(x, i, drop = FALSE) {
  if (is(i, "character")) {
    get_RLum(x, recordType = i, drop = drop)

  } else{
    get_RLum(x, record.id = i, drop = drop)

  }
}

#' @rdname methods_RLum
#' @export
`[.RLum.Results` <- function(x, i, drop = TRUE) {get_RLum(x, data.object = i, drop = drop)}


####################################################################################################
# methods for generic: `[[`
####################################################################################################
#' @rdname methods_RLum
#' @export
`[[.RLum.Analysis` <- function(x, i) {
  if (is(i, "character")) {
    get_RLum(x, recordType = i)

  } else{
    get_RLum(x, record.id = i)

  }
}

#' @rdname methods_RLum
#' @export
`[[.RLum.Results` <- function(x, i) {get_RLum(x, data.object = i)}

####################################################################################################
# methods for generic: `$`
####################################################################################################
#' @rdname methods_RLum
#'
#' @examples
#'
#' ##`$` access curves
#' IRSAR.RF.Data$RF
#'
#' @export
`$.RLum.Analysis` <- function(x, i) {get_RLum(x, recordType = i)}

#' @rdname methods_RLum
#' @export
`$.RLum.Results` <- function(x, i) {get_RLum(x, data.object = i)}
