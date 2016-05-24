#' Aggregate numeric,Date and categorical variables by an ID
#'
#' The \code{Aggregate} function (not to be confounded with aggregate) prepares a data frame for merging by computing the sum, mean and variance of all continuous (integer and numeric) variables by a given ID variable. For all categorical variabes (character and factor), it creates dummies and subsequently computes the sum by a given ID variable. For all Date variables, it computes recency and duration by a given ID and end date variable. This functions aims at maximum information extraction with a minimum amount of code.
#'
#' @param x A data frame without the ID. Categorical variables have to be of type character or factor and continuous variables have to be of type integer or numeric.
#' @param by A vector containing IDs.
#' @param end_ind A Date object, or something which can be coerced by \code{as.Date(origin, ...)} to such an object. If not specified,  we take the \code{Sys.Date()} as end date.
#' @param format A character string. If not specified, the ISO 8601 international standard which expresses a day "\%Y-\%m-\%d" is taken.
#' @param ... Extra parameters to be passed to the \code{dummy} function in the \code{dummy} package.
#'
#' @return A data frame with the aforementioned variables aggregated by the given ID variables
#'
#' @author Authors: Matthias Bogaert, Michel Ballings, Dirk Van den Poel, Maintainer: \email{matthias.bogaert@@UGent.be}
#' @examples
#' # Example
#' # Create some data
#' data <- data.frame(V1=as.factor(c('yes','no','no','yes','yes','yes','yes')),
#'                    V2=as.character(c(1,2,3,4,4,4,4)),V3=c(1:7),V4=as.numeric(c(7:1)),
#'                    V5 = as.Date(as.Date('2014-12-03'):as.Date('2014-12-09'), origin = "1970-01-01")
#'                    )
#' ID=as.character(c(1,1,1,1,2,2,2))

#' Aggregate(x=data,by=ID)
#'
#' # Examples of how to use the ... argument. See package dummy for details.
#' # library(dummmy)
#' # Aggregate(x=data,by=ID,object=categories(data))
#' # Aggregate(x=data,by=ID,p=2)

Aggregate <- function (x, by, end_ind = Sys.Date(), format = '%Y-%m-%d',  ...) {
    if (!is.data.frame(x)) stop("x need to be a data frame")

    categoricals <- sapply(x, is.factor) | sapply(x, is.character)
    if (any(categoricals == TRUE)) {
      dummies <- dummy(x[,categoricals,drop=FALSE], int=TRUE, ...)
      dummies_sum <- aggregate(dummies, list(ID = by), sum)
      dummies_tail <- aggregate(dummies, list(ID = by), tail, n = 1)
      ID <- dummies_tail$ID
      dummies_sum$ID <- dummies_tail$ID <- NULL
      names(dummies_sum) <- paste(names(dummies_sum), "_sum", sep = "")
      names(dummies_tail) <- paste(names(dummies_tail), "_last", sep = "")
      dummies_df <- data.frame(ID, dummies_sum, dummies_tail)
    }

    numerics <- sapply(x, is.numeric)
    if (any(numerics == TRUE)) {
      numerics_sum <- aggregate(x[,numerics,drop=FALSE], by = list(ID = by), sum)
      numerics_mean <- aggregate(x[,numerics,drop=FALSE], by = list(ID = by), mean)
      numerics_var <- aggregate(x[,numerics,drop=FALSE], by = list(ID = by), var)
      ID <- numerics_sum$ID
      numerics_sum$ID <- numerics_mean$ID <- numerics_var$ID <- NULL

      names(numerics_sum) <- paste(names(numerics_sum), "_sum", sep = "")
      names(numerics_mean) <- paste(names(numerics_mean), "_mean", sep = "")
      names(numerics_var) <- paste(names(numerics_var), "_var", sep = "")
      numerics_df <- data.frame(numerics_sum, numerics_mean, numerics_var)
    }

    dates <- sapply(x,function(z) is(z,"Date"))

    if(any(dates == TRUE)) {
      end_ind <- as.Date(end_ind, format = format)
      x[,dates] <- sapply(x[,dates], function (z) as.Date(z, format = format))
      dates_min <- aggregate(x[,dates, drop = FALSE], by = list(ID = by),min)
      ID <- dates_min$ID
      dates_min$ID <- NULL
      dates_duration <- data.frame(sapply(dates_min, function(z) as.numeric(end_ind - z)))
      names(dates_duration) <- paste(colnames(x)[dates], "_dura", sep = "")

      dates_max <- aggregate(x[,dates, drop = FALSE], by = list(ID = by),max)[,-1]
      dates_recency <- data.frame(sapply(dates_max,function(z) as.numeric(end_ind-z)))
      names(dates_recency) <- paste(colnames(x)[dates], "_rec", sep = "")

      dates_df <- data.frame(dates_duration, dates_recency)
    }

    if (any(dates == TRUE) && any(categoricals == TRUE) && any(numerics == TRUE)) {
      final <- data.frame(dummies_df, numerics_df, dates_df)
    } else if (any(categoricals == TRUE) && any(numerics == TRUE)) {
      final <- data.frame(dummies_df, numerics_df)
    } else if (any(categoricals == TRUE) && any(dates == TRUE)) {
      final <- data.frame(dummies_df,dates_df)
    } else if (any(dates == TRUE) && any(numerics == TRUE)) {
      final <- data.frame(ID, numerics_df, dates_df)
    } else if (any(categoricals == TRUE)) {
      final <- dummies_df
    } else if (any(numerics_df == TRUE)) {
      final <- data.frame(ID,numerics_df)
    } else if (any(dates_df == TRUE)) {
      final <- data.frame(ID, dates_df)
    }
    final
}
