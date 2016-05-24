# Functions for aggregate nodes within communities and communities within 
# collections.

# Many of these helper functions are not used. They are retained as private
# functions here in case the more flexible interface to .AggregateDataFrame
# is resurrected.
.MeanNaRm <- function(x)
{
    # The mean of x, NA removed. If all values are NA, NA is returned.
    if(all(is.na(x))) return (NA)
    else              return (mean(x, na.rm=TRUE))
}

.MedianNaRm <- function(x)
{
    # The median of x, NA removed. If all values are NA, NA is returned.
    if(all(is.na(x))) return (NA)
    else              return (median(x, na.rm=TRUE))
}

.SumNaRm <- function(x)
{
    # The sum of x, NA removed. If all values are NA, NA is returned.
    if(all(is.na(x))) return (NA)
    else              return (sum(x, na.rm=TRUE))
}

.MinNaRm <- function(x)
{
    # The min of x, NA removed. If all values are NA, NA is returned.
    if(all(is.na(x))) return (NA)
    else              return (min(x, na.rm=TRUE))
}

.MaxNaRm <- function(x)
{
    # The max of x, NA removed. If all values are NA, NA is returned.
    if(all(is.na(x))) return (NA)
    else              return (max(x, na.rm=TRUE))
}

.MeanNaRm <- function(x)
{
    # The mean of x, NA removed. If all values are NA, NA is returned.
    if(all(is.na(x))) return (NA)
    else              return (mean(x, na.rm=TRUE))
}

.JoinUnique <- function(x)
{
    # Joins together unique values in x with a comma.
    return (paste(unique(x), collapse=','))
}

.JoinUniqueEmptyRm <- function(x)
{
    # Raises an error if x contains more than one unique value, exclusing NA. 
    # Returns the single unique values otherwise.
    values <- unique(x)
    empty <- is.na(values) | ''==values
    if(all(empty))
    {
        return ('')
    }
    else if(any(empty))
    {
        # The test above is required because this line results in a vector 
        # of length 0 if no values are NA or ''.
        return (.JoinUnique(values[which(!empty)]))
    }
    else
    {
        return (.JoinUnique(values))
    }
}

.ExpectSingleUnique <- function(x)
{
    # Raises an error if x contains more than one unique value. Returns the 
    # single unique values otherwise.
    values <- unique(x)
    if(length(values)>1)
    {
        stop(paste('Contains more than one unique value [', 
                   paste(values, collapse=','), ']', sep=''))
    }
    else
    {
        return (values)
    }
}

.ExpectSingleUniqueEmptyRm <- function(x)
{
    # Raises an error if x contains more than one unique value, exclusing NA. 
    # Returns the single unique values otherwise.
    values <- unique(x)
    empty <- is.na(values) | ''==values
    if(all(empty))
    {
        return ('')
    }
    else if(any(empty))
    {
        # The test above is required because this line results in a vector 
        # of length 0 if no values are NA or ''.
        return (.ExpectSingleUnique(values[which(!empty)]))
    }
    else
    {
        return (.ExpectSingleUnique(values))
    }
}

.AggregateDataFrame <- function(data, aggregate.by, weight.by)
{
    # Aggregates data by the values in aggregate.by. 
    #   data:          a data.frame
    #   aggregate.by:  a vector of length nrow(data) to be used as INDEX param of tapply
    #   weight.by:     name of a column in data or NULL. 
    #                  If weight.by is not NULL then the weighted mean of 
    #                  all other numeric columns is computed and the 
    #                  arithmetic mean of weight.by is computed. 
    #                  If weight.by is NULL then the arithmetic mean of all     
    #                  numeric columns is computed.

    if(is.character(aggregate.by))
    {
        # tapply (used below) converts aggregate.by to a factor. 
        # If aggregate.by is a character, the resulting factor will be sorted 
        # alphabetically. Construct a factor here to retain the existing 
        # ordering.
        aggregate.by <- factor(aggregate.by, levels=unique(aggregate.by))
    }

    # class.behaviour and column.behaviour were arguments in a previous 
    # version of Cheddar.

    #  class.behaviour - a named list of functions. Names should be the names 
    #                    of R classes (character etc). Functions 
    #                    should accept a single vector of values and should 
    #                    return a single value, e.g. JoinUnique etc

    # We take the artihmetic mean of numerics and join together unique values 
    # of everything else, empty strings and NA ignored.
    class.behaviour  <- list(integer = mean, 
                             numeric = mean, 
                             character = .JoinUniqueEmptyRm, 
                             logical = .JoinUniqueEmptyRm)

    #   column.behaviour: NULL or a named list of functions. Overrides 
    #                     class.behaviour.
    column.behaviour <- NULL


    # Compute weighting if necessary
    if(!is.null(weight.by) && weight.by %in% colnames(data))
    {
        stopifnot(class(data[,weight.by]) %in% c('integer', 'numeric'))
        weighting <- data[,weight.by]

        # Compute the unweighted arithmetic mean of the weight.by column
        if(is.null(column.behaviour))
        {
            column.behaviour <- list()
        }
        column.behaviour[[weight.by]] <- mean
    }

    new.data <- lapply(colnames(data), function(n) 
    {
        if(!is.null(weight.by) && 
           !n %in% names(column.behaviour) && 
           class(data[,n]) %in% c('integer', 'numeric'))
        {
            # A numeric or integer column for which there is no 
            # column-specific override
            vals <- sapply(split(data[,c(n,weight.by)], aggregate.by), 
                           function(rows)
                           {
                                weighted.mean(rows[,1], rows[,2])
                           })
        }
        else
        {
            fn <- column.behaviour[[n]]
            if(is.null(fn)) fn <- class.behaviour[[class(data[,n])]]
            if(is.null(fn))
            {
                stop(paste('No function to aggregate node property [', n, 
                           '] of class [', class(data[,n]), ']', sep=''))
                     
            }

            vals <- tryCatch(as.vector(tapply(data[,n], aggregate.by, fn)), 
                             error=function(e) e)
        }

        if('simpleError' %in% class(vals))
        {
            # Raising this error inside tryCatch() results in an ugly 
            # traceback.
            stop(paste('Error aggregating [', n, ']: ', vals$message, sep=''))
        }
        else
        {
            stopifnot(length(vals)==length(unique(aggregate.by)))

            col <- data.frame(vals, stringsAsFactors=FALSE)
            colnames(col) <- n

            return (col)
        }
    })
    new.data <- do.call('cbind', new.data)
    return (new.data)
}

