# Stuff that is useful and doesn't belong anywhere else
FormatLM <- function(model, slope.95.ci=FALSE, ci.plus.minus.style=FALSE, 
                     r=FALSE, r.squared=TRUE, model.var.names=TRUE, 
                     dp=2)
{
    if(is.null(model))    return('')

    # An expression representation of a linear model equation

    # See "Automation of Mathematical Annotation in Plots", Uwe Liggs, R News, 
    # Vol. 2/3, December 2002

    fmt <- paste('%.', dp, 'f', sep='')

    variables <- as.character(attr(terms(model), 'variables'))
    y.name <- variables[1+attr(terms(model), 'response')]
    x.name <- names(model$coefficients)[2]

    # Slope confidence interavals
    confint <- ''
    if(slope.95.ci && !ci.plus.minus.style)
    {
        ci <- sprintf(fmt, confint(model, level=0.95)[x.name,])
        confint <- paste("~'(95% CI ", ci[1], ",", ci[2], ")'", sep='')
        confint <- parse(text=confint)[[1]]
    }
    else if(slope.95.ci && ci.plus.minus.style)
    {
        ci <- sprintf(fmt,diff(as.vector(confint(model, x.name, level=0.95)))/2)
        n <- nrow(model$model)
        confint <- paste("'' %+-%", ci, "~'(95% CI, n=", n, ")'", sep='')
        confint <- parse(text=confint)[[1]]
    }

    # Correlation coefficient
    r.val = ''
    if(r)
    {
        r.val <- sign(model$coefficients[x.name]) * summary(model)$r.squared^0.5
        r.val <- paste("','~r == '", sprintf(fmt, r.val), "'", sep='')
        r.val <- parse(text=r.val)[[1]]
    }

    # Coefficient of determination
    r2 <- ''
    if(r.squared)
    {
        r2 <- paste("','~r^2 == '", sprintf(fmt, summary(model)$r.squared), "'",
                    sep='')
        r2 <- parse(text=r2)[[1]]
    }

    # Coefficients
    co <- model$coefficients
    return (eval(substitute(expression(y == a~sign~b * x * confint * r * r2), 
                            list(y=ifelse(model.var.names, y.name, 'y'),
                                 a=sprintf(fmt, co[1]), 
                                 sign=ifelse(co[2]<0, '-', '+'),
                                 b=sprintf(fmt, abs(co[2])), 
                                 x=ifelse(model.var.names, x.name, 'x'),
                                 confint=confint, 
                                 r=r.val,
                                 r2=r2))))
}

.StripWhitespace <- function(v)
{
    # Returns v with leading and trailing whitespace removed
    return (gsub('^\\s+|\\s+$', '', v, perl=TRUE))
}

.SimpleRBindFill <- function(..., fill.with=NA)
{
    # A quick and dirty wrapper around rbind.data.frame() that fills 
    # missing columns with NA.
    allargs <- list(...)
    stopifnot(all(sapply(allargs, is.data.frame)))
    stopifnot(all(sapply(allargs, nrow)>0))

    allcnames <- unique(unlist(sapply(allargs, colnames)))

    # Add columns as required
    for(d in 1:length(allargs))
    {
        missing <- setdiff(allcnames, colnames(allargs[[d]]))
        allargs[[d]][,missing] <- fill.with
    }

    res <- do.call('rbind.data.frame', allargs)

    return (res)
}

PredationMatrixToLinks <- function(pm, link.property=NULL)
{
    # Returns a data.frame containing columns resource and consumer. 
    # pm should be a matrix or data.frame containing. Non-zero values indicate 
    # a trophic link between a consumer (column) and a resource (row).
    # If property is non-NULL then 

    if(2!=length(dim(pm)))
    {
        stop('pm is not a matrix')
    }

    # Check names
    if(is.null(colnames(pm)) || is.null(rownames(pm)))
    {
        stop('pm is missing either rownames or colnames')
    }

    resource <- which(pm!=0 & !is.na(pm)) %% nrow(pm)
    consumer <- 1+(which(pm!=0 & !is.na(pm)) %/% nrow(pm))

    # Fix the last row
    last.row <- which(resource==0)
    consumer[last.row] <- consumer[last.row]-1
    resource[last.row] <- nrow(pm)

    df <- data.frame(resource=rownames(pm)[resource], 
                     consumer=colnames(pm)[consumer], 
                     stringsAsFactors=FALSE)

    if(!is.null(link.property))
    {
        df <- cbind(df, pm[pm!=0 & !is.na(pm)])
        colnames(df) <- c('resource', 'consumer', link.property)
    }
    return (df)
}

.CallAssemblePropertiesFunction <- function(f, expected.size, args, 
                                            user.args=NULL)
{
    # Private helper for assembling computed properties
    # f - name of function
    # expected.size - if > 1, f should return either vector of length 
    #                 expected.size or a matrix / data.frame with 
    #                 expected.size rows. if == 1, f should return a vector 
    #                 of at least length 1
    # args - list of arguments
    # user.args - list of optional extra args provided by the user
    
    bad <- names(user.args) %in% names(args)
    if(any(bad))
    {
        stop(paste('The names [', 
                   paste(names(user.args)[bad], collapse=','), 
                   '] are not allowed in the call to', 
                   '[', f, ']', sep=''))
    }

    v <- do.call(f, args=c(user.args, args))

    # Convert results where length(v)>1 into a data.frame
    if(1==expected.size && is.null(dim(v)) && length(v)>1)
    {
        v <- t(v)
    }

    if(!is.null(dim(v)))
    {
        colnames(v)[colnames(v)==''] <- .UnnamedString()
    }

    # Check sizes
    v.length <- ifelse(is.null(dim(v)), length(v), nrow(v))
    if(expected.size!=v.length)
    {
        stop(paste('The function [', f, '] returned a vector, data.frame ', 
                   'or matrix of the incorrect size', sep=''))
    }

    return (as.data.frame(v, stringsAsFactors=FALSE, check.names=FALSE))
}

.AssembleProperties <- function(first.class, properties, ...)
{
    # Returns a data.frame of nrow(first.class) rows. 
    # first.class - a data.frame of the first-class properties 
    # properties - either a vector of characters or a list. 
    # If a vector of character, properties should contain the names of 
    # first-class properties or the names of functions that take a single 
    # community and return either a vector of length nrow(first.class) or a 
    # matrix of nrow nrow(first.class). 

    # If a list, elements should be either characters or lists. 
    # If a character, elements are treated as above. If a list, 
    # the first element be the name of a function and the following elements 
    # should be arguments to the function given in the first element.

    # Column names of the returned data.frame are taken from names(properties), 
    # if set, for first.class properties and functions that return vectors.

    # Functions should return either vectors of length expected.size or 
    # matrices / data.frames with expected.size rows.
    expected.size <- nrow(first.class)

    res <- as.data.frame(matrix(NA, ncol=0, nrow=expected.size), 
                         stringsAsFactors=FALSE)

    for(index in 1:length(properties))
    {
        if(is.list(properties)) p <- properties[[index]]
        else                    p <- properties[index]

        if(is.list(p))
        {
            # A function + args
            if(length(p)>1 && is.character(p[[1]]))
            {
                f <- p[[1]]
                args <- p[2:length(p)]
                v <- .CallAssemblePropertiesFunction(f, expected.size, 
                                                     list(...), args)
            }
            else
            {
                stop(paste('Badly formed computed property given in [', index, 
                           ']', sep=''))
            }
        }
        else if(p %in% colnames(first.class))
        {
            # A first-class property
            v <- first.class[,p,drop=FALSE]
        }
        else if(tryCatch(is.function(eval(parse(text=p))),
                         error=function(e) FALSE) && 'category'!=p)
        {
            # A function

            # WARNING: Nasty hack in the above test
            # category is a defunct R function. 'category' is a 
            # commonly-requested node property. If this community does not 
            # have a first-class property called 'category', we avoid 
            # the nasty error message raised by trying to call the defunct 
            # category function.
            v <- .CallAssemblePropertiesFunction(p, expected.size, list(...))
        }
        else
        {
            # A bad name - use NA
            v <- as.data.frame(rep(NA, expected.size))
        }

        # Set results name(s)
        # Use the function or property name as the column name
        if(is.list(p))
        {
            v.name <- p[1]
        }
        else
        {
            v.name <- p
        }

        if(!is.null(names(properties)) && ''!=names(properties)[index])
        {
            v.name <- names(properties)[index]
            if(1<ncol(v))  colnames(v) <- paste(v.name, colnames(v), sep='.')
        }

        if(1==ncol(v)) colnames(v) <- v.name

        res <- cbind(res, v)
    }
    rownames(res) <- rownames(first.class)
    return (res)
}

.CAdjacencyList <- function(community, values)
{
    # Returns an adjacency list suitable for passing to C functions. 
    # values should be a list like those returned by ConsumersByNode() and 
    # ResourcesByNode().

    # The returned matrix will have NumberOfNodes() rows in which first col is 
    # number of consumers (resources) of that node and subsequent cols are ids 
    # of consumers (resource). Elements are 0-indexed node indices.
    # Returns an adjacency list. 
    # values should be a list like those returned by ConsumersByNode() and 
    # ResourcesByNode().

    # The returned matrix will have NumberOfNodes() rows in which first col is 
    # number of consumers (resources) of that node and subsequent cols are ids 
    # of consumers (resource).
    values <- lapply(values, NodeNameIndices, community=community)
    values <- lapply(values, unname)
    a <- matrix(NA, nrow=NumberOfNodes(community), 
                    ncol=1+max(sapply(values, length)))
    for(i in 1:length(values))
    {
        n <- length(values[[i]])
        a[i,1] <- n
        if(n>0)
        {
            a[i, 2:(1+n)] <- values[[i]] - 1
        }
    }

    rownames(a) <- unname(NP(community, 'node'))
    return (a)
}

.StringStartsWith <- function(string, starts.with)
{
    return (starts.with==substr(string, 1, nchar(starts.with)))
}

.WriteForRob <- function(community, path)
{
    # Writes community in a format required by Rob's code that computes chains
    if(!is.Community(community)) stop('Not a Community')

    x <- .CAdjacencyList(community, ConsumersByNode(community))
    # First col is the number of consumers. Rob's code requires 1-indexed node numbers.
    x[,1] <- 1:nrow(x)

    # Other columns are 0-indexed node number. Rob's code requires 1-indexed node numbers.
    x[,2:ncol(x)] <- x[,2:ncol(x)]+1

    if(missing(path))
    {
        path <- paste(gsub('[ \\.]', '_', CPS(community)$title), '_input.txt', sep='')
    }

    write.table(x, file=path, na='', row.names=FALSE, col.names=FALSE, sep=' ', 
                quote=FALSE)
}
