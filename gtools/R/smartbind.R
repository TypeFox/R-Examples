##
## Function to do rbind of data frames quickly, even if the columns don't match
##

smartbind <- function(..., fill=NA, sep=':', verbose=FALSE)
  {
    data <- list(...)
    if(is.null(names(data)))
      names(data) <- as.character(1:length(data))
    data <- lapply(data,
                   function(x)
                   if(is.matrix(x) || is.data.frame(x))
                     x
                   else
                     data.frame(as.list(x), check.names=FALSE)
                   )

    #retval <- new.env()
    retval <- list()
    rowLens <- unlist(lapply(data, nrow))
    nrows <- sum(rowLens)

    rowNameList <- unlist(lapply( names(data),
                                 function(x)
                                   if(rowLens[x]<=1) x
                                   else paste(x, seq(1,rowLens[x]),sep=sep))
                          )

    colClassList     <- vector(mode="list", length=length(data))
    factorColumnList <- vector(mode="list", length=length(data))
    factorLevelList  <- vector(mode="list", length=length(data))


    start      <- 1
    blockIndex <- 1
    for(block in data)
      {
        colClassList    [[blockIndex]] <- list()
        factorColumnList[[blockIndex]] <- character(length=0)
        factorLevelList [[blockIndex]] <- list()

        if(verbose) print(head(block))
        end <- start+nrow(block)-1
        for(col in colnames(block))
          {
            classVec <- class(block[,col])

            ## store class and factor level information for later use
            colClassList[[blockIndex]][[col]] <- classVec
            if("factor" %in% classVec)
              {

                factorColumnList[[blockIndex]] <-
                  c(factorColumnList[[blockIndex]], col)

                factorLevelList[[blockIndex]][[col]] <-
                  levels(block[,col])
            }

            if( !(col %in% names(retval)))
              {
                if(verbose) cat("Start:", start,
                                "  End:", end,
                                "  Column:", col,
                                "\n", sep="")

                if ("factor" %in% classVec)
                  {
                    newclass <- "character"
                  }
                else
                  newclass <- classVec[1]

                ## Coerce everything that isn't a native type to character
                if(! (newclass %in% c("logical", "integer", "numeric",
                                     "complex", "character", "raw") ))
                    {
                        newclass <- "character"
                        warning("Converting non-atomic type column '", col,
                                "' to type character.")
                    }

                retval[[col]] <- as.vector(rep(fill,nrows), mode=newclass)
              }

            mode <- class(retval[[col]])
            if(mode=="character")
                vals <- as.character(block[,col])
            else
                vals <- block[,col]

            retval[[col]][start:end] <- as.vector(vals, mode=mode)
          }
        start <- end+1
        blockIndex <- blockIndex+1
      }

    all.equal.or.null <- function(x,y,...)
      {
        if(is.null(x) || is.null(y) )
          return(TRUE)
        else
          return(all.equal(x,y,...))
      }

    ## Handle factors, merging levels
    for( col in unique(unlist(factorColumnList)) )
      {
        ## Ensure column classes match across blocks
        colClasses <- lapply(colClassList, function(x) x[[col]])
        firstNotNull <- which(!sapply(colClasses, is.null))[1]
        allSameOrNull <- all(sapply(colClasses[-firstNotNull],
                              function(x) isTRUE(all.equal.or.null(colClasses[[firstNotNull]], x))
                              )
                       )

        if(allSameOrNull)
          {
            # grab the first *non-NULL* class information
            colClass <- colClasses[[firstNotNull]]
          }
        else
          {
            warning("Column class mismatch for '", col, "'. ",
                    "Converting column to class 'character'.")
            next()
          }


        ## check if factor levels are all the same
        colLevels <- lapply(factorLevelList, function(x) x[[col]])
        firstNotNull <- which(!sapply(colLevels, is.null))[1]
        allSameOrNull <- all(sapply(colLevels[-firstNotNull],
                                    function(x) isTRUE(all.equal.or.null(colLevels[[firstNotNull]], x))
                                    )
                             )


        if(allSameOrNull)
          {
            if("ordered" %in% colClass)
              retval[[col]] <- ordered(retval[[col]], levels=colLevels[[firstNotNull]] )
            else
              retval[[col]] <- factor(retval[[col]], levels=colLevels[[firstNotNull]] )
          }
        else
          {
            ## Check if longest set of levels is a superset of all others,
            ## and use that one
            longestIndex  <- which.max( sapply(colLevels, length) )
            longestLevels <- colLevels[[longestIndex]]
            allSubset <- sapply(colLevels[-longestIndex],
                                function(l) all(l %in% longestLevels)
                                )
            if(allSubset)
              {
                if("ordered" %in% colClass)
                  retval[[col]] <- ordered(retval[[col]], levels=longestLevels )
                else
                  retval[[col]] <- factor(retval[[col]], levels=longestLevels )
              }
            else
              {
                # form superset by appending to longest level set
                levelSuperSet <- unique(c(longestLevels, unlist(colLevels)))
                retval[[col]] <- factor(retval[[col]], levels=levelSuperSet )

                if(length(colClass)>1) # not just plain factor
                   {
                    warning( "column '", col, "' of class ",
                            paste("'", colClass, "'", collapse=":",
                                  sep="'"),
                            " converted to class 'factor'. Check level ordering." )
                  }

              }
          }
      }

    attr(retval,"row.names") <- rowNameList
    class(retval) <- "data.frame"
    return(retval)
  }

