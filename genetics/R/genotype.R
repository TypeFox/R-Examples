# $Id: genotype.R 1337 2008-04-30 00:54:56Z warnes $

genotype  <- function(a1, a2=NULL, alleles=NULL, sep="/",
                      remove.spaces=TRUE,
                      reorder=c("yes", "no", "default", "ascii", "freq"),
                      allow.partial.missing=FALSE,
                      locus=NULL, genotypeOrder=NULL)
{
    if(missing(reorder))
      reorder  <- "freq"
    else
      reorder <- match.arg(reorder)

    if(is.genotype(a1)){
        a1  <-  as.character(a1)
        ## ignore a2
        a2 <- NULL
    }
    else
      {
        a1.d <- dim(a1)
        a1 <- as.character(a1)
        dim(a1) <- a1.d
        a1[is.na(a1)] <- ""   # necessary because of bug in grep & friends,
                              # will be fixed in 1.7.1
      }

    if(!is.null(a2))
      {
        a2.d <- dim(a2)
        a2 <- as.character(a2)
        dim(a2) <- a2.d
        a1[is.na(a1)] <- ""  # necessary because of bugs in grep & friends
                             # will be fixed in 1.7.1
      }

    if(remove.spaces)
    {
        a1dim <- dim(a1)
        a1  <-  gsub("[ \t]", "", a1)
        dim(a1) <- a1dim
        if(!is.null(a2))
            a2  <-  gsub("[ \t]", "", a2)
    }

    if(!is.null(dim(a1)) && ncol(a1) > 1)
        parts <- a1[,1:2]
    else if(!is.null(a2))
        parts  <- cbind(a1,a2)
    else
      {
        # if sep is empty, assume allele names are single characters
        # pasted together
        if(sep=="")
          sep  <- 1

        # Based on the value of sep, reformat into our standard
        # name-slash-name format
        if (is.character(sep) )
          {
            part.list   <- strsplit(a1,sep)
            part.list[ sapply(part.list, length)==0] <- NA

            ## Handle missing / empty values correctly. 
            ## Without this, empty elements are silently dropped
            ## and/or cause errors

            # only first field was given
            half.empties  <- lapply(part.list, length)==1
            part.list[half.empties]  <-  lapply(part.list[half.empties],c,NA)

            # neither field was given
            empties  <- is.na(a1) | lapply(part.list, length)==0
            part.list[empties]  <- list(c(NA,NA))

            parts <- matrix(unlist(part.list),ncol=2,byrow=TRUE)

          }
        else if (is.numeric(sep))
          parts  <- cbind( substring(a1,1,sep), substring(a1,sep+1,9999))
        else
          stop(paste("I don't know how to handle sep=",sep))
      }

    mode(parts) <- "character"  # needed for bare NA's o

    # convert entirely whitespace alleles to NAs
    temp  <- grep("^[ \t]*$", parts)
    parts[temp]  <-  NA

    #parts[parts=="NA"]  <-  NA

    if(missing(alleles) || is.null(alleles))
      alleles <- unique(c(na.omit(parts)))
    else
      {
        which.alleles  <- !(parts %in% alleles)
        ## Skipping NA's
        which.alleles <- which.alleles & !is.na(parts)
        if(any(which.alleles))
          {
            warning("Found data values not matching specified alleles. ",
                    "Converting to NA.")
            parts[which.alleles] <- NA
          }
      }

    if(!allow.partial.missing)
      parts[is.na(parts[,1]) | is.na(parts[,2]),]  <- c(NA,NA)

    if(reorder!="no")
    {
        if(reorder=="ascii")
        {
            alleles <-  sort(alleles)
        }
        else if(reorder=="freq")
        {
            ## get reordering of alleles by frequency
            tmp  <- names(rev(sort(table(parts))))
            alleles  <- unique(c(tmp,alleles))
        }

        reorder  <- function( x, alleles)
        {
            tmp <- match( x, alleles )
            x[order(tmp)]
        }

        parts  <- t(apply(parts,1, reorder, alleles))

      }

    tmp  <-  ifelse( is.na(parts[,1]) & is.na(parts[,2]),
                    NA,
                    apply(parts,1,paste,collapse="/") )

    object  <- factor( tmp )

    # force "NA" not to be a factor level
    ll  <- levels(object)  <-  na.omit(levels(object))

    class(object)  <-  c("genotype","factor")
    attr(object,"allele.names")  <- alleles
    attr(object,"allele.map")  <- do.call("rbind", strsplit(ll, "/"))

    genotypeOrder(object) <- genotypeOrder

    if(is.null(locus) || is.locus(locus)  )
      attr(object,"locus")  <- locus
    else
      stop("parameter locus must be of class locus")
    return(object)
  }

is.genotype  <- function(x)
    inherits(x, "genotype")

is.haplotype  <- function(x)
    inherits(x, "haplotype")

###
### Haplotype -- differs only in that order of a1,a2 is considered siginificant
###
haplotype <- function (a1, a2 = NULL, alleles = NULL, sep = "/",
                       remove.spaces = TRUE, reorder = "no",
                       allow.partial.missing = FALSE, locus = NULL,
                       genotypeOrder=NULL) 
{
    retval <- genotype(a1 = a1, a2 = a2, alleles = alleles, sep = sep, 
                       remove.spaces = remove.spaces, reorder = reorder,
                       allow.partial.missing = allow.partial.missing, 
                       locus = locus, genotypeOrder=genotypeOrder)
    class(retval) <- c("haplotype", "genotype", "factor")
    retval
}


as.haplotype  <- function(x,...)
{
    retval <- as.genotype(x,...,reorder="no")
    class(retval)  <- c("haplotype","genotype","factor")
    retval
}

###
### Display by giving values plus list of alleles
###

print.genotype  <-  function(x,...)
  {
    if(!is.null(attr(x,"locus")))
        print(attr(x,"locus"))
    print(as.character(x))
    cat("Alleles:", allele.names(x), "\n" )
    invisible(x)
  }

###
### Conversion Functions
###

as.genotype  <- function (x,...) 
  UseMethod("as.genotype")

# Do we want to do this?
as.genotype.default  <-  function(x,...)
  genotype(x,...)

#  stop("No method to convert this object to a genotype")

# for characters, and factors, just do the standard thing (factors get
# implicitly converted to characters so both have the same effect.
as.genotype.character  <-  function(x,...)
  genotype(x,...)

as.genotype.factor  <-  function(x,...)
  genotype(as.character(x),...)

as.genotype.genotype  <- function(x,...)
  return(x)

as.genotype.haplotype  <- function(x,...)
  return(x)


## genotype.allele.counts give the count of each allele type as a
## matrix.  Collapse back into the form we need

as.genotype.allele.count  <- function(x, alleles=c("A","B"), ...)
  {
    if(!is.matrix(x) & !is.data.frame(x) )
      {
        x  <- cbind(x, 2-x)
        colnames(x)  <- alleles
      }

    if(any(x > 2, na.rm=TRUE) || any( x < 0, na.rm=TRUE ) )
      stop("Allele counts must be in {0,1,2}")

    allele.names  <-  colnames(x)
    tmp  <-  apply(x, 1, function(y)
                    rep( colnames(x), ifelse(is.na(y), 0, y) ))

    if(!is.matrix(tmp))
      retval  <-  genotype(sapply(tmp,paste,collapse="/"), alleles=alleles, ...)
    else
      retval  <- genotype(a1=tmp[1,], a2=tmp[2,], ... )
    return(retval)
  }

allele.count.2.genotype  <-  function(...)
  as.genotype.allele.count(...)



as.genotype.table <- function(x, alleles, ...)
  {
    #if(missing(alleles)) alleles <- unique(unlist(dimnames(x)))
    tmp <- outer( rownames(x), colnames(x), paste, sep="/")
    retval <- genotype( rep(tmp,x), alleles=alleles )
    retval
  }


###
### Equality test for genotype, assumes allele order is _not_ significant
###
"==.genotype"  <-  function(x,y)
  {
    if(!is.genotype(y))
      y <- as.genotype(y)

    x.a1  <- allele(x,1)
    x.a2  <- allele(x,2)

    y.a1  <- allele(y,1)
    y.a2  <- allele(y,2)

    return( (x.a1==y.a1 & x.a2==y.a2) | (x.a1==y.a2 & x.a2==y.a1) )
  }

###
### Equality test for haplotype, assumes allele order _is_ significant
###
"==.haplotype"  <-  function(x,y)
  {
    if(!is.genotype(y))
      y <- as.haplotype(y)

    x.a1  <- allele(x,1)
    x.a2  <- allele(x,2)

    y.a1  <- allele(y,1)
    y.a2  <- allele(y,2)

    return( x.a1==y.a1 & x.a2==y.a2 )
  }

###
### is.element i.e. %in%
###

"%in%" <- function(x, table)
  UseMethod("%in%")

## Get default method for %in% from base package
"%in%.default" <- get("%in%", pos="package:base")

"%in%.genotype" <- function(x, table)
{
  xA1  <- allele(x, 1)
  xA2  <- allele(x, 2)

  x1 <- paste(xA1, xA2, sep="/")
  x2 <- paste(xA2, xA1, sep="/")

  ## Return
  ((x1 %in% table) | (x2 %in% table))
}

"%in%.haplotype" <- function(x, table)
  as.character(x) %in% as.character(table)

###
### Extract the first and/or second allele.
###
### By default, return a 2 column matrix containing both alleles
###

#allele  <- function (x,...) 
#  UseMethod("allele")


#allele.genotype  <-  function(x, which=c(1,2) )
allele  <-  function(x, which=c(1,2) )
  {
    alleles.x  <- attr(x,"allele.map")
    retval  <- alleles.x[as.integer(x),which]
    attr(retval,"locus")  <- attr(x,"locus")
    attr(retval,"which")  <- which
    attr(retval,"allele.names")  <- allele.names(x)    
    #class(retval)  <- c("allele.genotype", class(retval))
    return( retval)
  }

as.factor  <- function(x, ...)
  UseMethod("as.factor")

as.factor.default  <- get("as.factor",pos="package:base")
formals(as.factor.default) <- c(formals(as.factor.default),alist(...= ))

as.factor.genotype <- function(x, ...)
  {
    attr(x,"class") <- "factor"
    attr(x,"allele.names") <- NULL
    attr(x,"allele.map") <- NULL
    attr(x,"locus") <- NULL
    attr(x,"genotypeOrder") <- NULL
    x
  }

as.factor.allele.genotype  <-  function(x,...)
  factor(x,levels=allele.names(x))

print.allele.genotype  <- function(x,...)
  {
    if(!is.null(attr(x,"locus")))
      print(attr(x,"locus"))
    cat("Allele(s):", attr(x,"which"), "\n")
    attr(x, "which")  <-  attr(x, "class") <- attr(x,"locus") <- attr(x,"allele.names")  <- NULL
    NextMethod("print",x)
  }


###
### Obtain the count of the number of copies of alleles for each individual
###
### By default, return a matrix containing the counts for all possible allele values.
###

#allele.count  <- function (x,...) 
#  UseMethod("allele.count")

#allele.count.default <- function (x, ... )
#  {
#    x <- as.genotype(x)
#    allele.count(x, ...)
#  }

#allele.count.genotype  <- function(x, allele.name=allele.names(x),

allele.count  <- function(x, allele.name=allele.names(x),
                          any=!missing(allele.name), na.rm=FALSE)
{
  if(!missing(allele.name) && length(allele.name)==1)
    {
      a.1  <- allele(x,1)
      a.2  <- allele(x,2)

      retval  <- ifelse(is.na(a.1) | is.na(a.2),
                        ifelse(na.rm, 0, NA),
                        (a.1==allele.name) + (a.2==allele.name) )
#      class(retval)  <- "allele.count"
      attr(retval,"allele") <- allele.name
      attr(retval,"locus")  <- attr(x,"locus")
      return(retval)
    }
  else
    {
      retval  <- sapply( allele.name, function(y) allele.count(x,y))
      if(any==TRUE && is.matrix(retval)  )
      retval  <- apply(retval,1,sum,na.rm=na.rm)
      if(na.rm) retval[is.na(retval)]  <- 0
#      class(retval)  <- "allele.count"
      attr(retval,"locus")  <- attr(x,"locus")
      return(retval)
    }

}


#print.allele.count  <- function(x,...)
#  { 
#    if(!is.null(attr(x,"locus")))
#        print(attr(x,"locus"))
#    
#    if(is.null(attr(x,"allele")))
#      cat("Allele Counts:\n")
#    else
#      cat("Allele Count (", attr(x,"allele"), " allele):\n", sep="")
#    val  <- x
#    attr(val,"class")  <- NULL
#    attr(val,"allele")  <- NULL
#    print(val)
#    invisible(x)
#  }

###
### Check for the presence of alleles for each individual
###
### By default, return a matrix containing indicators for all possible
### allele values except the last.
###
#
#allele.ind  <-  function(x,allele)
#  {
##    if(missing(allele))
##      stop("Alleles to test must be specified")
##    if(length(allele)==1)
#      retval  <- allele.count(x,allele) > 0
##    else
##      retval  <- apply(allele.count(x,allele) ,1,sum) > 0
#
#      if(missing(allele))
#          allele  <-  colnames(retval)
#      attr(retval,"allele")  <- allele
#      attr(retval,"locus")  <- attr(x,"locus")
#      class(retval)  <-  "allele.ind"
#      return(retval)
#  }

#print.allele.ind  <- function(x,...)
#  {
#    if(!is.null(attr(x,"locus")))
#      print(attr(x,"locus"))
#    
#    cat("Indicator(s) for allele(s):", attr(x,"allele"), "\n")
#    attr(x,"locus")  <-  attr(x,"class")  <- attr(x,"allele")  <-  NULL
#    NextMethod("print",x)
#  }

###
### Methods for creating subsets based on a genotype
###

homozygote  <- function (x,allele.name,...) 
  UseMethod("homozygote")

homozygote.genotype  <-  function(x,allele.name,...)
  {
    a1  <- allele(x,1)
    a2  <- allele(x,2)
    if(missing(allele.name))
      retval  <- ifelse( is.na(a1) | is.na(a2), NA, a1==a2 )
    else
      retval  <- ifelse( is.na(a1) | is.na(a2), NA,
                         a1==allele.name & a2==allele.name )
    attr(retval,"locus")  <-  attr(x,"locus")
#    class(retval)  <-  "homozygote"
    return(retval)
  }

#print.homozygote  <- function(x,...)
#  {
#    if(!is.null(attr(x,"locus")))
#      print(attr(x,"locus"))
#    
#    cat("Homozygote Indicators:\n")
#    attr(x,"locus")  <-  attr(x,"class")  <- attr(x,"allele")  <-  NULL
#    NextMethod("print",x)
#  }


heterozygote  <- function (x,allele.name,...) 
  UseMethod("heterozygote")

heterozygote.genotype  <-  function(x,allele.name,...)
  {
  {
    a1  <- allele(x,1)
    a2  <- allele(x,2)
    if(missing(allele.name))
      retval  <- ifelse( is.na(a1) | is.na(a2), NA, !a1==a2 )
    else
      retval  <- ((a1 %in% allele.name) | (a2 %in% allele.name)) &
                 (a1 != a2)
    attr(retval,"locus")  <-  attr(x,"locus")
#    class(retval)  <-  "homozygote"
    return(retval)
  }
  }

#print.heterozygote  <- function(x,...)
#  {
#    if(!is.null(attr(x,"locus")))
#      print(attr(x,"locus"))
#    
#    cat("Heterozygote Indicators:\n")
#    attr(x,"locus")  <-  attr(x,"class")  <- attr(x,"allele")  <-  NULL
#    NextMethod("print",x)
#  }

carrier <- function (x,allele.name,...) 
  UseMethod("carrier")

carrier.genotype  <-  function(x, allele.name=allele.names(x),
                                   any=!missing(allele.name), na.rm=FALSE, ...)
{
  retval  <- allele.count(x,allele.name=allele.name,any=any,na.rm=na.rm) > 0

  attr(retval,"allele")  <- retval["allele"]
  attr(retval,"locus")  <-  attr(x,"locus")
#  class(retval)  <- "carrier"
  return(retval)
}


#print.carrier  <- function(x,...)
#  {
#    if(!is.null(attr(x,"locus")))
#      print(attr(x,"locus"))
#    
#    cat("Carrier Indicator(s) for allele(s):", attr(x,"allele"), "\n")
#    attr(x,"locus")  <-  attr(x,"class")  <- attr(x,"allele")  <-  NULL
#    NextMethod("print",unclass(x))
#  }


###
###
###

allele.names<- function(x)
  {
    retval  <- attr(x,"allele.names")
    if(is.null(retval))
      retval  <- x$allele.names
    return(retval)
  }

###
### Subset method
###

"[.genotype"  <-  function(x, i, drop=FALSE)
  {
    allelesOld <- attr(x, "allele.names")
    retval  <- NextMethod("[")

    # force "NA" not to be a factor level
    ll  <- levels(retval)  <-  na.omit(levels(retval))

    class(retval)  <-  c("genotype","factor")

    if(drop) {
      alleles <- unique( unlist(strsplit(ll, "/") ) )
    } else {
      alleles <- attr(x, "allele.names")
    }
    attr(retval,"allele.names")  <- alleles
    attr(retval,"allele.map")  <- do.call("rbind", strsplit(ll, "/"))
    attr(retval,"locus")  <- attr(x,"locus")
    attr(retval,"label")  <-  attr(x,"label")
    goCur <- attr(x, "genotypeOrder")
    if(drop) {
      ## Removing genotype names having dropped alleles
      allelesOld <- allelesOld[!(allelesOld %in% alleles)]
      tmp <- allele(as.haplotype(goCur))
      test <- tmp %in% allelesOld
      test <- rowSums(matrix(test, ncol=2)) > 0
      attr(retval, "genotypeOrder") <- goCur[!test]
    } else {
      attr(retval, "genotypeOrder") <- goCur
    }
    return(retval)
  }

"[.haplotype"  <-  function(x, i, drop=FALSE)
  {
    retval  <- NextMethod("[")
    class(retval) <- c("haplotype","genotype","factor")
    retval
  }

###
### Subset Assigment method
###

"[<-.genotype"  <-  function(x, i, value)
  {
    ## Special case for insertion of NA and "" values
    if(all( is.na(value) | as.character(value)<="" ) )
      {
        x.class <- class(x)
        x <- unclass(x)
        x[i] <- NA
        class(x) <- x.class
        return(x)
      }
      
    if(!is.genotype(value))
      {
        value <- genotype(value)
      }

    lx <- levels(x)
    lv <- levels(value)
    ax <- allele.names(x)
    av <- allele.names(value)

    m  <- is.na(match(av,ax) )
    if( any( m  )  )
       warning(paste("Adding new allele name:", av[m], "\n"))

    la <- unique(c(lx,lv))
    aa <- unique(c(ax,av))

    cx <- class(x)
    nas <- is.na(x)

    data  <-  match(levels(value)[value],la)

    class(x) <- NULL
    x[i] <- data
    attr(x, "levels") <- la
    map  <- attr(x, "allele.map")  <- do.call("rbind", strsplit(la, "/"))
    attr(x, "allele.names")  <- aa
    goCur <- attr(x, "genotypeOrder")
    goAll <- expectedGenotypes(alleles=aa, haplotype=TRUE)
    attr(x, "genotypeOrder") <- c(goCur, goAll[!(goAll %in% goCur)])
    class(x) <- cx
    x
  }

"[<-.haplotype"  <-  function(x, i, value)
  {
    if(!is.haplotype(value))
      stop("Assigned value must be of class haplotype.")
    NextMethod("[<-")
  }

nallele <- function(x)
  length(allele.names(x))

genotypeOrder <- function(x)
  attr(x, "genotypeOrder")

"genotypeOrder<-" <- function(x, value)
{
  if(!is.genotype(x)) stop("'x' must be of a genotype class")
  alleles <- allele.names(x)

  goAll <- expectedGenotypes(alleles=alleles, haplotype=TRUE)
  goDef <- unique(sort(as.character(x)))

  if(is.null(value)) {
    attr(x, "genotypeOrder") <- goAll
  } else {
    value <- unique(value)

    ## Stop msg says all
    parts <- strsplit(x=value, split="/")
    parts <- sapply(parts, c)
    test <- !(parts %in% alleles)
    if(any(test))
      stop("adding genotype names with alleles that are not in the data")

    ## Any genotypes in the data that are not in value?
    test <- !(goDef %in% value)
    if(any(test)) {

      ## These values are in all possible genotypes/haplotypes
      testDefinAll <- goDef[test] %in% goAll
      ## but not in value
      testDefinAllnotVa <- !(goDef[testDefinAll] %in% value)
      goPos <- goDef[testDefinAllnotVa]

      ## We could simply add goPos to value now. However, A/B in goPos
      ## should also match B/A, since genotype() allows reordering of
      ## original data and additionally we want this to work also for
      ## haplotype.

      ## Extend value first. We do not do this before, since one
      ## might not necessarily like to have B/A together with A/B in
      ## first place.
      value <- .genotype2Haplotype(x=value)
      ## Remove heterozygos matches in goPos
      test <- !(goPos %in% value)
      goPos <- goPos[test]

      ## Add goPos to the end of value
      if(any(test)) value <- c(value, goPos)

      ## If there are still some values in all, but not in value
      ## now, we just add them at the end
      testGOnotAll <- !(goAll %in% value)
      if(any(testGOnotAll)) value <- c(value, goAll[testGOnotAll])
    } else {
      value <- .genotype2Haplotype(x=value)
    }
    attr(x, "genotypeOrder") <- value
  }
  x
}

.genotype2Haplotype <- function(x)
{
  ## Internal function
  ##
  ## Returns a character vector of possible haplotypes for given genotypes
  ## in such a way that for say c("A/A", "A/B", "B/B") you get c("A/A",
  ## "A/B", "B/A", "B/B") i.e. "B/A" comes directly after "A/B"!
  ##
  ## x - character, vector of genotype values in form allele1/allele2
  ##
  ## Details
  ## Unique values of x are taken i.e. first occurrence prevails
  ##
  ## Example
  ## .genotype2Haplotype(c("A/A", "A/B", "B/B"))
  ## "A/A" "A/B" "B/A" "B/B"
  ## .genotype2Haplotype(c("B/B", "A/B", "A/A"))
  ## "B/B" "A/B" "B/A" "A/A"

  x <- unique(x)
  N <- length(x)
  parts <- .genotype2Allele(x=x)
  parts <- rbind(parts, parts[, 2:1])
  ind <- rep(1:N, each=2) + c(0, N)
  parts <- parts[ind, ]
  parts <- unique(paste(parts[, 1], parts[, 2], sep="/"))
  parts
}

.genotype2Allele <- function(x)
{
  ## Internal function
  ##
  ## Returns a matrix of alleles from a character vector of genotype names
  ##
  ## x - character, vector of genotype values in form allele1/allele2
  ##
  ## Details:
  ## Coercing to character is done for x.
  ##
  ## Example
  ## .genotype2Allele(c("A/A", "A/B", "B/B"))
  ##      [,1] [,2]
  ## [1,] "A"  "A"
  ## [2,] "A"  "B"
  ## [3,] "B"  "B"

  parts <- strsplit(x=as.character(x), split="/")
  parts <- t(sapply(parts, c))
  parts
}
