genotype<-function(a1, a2=NULL, alleles=NULL, sep="/",
                      remove.spaces=TRUE,
                      reorder=c("yes", "no", "default", "ascii", "freq"),
                      allow.partial.missing=FALSE,
                      locus=NULL)                    
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

    if(!allow.partial.missing)
      parts[is.na(parts[,1]) | is.na(parts[,2]),]  <- c(NA,NA)
    
    if(missing(alleles) || is.null(alleles))
      alleles <- unique(c(na.omit(parts)))
    else
      {
        which.alleles  <- !(parts %in% alleles)
        if(any(which.alleles))
          {
            warning("Found data values not matching specified alleles. ",
                    "Converting to NA.")
            parts[which.alleles] <- NA
          }
      }

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
    if(is.null(locus) || is.locus(locus)  )
      attr(object,"locus")  <- locus
    else
      stop("parameter locus must be of class locus")
    return(object)
  }


