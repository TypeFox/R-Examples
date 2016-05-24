# $Id: makeGenotypes.R 1340 2008-08-20 19:04:32Z warnes $

#
# convert all genotype-compatible variables in a dataframe to genotypes
#
makeGenotypes <- function( data, convert, sep="/", tol=0.5, ...,
                           method=as.genotype )
                          
  {
    data <- as.data.frame(data)
    
    if(missing(convert))
      {
        fun <- function(x) length(unlist(grep(sep, as.character(x) )))
        convert <- sapply( data,  fun )/nrow(data) > tol
      }

    #cat("Convert:");print(convert);

    if(is.list(convert))
      {
        if( !all(sapply(convert,length)==2) )
          stop("When convert is a list, each element must be a 2-vector.")

        namelist <- names(data)
        
        for(pair in convert)
           {
             if(!is.character(pair))
               pair <- namelist[pair]
             # replace first column in pair with new data,
             index <- which(colnames(data)==pair[1])
             data[[ index ]] <- method(data[[ pair[1] ]],
                                       data[[ pair[2] ]], sep=sep,
                                                ... )
             colnames(data)[index] <- paste(pair,collapse=sep)
             data[[ pair[1] ]] <- data[[ pair[2] ]] <- NULL
           }
      }
    else
      {
        if(is.character(convert))
          namelist <- convert
        else
          namelist <- colnames(data)[convert]

        for(col in namelist)
          data[[col]] <- method(data[[col]], sep=sep, ... )
      }

    data
}

makeHaplotypes <- function( data, convert, sep="/", tol=0.9, ... )
  {
    makeGenotypes( data=data, convert=convert, sep=sep, tol=tol,
                  method=as.haplotype, ... )    
  }

