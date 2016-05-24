`expected.segRatio` <-
function(ploidy.level=stop("No ploidy level set"),
                              type.parents=c("heterogeneous","homozygous"))
{
  ## Description: Generates expected proportions for various dosages for
  ## dominant markers in regular autopolyploids
  ##
  ## Argument:
  ## ploidy.level: the number of homologous chromosomes, either as numeric
  ## or as a character string
  ##
  ## Value:
  ## ratio: vector of proportions for each dosage
  ## ploidy.level
  ## ploidy.name
  ##
  ## NB: needs better error checking if ploidy.level is numeric

  type <- match.arg(type.parents)
  
  p.names <- c("Diploid","Tetraploid","Hexaploid","Octaploid",
               "Decaploid","Dodecaploid","Heccaidecaploid","Unknown")
  
  if(length(ploidy.level)>1)
    stop("Error: only one ploidy level allowed")

  if(!(is.numeric(ploidy.level)|is.character(ploidy.level)))
    stop(cat("Error: illegal ploidy level specified\n",
             " 'ploidy.level' must be numeric or one of \n",
             p.names[-length(p.names)],"\n"))
  
  if(is.numeric(ploidy.level)) {
    pl <- ploidy.level
  }

  ## if character, trim ploidy level to first 4 characters and return
  ## numeric ploidy level
  if(is.character(ploidy.level)) {
    pl <- switch(tolower(strtrim(ploidy.level,4)),
                 dipl=2, tetr=4, hexa=6, octa=8, octo=8, decc=10,
                 dode=12, hecc=14, -1)
  }

  if(pl<2) stop(cat("Error: illegal ploidy level specified\n",
                    " 'ploidy.level' must be numeric or one of \n",
                    p.names[-length(p.names)],"\n"))
  
  n.doses <- pl/2

  if (type == "heterogeneous") {
##    r.names <- c("SDRF","DDRF","TDRF","QRDF","5RDF","6RDF","7RDF","8RDF","9RDF")
    r.names <- c("SD","DD","TD","QD","5D","6D","7D","8D","9D")
    ratio <- rep(NA, n.doses)
    names(ratio) <- r.names[1:n.doses]
    for(k in 1:n.doses){ # ripol et al (1999) equation 1, NB: choose(k,0)=1
      ratio[k] <- 1 - choose(pl-k,n.doses)/choose(pl,n.doses)
    }
  } else {
    r.names <- c("SDxSD","SDxDD","DDxDD","DDxTD","TDxTD","TDxQD","QDxQD",
                 "QDx5D","5Dx5D","5Dx6D","6Dx6D","6Dx7D","7Dx7D",
                 "7Dx8D","8Dx8D")
    tot.doses <- pl - 1
    ratio <- rep(NA, tot.doses)
    names(ratio) <- r.names[1:tot.doses]
    for(j in 1:n.doses){ # deduced for octaploids only but guess its OK
      for (k in 1:n.doses){
        ratio[j+k-1] <- 1 -
          choose(pl-j,n.doses)*choose(pl-k,n.doses)/choose(pl,n.doses)^2
      }
    }
  }

  ## if ploidy level even thats OK but otherwise print a warning message
  if ( (pl %% 2) != 0 ) {
    warn <- "Warning: ploidy level not even - results may be unexpected"
    cat(warn,"\n")
    return(list(ratio=ratio, ploidy.level=pl, ploidy.name="Unknown",
                type.parents=type, warning=warn))
  } else {
    return(list(ratio=ratio, ploidy.level=pl, ploidy.name=p.names[n.doses],
                type.parents=type))
  }
}

