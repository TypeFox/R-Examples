#
# function to read simcoal output and integrate into rmetasim
#
#

#
# Mark's functions to parse arlequin files (modified by KKM)
#

`%upto%` <- function (from, to) 
if (from <= to) from:to else numeric(0)

coal2rmet <- function( file, norm=TRUE){
#  print(paste("enter coal2rmet",file))
  dat <- parse.arlequin( file)
  pops <- unique( dat$pop)
  is.nuc <- !is.factor( dat[,2])
  dat[,-1] <- lapply( dat[,-1,drop=F], as.character)
  if( is.nuc) {
    nr <- nrow( dat)
###    dat <- dat[ rep( 1:nr, each=2),]
###    nl <- (ncol( dat)-1)/2
###    dat[ 2*(1:nr), (1:nl)*2] <- dat[ 2*(1:nr)-1, (1:nl)*2+1]
###    dat <- dat[ ,-2*(1:nl)-1]
### deleted by AES 12/12/07
###  
### all that happens with this object is ultimately the frequencies of alleles
### at the different loci. The above code, in an attempt to 'diploidize' simcoal
### output reduces the number of loci simulated.  Eliminating it cuts the population
### size in half, but that is one half of the simcoal haploid number of individuals.  It
### also uses the correct number of loci
### 
    
  }

  dat[,-1] <- lapply( dat[,-1,drop=F], as.factor)
  n.locs <- ncol( dat)-1
  loctable <- lapply( dat[,-1,drop=F], function( x) {
    m <- t( table( dat[,1], as.integer( x)))
    dimnames( m)[[1]] <- levels( x)
    m
  })

  if( norm)
    loctable <- lapply( loctable, function( x) x / rep( colSums( x), each=nrow( x)))
#  print(paste("exit coal2rmet",file))
  loctable

}

parse.arlequin <- function(file)
{
  dat <- readLines( file)
  data <- NULL
  i.pop <- 0
  while( !is.na( droppo <- grep( 'SampleData *=', dat)[1])) {
    i.pop <- i.pop + 1
    dat <- dat[ -(1:droppo)]
    endo <- grep( '^ *\\} *$', dat)[1]
    dati <- dat[ 1 %upto% (endo-1)]
    blanks <- grep( '^( ||t)*$', dati)
    if( length( blanks))
      dati <- dati[ -blanks]
    is.indiv <- regexpr( '_', dati)>0
    if( is.nuc <- !all( is.indiv))
      dati <- paste( dati[ is.indiv], dati[ !is.indiv], sep=' ')

    fo <- textConnection( dati)
    on.exit( close( fo))
    data.i <- cbind( pop=i.pop, read.table( fo, header=FALSE, row.names=NULL)[ ,-(1:2),drop=FALSE])
    close( fo)
    on.exit()

    if( is.null( data))
      data <- data.i
    else
      data <- rbind( data, data.i)
  }

  if( is.nuc) {
    nl <- (ncol( data)-1)/2
    data <- data[ ,c (1, 1+c( matrix( 1:(2*nl), nrow=2, byrow=T)))]
  }
  data
}

# EG: coal2rmet( 'tossm_0.arp')




#
#
# Allan's function to glue the parse2d data into rmetasim format
# This function will clobber existing individuals and loci in rland object
#
#
#KKM remove 'seqmut' and 'msmut' from arguments list and replace with
#KKM 'mut.rates', which should be an array equal in length to the number of loci.


landscape.coalinput <- function (rland, npp = 200, arlseq = NULL, arlms = NULL, seqsitemut = 1e-07, 
    msmut = 5e-04, mut.rates = NULL) 
{
    if (is.null(arlseq) & is.null(arlms) & is.null(mut.rates)) {
        print("you must specify some type of coalescent input")
        rland
    } else
        {
            if (!is.null(arlseq))
                {
                    clocseq <- coal2rmet(arlseq)
                    ###need some logic to check for monomorphism
                    ###check for all N; if so add a random sequence.  otherwise do nothing
                    clocseq <- lapply(clocseq,function(x)
                                      {
                                          if (length(rownames(x))==1) #this is probably the same check as below
                                          if (nchar(gsub("N","",rownames(x)))==0) #check if all Ns
                                              {
                                                  rownames(x) <- paste0(sample(c("A","G","T","C"),nchar(rownames(x)),
                                                                               replace=T),collapse='')
                                              }
                                          x
                                      })
                }  else clocseq <- NULL
            if (!is.null(clocseq)) {
                clocseq <- clocseq[1]
            }
            if (!is.null(arlms)) {
                clocms <- lapply(arlms, function(x) coal2rmet(x))
                clocms <- do.call(c, clocms)
            } else clocms <- NULL
        cloc <- c(clocseq, clocms)
        if (is.null(mut.rates)) {
            mut.rates = rep(NA, length(cloc))
        }
        for (loc in 1:length(cloc)) {
            states <- as.character(rownames(cloc[[loc]]))
            attr(states, "names") <- NULL
            
            
            ltype <- c(1, 2)[length(grep("T|A|G|C", states[1])) + 1] #1=ssr 2=seq
            freqs <- as.numeric(apply(as.matrix(cloc[[loc]]), 
                1, mean))
            attr(freqs, "names") <- NULL
            if (ltype == 1) {
                if (is.na(mut.rates[loc])) {
                  mut.rates[loc] <- msmut
                }
                rland <- landscape.new.locus(rland, type = ltype, 
                  ploidy = 2, mutationrate = mut.rates[loc], 
                  numalleles = length(states), frequencies = freqs, 
                  states = as.numeric(states), transmission = 0)
            } else {
                if (is.na(mut.rates[loc])) {
                  mut.rates[loc] <- seqsitemut
                }
                rland <- landscape.new.locus(rland, type = ltype, 
                  ploidy = 1, mutationrate = mut.rates[loc], 
                  numalleles = length(states), frequencies = freqs, 
                  states = states, allelesize = nchar(states[1]), 
                  transmission = 1)
            }
        }
        S <- rland$demography$localdem[[1]]$LocalS
        R <- rland$demography$localdem[[1]]$LocalR
        ev <- eigen((R + diag(dim(R)[1])) %*% S)$vectors[, 1]
	##this line here, to standardize the eigenvector--DPG 26feb08
	  ev <- ev/sum(ev)
        if (length(npp) == 1) 
            NperPop <- rep(npp, rland$intparam$habitats)
        else NperPop <- npp
        indlist <- vector("list", length(NperPop))
        for (p in 1:length(NperPop)) {
            im <- matrix(NA, nrow = NperPop[p], ncol = landscape.democol() + 
                sum(landscape.ploidy(rland)))
            colcnt <- landscape.democol() + 1
            for (loc in 1:length(landscape.ploidy(rland))) {
                probs <- cloc[[loc]][, p]
                indices <- sapply(rland$loci[[loc]]$alleles, 
                  function(x) {
                    x$aindex
                  })
                for (al in 1:landscape.ploidy(rland)[loc]) {
                  im[, colcnt] <- sample(indices, NperPop[p], 
                    replace = T, prob = probs)
                  colcnt <- colcnt + 1
                }
            }
            im[, 1] <- sample((0:(rland$intparam$stages - 1)) + 
                ((p - 1) * rland$intparam$stages), NperPop[p], 
                replace = T, prob = ev)
            im[, 2] <- rep(0, NperPop[p])
            im[, 3] <- im[, 2]
            im[, 5:landscape.democol()] <- 0
            indlist[[p]] <- im
        }
        individuals <- do.call("rbind", indlist)
        individuals <- individuals[order(individuals[, 1]), ]
        individuals[, 4] <- 1:dim(individuals)[1]
        rland$individuals <- matrix(as.integer(individuals), 
            nrow = dim(individuals)[1])
        rland$intparam$nextid <- dim(individuals)[1] + 
            1
        rland
    }
}
