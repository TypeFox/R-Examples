`sim.autoMarkers` <-
function(ploidy.level, dose.proportion, n.markers=500, n.individuals=200,
         seg.ratios, no.dosage.classes,
         type.parents=c("heterogeneous","homozygous"),
         marker.names=paste("M",1:n.markers,sep="."),
         individual.names=paste("X",1:n.individuals,sep="."),
         overdispersion=FALSE, shape1=50, seed)
{
  
  ## Description: Simulates dominant markers from an autopolyploid
  ## cross given the ploidy level or expected segregation ratios and
  ## the proportions in each dosage marker class. This may be chosen
  ## from tetraploid to heccaidecaploid or instead the segregation
  ## ratios may be specified explicitly. For later use in simulation
  ## studies, other parameters such as the true dosage of each
  ## marker are also returned.

  ## Arguments:
  ## ploidy.level: the number of homologous chromosomes, either as
  ##               numeric (single value) or as a character string
  ##               containing type
  ##               tetraploid, hexaploid, octoploid, ....
  ## dose.proportion:  the proportion of markers to be simulated in each
  ##                   dosage class. Note that the exact number will be
  ##                   randomly generated from the multinomial distribution
  ## n.markers:   number of markers (Default: 500)
  ## n.individuals: number of individuals in the cross (Default: 200)
  ## seg.ratios: numeric vector containing segregation proportion to be
  ##             supplied if you wish to overide automatic calculations using
  ##             ploidy.level
  ## no.dosage.classes: only generate markers for the first
  ##                    'no.dosage.classes' classes (if set)
  ## type.parents: heterogeneous for (1,0) or (0,1) homozygous for (1,1)
  ##               (default: heterogeneous)
  ## marker.names: labels for markers (Default: M.1 ... M.n.markers)
  ## individual.names: labels for offspring (Default: ... X.j ... )
  ## overdispersion: logical indicating overdispersion (Default: =FALSE)
  ## shape1: shape1 parameter(s) for the beta distribution used to generate
  ##         the Binomial parameter p, either of length 1 or no.dosage.classes
  ##         Default: a small amount of overdispersion
  ## shape2: CALCULATED - NOT SPECIFIED because once expected segratio
  ##         generated, this is fixed
  ## seed: integer used to set seed for random number generator (RNG) 
  ##               which (if set) may be used to reproduce results

  ## Values:
  ## markers: matrix of 0,1 dominant markers with individuals as cols and
  ##          rows as markers
  ## E.segRatio: expected segregation porportions, list with components
  ##             ratio:    segregation proportions
  ##             ploidy.level:   level of ploidy 4,6,8, ....
  ##             ploidy.name: tetraploid, ..., unknown
  ## type.parents: heterogeneous for (1,0) or (0,1) homozygous for (1,1)
  ## dose.proportion: proportions of markers set for each dosage class
  ## n.markers:   number of markers (Default: 500)
  ## n.individuals: number of individuals in the cross (Default: 200)
  ## true.doses: list containing
  ##             dosage: doses generated for each marker for simulation
  ##             table.dosages: summary of no.s in each dosage
  ##             names: names SDRF, DDRF, etc
  ## seg.ratios: object of class segRatio containing segregation ratios
  ## time.generated: time/date when data set generated
  ## seed:     seed for random number generator seed which could be used to
  ##           reproduce results (I hope)
  ## overdispersion: either a list with overdispersion=logical depending
  ##                 on whether overdispersion is set and if set then two
  ##                 extra numeric variables 'shape1' and 'shape2' for the
  ##                 beta distrubution for generating probs
  ## call:          matches arguments when function called

  ## params for writing/testing:
  ##   n.markers=20; n.individuals=10; ploidy.level=4;  dose.proportion=c(.8,.2)
  ##   type.parents="heterogeneous"; na.proportion=0; misclass=0
  ##   marker.names=paste("M",1:n.markers,sep=".")
  ##   individual.names=paste("X",1:n.individuals,sep=".")
  ##   type <- type.parents
  ##   na.proportion <- list(indiv=c(0.1,0.1),marker=c(0.1,0.2))

  if (missing(seed)) {
    ## save.seed <- .Random.seed # if same results required then set this
  } else {
    set.seed(seed)
  }
  
  type <- match.arg(type.parents)
  E.segRatio <- expected.segRatio(ploidy.level, type.parents=type)

  if (!missing(seg.ratios)){
    E.segRatio$ratio <- seg.ratios
    cat("Warning: default segregation proportions set to")
    print(E.segRatio$ratio); cat("\n")
  }
  ##print(E.segRatio)

  if (length(dose.proportion) > length(E.segRatio$ratio))
    stop("Too many dosage class proportions specified")
  if (length(dose.proportion) < length(E.segRatio$ratio)) {
    cat("Warning: segregation proportions truncated to length",
        length(dose.proportion),"\n")
    E.segRatio$ratio  <- E.segRatio$ratio[1:length(dose.proportion)]
  }

  ## overspersion set

  if (overdispersion & length(shape1)==1) {

  ##  if (overdispersion & shape1==0)
  ##    stop("Overdispersion: Beta parameter 'shape1' must be specified")
  
    shape1 <- rep(shape1,length(E.segRatio$ratio))
  }
  
  ## if set, only use first no.dosage.classes
  if (!missing(no.dosage.classes))
    {
      if (no.dosage.classes>length(E.segRatio$ratio)) {
        cat("Warning no.dosage.classes too large so ignored\n")
      } else {
        if (no.dosage.classes < length(E.segRatio$ratio)) {
          cat("Warning: segregation proportions truncated to length",
              no.dosage.classes,"and segregation proportions adjusted\n")
          E.segRatio$ratio  <- E.segRatio$ratio[1:no.dosage.classes]
          E.segRatio$ratio  <- E.segRatio$ratio/sum(E.segRatio$ratio)
          dose.proportion <- dose.proportion[1:no.dosage.classes]

          if (overdispersion & length(shape1)==1) {
            shape1 <- shape1[1:no.dosage.classes]

          }
        }
      }
    }

  ## print(shape1)
  
  ## dose.proportions - make sure add to 1 and give names if necessary
  
  if (sum(dose.proportion) != 1) {
    cat("Warning: proportions rescaled to sum to 1\n")
    dose.proportion <- dose.proportion/(sum(dose.proportion))
  }

  if (length(names(dose.proportion))==0) {
    names(dose.proportion) <- names(E.segRatio$ratio)
  }

  param.overdispersion <- list(overdispersion=overdispersion)
  if (overdispersion) {
    names(shape1) <- names(E.segRatio$ratio)
    param.overdispersion$shape1 <- shape1
    shape2 <- shape1/E.segRatio$ratio - shape1
    param.overdispersion$shape2 <- shape2
  }

  ## print(shape2)
  
  ## generate number of each type of marker (
  ## NB: maybe can do all this in one hit by multiplying class proportions
  ##     by seg proportions and even for different experiments - hmmm

  no.in.doses <- rmultinom(1, n.markers, dose.proportion)
  dimnames(no.in.doses) <-  list(names(E.segRatio$ratio), "No.markers")

  ## generate values for each marker given dosage

  markers <- matrix(NA, ncol=n.individuals, nrow=n.markers, 
                    dimnames=list(marker.names,individual.names))
  true.dosage <- rep(1:length(E.segRatio$ratio),times=no.in.doses)
  true.dosage.names <- names(E.segRatio$ratio)[true.dosage]
  if (overdispersion) {
    for (i in 1:n.markers) {
      markers[i,] <- rbinom(n.individuals, 1,
                            prob=rbeta( n=1, shape1=shape1[true.dosage[i]],
                                        shape2=shape2[true.dosage[i]] ))
    }
  } else {
    for (i in 1:n.markers) {
      markers[i,] <- rbinom(n.individuals, 1,
                            E.segRatio$ratio[true.dosage[i]])
    }
  }

  time.generated <- date()

  seg.ratios <- segregationRatios(markers)
  
  res <- list (markers=markers, E.segRatio=E.segRatio,
               type.parents=type, ploidy.level=ploidy.level,
               n.markers=n.markers,n.individuals=n.individuals,
               dose.proportion=dose.proportion,
               true.doses=list(dosage=true.dosage, table.doses= no.in.doses, 
                 names=true.dosage.names), seg.ratios=seg.ratios,
               time.generated=time.generated,
               overdispersion=param.overdispersion, call=match.call())

  if (!missing(seed)) res$seed <-  seed

  oldClass(res) <- "simAutoMarkers"
  return(res)
}

