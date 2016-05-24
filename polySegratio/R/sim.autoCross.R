`sim.autoCross` <-
function(ploidy.level,
           prop.par.type=structure(c(0.4,0.4,0.2),names=c("p10","p01","p11")),
           n.markers=500, n.individuals=200,
           dose.proportion, 
           true.seg.ratios,
           no.dosage.classes,
           marker.names=paste("M",1:n.markers,sep="."),
           individual.names=paste("X",1:n.individuals,sep="."),
           parent.names=c("P.1","P.2"), seed)
{
  ## Description: Simulates dominant markers from an autopolyploid
  ## cross given the ploidy level and/or expected segregation ratios
  ## and the proportions in each dosage marker class. This is a
  ## wrapper to sim.automarkers to generate markers for '10', '01' and
  ## '11' parents and so arguments for each parental type may be set
  ## as a list with components 'p10','p01','p11'. All parameters
  ## except the proportions of marker dosage types can be left at the
  ## default. If only one value is set, then individual list
  ## components will be assumed to be equal. The marker matrix is
  ## prepended with parental marker alleles. An alternative is to
  ## simply create each group using sim.automarkers and cbind them.

  ## Arguments:
  ## ploidy.level: the number of homologous chromosomes, either as
  ##               numeric (single value) or as a character string
  ##               containing type
  ##               tetraploid, hexaploid, octoploid, ....
  ## prop.par.type: the proportion of markers generated from
  ##                each parental type '10', '01' and '11'. Note that
  ##                the exact number will be
  ##                randomly generated from the multinomial distribution
  ##                (Default: c(0.4,0.4,0.2))
  ## n.markers:   number of markers (Default: 500)
  ## n.individuals: number of individuals in the cross (Default: 200)
  ## dose.proportion:  the proportion of markers to be simulated in each
  ##                   dosage class. Note that the exact number will be
  ##                   randomly generated from the multinomial distribution
  ##                   NB: If a vector is supplied the dose.proportion is same
  ##                   for each parental type otherwise as list with
  ##                   components 'p01', 'p10' and 'p11'
  ## no.dosage.classes: vector containing the number of dosage classes
  ## true.seg.ratios: numeric vector containing segregation proportion to be
  ##             supplied if you wish to overide automatic calculations using
  ##             ploidy.level
  ## parent.names: default for first 2 rows of marker matrix containing
  ##               parental markers
  ## seed:     random number generator (RNG) state for random number
  ##               which will be set at start to reproduce results

  ## Values:
  ## markers: matrix of 0,1 dominant markers with indiviuals as cols and
  ##          rows as markers with parents prepended 
  ## true.doses: list containing
  ##             doseage: doses generated for each marker for simulation
  ##             table.dosages: summary of no.s in each dosage
  ##             names: names SDRF, DDRF, etc
  ## p10, p01, p11: objects of class simAutoMarkers for each parental type
  ## ploidy.level: the number of homologous chromosomes, either as
  ##               numeric (single value) or as a character string
  ##               containing type
  ##               tetraploid, hexaploid, octoploid, ...
  ## prop.par.type: proportion of markers for each parental type p10, p01, p11
  ## n.markers: no. of markers
  ## n.individuals: no.individuals
  ## dose.proportion: proportion in exach dose- if numeric is the same for
  ##                  p10, p01, p11, else a list with components p10, p01, p11
  ## seg.ratios: object of class segRatio containing segregation ratios
  ## no.dosage.classes: no. in each dosage class
  ## no.parType: no. in each parental type
  ## time.generated: time/date when data set generated
  ## seed:     seed for random number generator seed which could be used to
  ##           reproduce results (I hope)
  ## call:          matches arguments when function called


  if (!missing(seed)) {
    set.seed(seed)
  }

  ## set up no.s for each parental cross
  par.names <- c("p10","p01","p11")
  type <- c("heter","heter","homo")

  ## proportions of each parental cross type
  if (length(prop.par.type) != 3)
    stop("'prop.par.type' is proportions of '10', '01' and '11' markers")
  if (length(names(prop.par.type))==0)
    names(prop.par.type) <- par.names
  if (sum(prop.par.type) != 1) {
    cat("Warning: proportions 'prop.par.type' adjusted to sum to 1")
    prop.par.type <- prop.par.type/sum(prop.par.type)
  }

  no.parType <- rmultinom(1, n.markers, prop.par.type)
  colnames(no.parType) <- "No.markers"

  ## for each parental type determine no. of markers for each dosage

  ## first determine dose.proportion for each parental type
  if (mode(dose.proportion) == "numeric") {
    dose.proportion <- list(p10=dose.proportion, p01=dose.proportion,
                            p11=dose.proportion)
    warning("Dose proportions set the same for all parental types - dubious")
  } else {
    if (mode(dose.proportion) != "list")
      stop("'dose.proportion' must be numeric or list")
    if (length(pmatch(names(dose.proportion),par.names))!=3) { ## |
      ##        any(is.na(pmatch(names(dose.proportion),par.names))) )
      stop("'dose.proportion' must contain proportions for each dose p10, p01 and p11")
    }
  }
    
  par.types <- pmatch(par.names,names(dose.proportion))
  
  ## first determine dose.proportion for each parental type
  if (!missing(true.seg.ratios)) {
    if (mode(true.seg.ratios) != "list")
      stop("'true.seg.ratios' must be list with at least one of p10, p01 or p11")
    ## check names - no strange ones
    if (length(pmatch(names(true.seg.ratios),par.names[par.types]))<1 |
        any(is.na(pmatch(names(true.seg.ratios),par.names[par.types]))) )
      stop("'dose.proportion' must contain proportions for each dose for at least one of p10, p01 or p11")
  }

  ## set some parameters if not set already
  if (missing(no.dosage.classes)){
    no.dosage.classes <- dose.proportion
    for (i in 1:length(no.dosage.classes)) {
      no.dosage.classes[[i]] <- length(dose.proportion[[i]])
    }
  } else {
    if (!mode(no.dosage.classes) == "list"|mode(no.dosage.classes) == "numeric")
      stop("'no.dosage.classes' must be list or numeric with components p10, p01 and p11")
    ## check names - no strange ones
    if (length(no.dosage.classes) != 3)
      stop("Three dosage classes are required")
    if (length(pmatch(names(no.dosage.classes),par.names[par.types]))!=3 |
        any(is.na(pmatch(names(no.dosage.classes),par.names[par.types]))) )
      stop("'no.dosage.classes' must contain no. of proportions for each dose for  p10, p01 and p11")
    ## now shorten dose.proportion if necessary
    if (1 %in% par.types){
      if (no.dosage.classes$p10 < length(dose.proportion$p10))
        dose.proportion$p10 <- dose.proportion$p01[1:no.dosage.classes$p10]
    }
    if (2 %in% par.types){
      if (no.dosage.classes$p01 < length(dose.proportion$p01))
        dose.proportion$p01 <- dose.proportion$p01[1:no.dosage.classes$p01]
    }
    if (3 %in% par.types){
      if (no.dosage.classes$p11 < length(dose.proportion$p11))
        dose.proportion$p11 <- dose.proportion$p01[1:no.dosage.classes$p11]
    }
  }
  
  ## now generate the markers separately for each parental type, put
  ## them together and store various info suitable for simulation
  ## studies
  ## NB: real data may be coreced into similar form using split.markers
  ##     although no true values will be stored

  sim.data <- list(markers=NA, true.dosage=NA, name.true.dose=NA)
  
  m.first <- 0
  j <- length(sim.data)
  markers <- NULL
  dosage <- NULL
  name.dose <- NULL
  for (i in par.types) {

    j <- j+1
    m.no <- m.first + c(1:no.parType[i])
    m.first <- m.first+no.parType[i]

 ##   if (missing(true.seg.ratios)) {
    sim.data[[j]] <-
      sim.autoMarkers(ploidy.level, dose.proportion=dose.proportion[[i]],
                      n.markers=no.parType[i], n.individuals=n.individuals,
                      type.parents=type[i], marker.names=marker.names[m.no],
                      no.dosage.classes=no.dosage.classes[[i]],
                      individual.names=individual.names)
    markers <- rbind(markers,sim.data[[j]]$markers)
    dosage <- c(dosage,sim.data[[j]]$true.dose$dosage)
    name.dose <- c(name.dose,sim.data[[j]]$true.dose$names)
    names(sim.data)[j] <- par.names[i]
  }

  sim.data$true.dosage <- dosage
  sim.data$name.true.dose <- name.dose

  # add.in parental markers

  tmp <- matrix(c(rep(c(1,0),no.parType[1]) , rep(c(0,1),no.parType[2]),
         rep(c(1,1),no.parType[3])), ncol=2, nrow=n.markers, byrow=TRUE)
  colnames(tmp) <- parent.names
  sim.data$markers <- cbind(tmp,markers)
  sim.data$seg.ratios <- segregationRatios(markers)

  if (! missing(seed)) sim.data$seed <- seed

  res <- c( sim.data,
                list(ploidy.level=ploidy.level,
                     prop.par.type=prop.par.type, n.markers=n.markers,
                     n.individuals=n.individuals,
                     dose.proportion=dose.proportion,
                     no.dosage.classes=no.dosage.classes,
                     no.parType=no.parType, time.generated=date(),
                     call=match.call() ))

  
  oldClass(res) <- "simAutoCross"
  return(res)

}

