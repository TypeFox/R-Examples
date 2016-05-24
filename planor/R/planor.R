## "planor" package 
## A package dedicated to automatic generation of regular fractional factorial designs,
## including fractional designs, orthogonal block designs, row-column designs and split-plots.
##  Research Unit MaIAGE, INRA, Jouy en Josas, France.
## EXAMPLES
## DESIGN SPECIFICATIONS
## Treatments: four 3-level factors A, B, C, D
## Units: 27 in 3 blocks of size 9
## Non-negligible factorial terms:
##   block + A + B + C + D + A:B + A:C + A:D + B:C + B:D + C:D
## Factorial terms to estimate:
##   A + B + C + D
## 1. DIRECT GENERATION, USING 'regular.design'
## >  mydesign <- regular.design(factors=c("block", LETTERS[1:4]),
## >    nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
## >    nunits=3^3, randomize=~block/UNITS)

## DUMMY ANALYSIS
## Here we omit two-factor interactions from the model, so they are 
## confounded with the residuals (but not with ABCD main effects)
## >  set.seed(123)
## >  mydesigndata=mydesign@design
## >  mydesigndata$Y <- runif(27)
## >  mydesign.aov <- aov(Y ~ block + A + B + C + D, data=mydesigndata)
## >  summary(mydesign.aov)
## 2. STEP-BY-STEP GENERATION, USING 'planor.designkey'
## >  F0 <- planor.factors(factors=c( "block", LETTERS[1:4]), nlevels=rep(3,5),
## >    block=~block)
## >  M0 <- planor.model(model=~block+(A+B+C+D)^2, estimate=~A+B+C+D) 
## >  K0 <- planor.designkey(factors=F0, model=M0, nunits=3^3, max.sol=2)
## >  summary(K0)
## >  mydesign.S4 <- planor.design(key=K0, select=2)







  FREE <- TRUE

##---------------------------------------------------------------------------
##  MAIN STRUCTURES (specific classes or not)
## 1. Names for specific objects or arguments (see also INTERNAL NOTATION CONVENTIONS)
## factors: a data.frame (usually with 0 row); can be created by planor.factors
## modlist: a list of (model-formula, estimate-formula) pairs,
##              usually declared in planor.model
## modmat: a list of (model-terms-matrix, estimate-terms-matrix) pairs,
##              usually declared in planor.model
## 2. Classes
## designfactors : an S4 class, typically an output from planor.factors
##         fact.info : a dataframe with one row per factor and columns
##                     'names', 'nlev', 'ordered'
##         pseudo.info : a dataframe with one row per pseudofactor and columns
##                  'parent', 'nlev', 'block', 'ordered', 'model', 'basic'
##         levels : a list giving the levels of the factors
## designkey : an S4 class, typically an output from pick
##         main: a single design-key solution = a list with one key matrix for each prime. Each one is a 'keymatrix' object.
##         factors: the 'designfactors' object that defines the factors
##         model: the list of components of type c(model,estimate)
##                containing the model and estimate specifications
##         nunits: the number of units in the design.
##         recursive: logical, TRUE if the design has been constructed recursively
## keymatrix : an S4 class, a component of designkey or keyring
##         main: a solution  for a prime
##         p: value of the prime
## keyring : an S4 class, a component of listofkeyrings
##         main: a set of design-key solutions = a list of solutions for a prime. Each one is a keymatrix.
##         p: value of the prime
##         LIB: list of the rownames and colnames of the key matrices
## listofdesignkeys : an S4 class, typically an output from planor.designkey when the research is recursive
##         main: a list of design-key solutions; each component
##                  of main is a whole solution list across the different primes. It is an object of class 'designkey'
##         factors: the 'designfactors' object that defines the factors
##         model: the list of components of type c(model,estimate)
##                containing the model and estimate specifications
##         nunits: the number of units in the design.
## listofkeyrings : an S4 class, typically an output from planor.designkey when the research is not recursive
##         main: a list of design-key solutions; each component
##                  of main is  a list of design
##                  key solutions for a given prime. It is an object of class 'keyring'
##         factors: the 'designfactors' object that defines the factors
##         model: the list of components of type c(model,estimate)
##                containing the model and estimate specifications
##         nunits: the number of units in the design.
## planordesign : an S4 class, typically an output from planor.design.designkey
##         design: a dataframe containing the final design
##         factors: the 'designfactors' object that defines the factors
##         model: the modlist containing the model and estimate specifications
##         designkey: a list which contains the designkey matrices used to create the object
##         nunits: the number of units in the design.
##         recursive: a "logical" equal to TRUE if the design has been constructed recursively
##---------------------------------------------------------------------------
##  FUNCTIONS IN THIS FILE
##---------------------------------------------------------------------------
## 1. DESIGN INPUT AND SPECIFICATIONS
##---------------------------------------------------------------------------
## planor.factors <- function(factors=NULL, nlevels=NULL,
##                            block=NULL,
##                            ordered=NULL,
##                            hierarchy=NULL,
##                            dummy=FALSE)
## planor.model <- function(model, estimate, listofmodels,
##                          resolution, factors)
## planor.designkey <- function(
##                              ## arguments for planor.factors
##                              factors=NULL, nlevels=NULL,
##                              block=NULL,
##                              ordered=NULL,
##                              hierarchy=NULL,
##                              ## arguments for planor.model
##                              model=NULL, estimate=NULL,
##                              listofmodels=NULL,
##                              ## or resolution for automatic generation of the model
##                              resolution=NULL,
##                              ## other arguments
##                              nunits=NULL, base=NULL,
##                              max.sol=1, randomsearch=FALSE,
## 			     verbose=TRUE)
##---------------------------------------------------------------------------
## 2. MAJOR CALCULUS FUNCTIONS
##---------------------------------------------------------------------------
## regular.design <- function(
##                    ## arguments for planor.factors
##                    factors=NULL, nlevels=NULL,
##                     block=NULL, ordered=NULL, hierarchy=NULL,
##                     ## arguments for planor.model
##                     model=NULL, estimate=NULL, listofmodels=NULL,
##                     ## or resolution for automatic generation of the model
##                     resolution=NULL,
##                     ## other arguments
##                     nunits=NULL, base=NULL,
##                     ## design stage
##                     randomize=NULL,
##                     randomsearch=FALSE,
##                     ## information
##                     output="planordesign",
##                     verbose=FALSE,
##                     ...)
## planor.hierarchy  <- function(factors, hierarchy=NULL)
## planor.pseudofactors <- function(factors)
## planor.harmonize <- function(
##                              ## arguments for planor.factors
##                              factors=NULL, nlevels=NULL,
##                              ordered=NULL,
##                              hierarchy=NULL,
##                              ## arguments for planor.model
##                              model=NULL, estimate=NULL,
##                              listofmodels=NULL,
##                              ## other arguments
##                              base=NULL)
## planor.modelterms <- function(modlist)
## planor.ineligibleterms <- function(modmat)
## planor.ineligibleset <- function(factors, ineligible)
## planor.designkey.basep <- function(p, r, ineligible, hierarchy, predefined=NULL,
##                                   max.sol=1, randomsearch=FALSE,
##				    verbose=FALSE){
## planor.designkey.recursive <- function(k,nb.sol,PVuqf,NPuqf,
##                                        b.INELIGtpf,NIVtpf,
##                                        b.INELIGuqf,PREDEF,
##                                        max.sol=1)
##---------------------------------------------------------------------------
## 3. SECONDARY CALCULUS FUNCTIONS
##---------------------------------------------------------------------------
## planor.kernelcheck.basep <- function(PhiStar, admissible, IneligibleSet, p)
## weightorder.basep <- function(mat,p,factnum,blocklog)

## printpower <- function(p.k)
## wprofile <- function(x)

##---------------------------------------------------------------------------
## 4. OUTPUT FUNCTIONS
##---------------------------------------------------------------------------
## planor.design.levels <- function(key,start=1)
##---------------------------------------------------------------------------
## INTERNAL NOTATION CONVENTIONS (not used in all functions any more)
## Many scalar and vector objects are named by combining an uppercase prefix
## and a lowercase suffix. This is inspired by, but different from, the
## conventions used in Planor or in Kobilinsky 2000. For lists or vectors, the
## lowercase suffix is directly related to the object length or row number
## (or occasionally column number).
## (Occasionnally, the lower case part is put below the upper case part for a
## matrix. In that case the lower case part refers to columns.)
## Uppercase prefix:
##  N : number of
##  LIB : names of
##  LAB : level names of
##  NIV : level numbers of
##  BLOC : block-or-not-block indicator of
##  NQ : number of distinct primes in the level-number decomposition of (Q for quasifactor)
##  NP : number of primes in the level-number decomposition of (P for pseudofactor)
##  PV : prime values
##  TFACT : treatment factor at the origin of
## Uppercase prefix for matrices:
##  INELIG : ineligible terms (indicator of presence of ... in the terms)
## Lowercase suffix:
##  tf : treatment factors given by the user (ft in Kobilinsky 2000)
##  tqf : treatment quasifactors (one quasifactor per distinct prime per factor)
##  tpf : treatment pseudofactors (one pseudofactor per prime per factor)
##  bf : basic factor
##  btf : basic-and-treatment factor (cf. some basic factors may be absent from
##        the models and estimate formulae, then there are fewer btf than bf))
##  uf : unit or base factor
##  uqf : unit quasifactors (one quasifactor per distinct prime in the unit number)
##  upf : unit pseudofactor
##  prime : prime (remplacer par Sylow?)
##  ipft : ineligible pseudofactorial term
##---------------------------------------------------------------------------


##---------------------------------------------------------------------------
## MAIN FUNCTIONS
##---------------------------------------------------------------------------
## 1. DESIGN INPUT AND SPECIFICATIONS
##---------------------------------------------------------------------------


## "planor.factors"
##   A user-friendly function to create an
##   object of class 'designfactors',
##    either by giving the factor names and level numbers, or by giving
##    a list of factor levels. Both ways can be used in the same call
##
## ARGUMENTS
##   - factors:  a character vector of factor names, or possibly a scalar, a dataframe
##                     or a list (see DETAILS)
##   - nlevels:  a vector of level numbers for each factor name (see DETAILS)
##   - block:     an additive model formula to indicate the block factors
##   - ordered:   an additive model formula to indicate the ordered factors
##   - hierarchy:   a formula or a list of formulae to indicate hierarchy relationships between factors
## NOTE
##     The basic usage is to specify the names of the factors by a character vector of length 'n' in argument 'factors'
## and their numbers of levels by a numeric vector of length 'n' in argument 'nlevels'.
## Alternatively, the 'factors' argument can be an integer 'n',
## in which case the first 'n' capital letters of the alphabet are used as factor names.
## If 'nlevels' is a scalar, it is considered that all factors have this number of levels.
## There are two more possibilities.
## If 'factors' is a dataframe, the factors in this dataframe are extracted together with their levels.
## Finally 'factors' can be a list of vectors (but not a dataframe),
## where each vector in the list has the name of a factor and contains the levels of that factor.
## Note that 'nlevels' is ignored in these latter two cases. See the examples.
##     The argument 'block' is used by the functions that give the properties of the design keys.
##     The argument 'ordered' must be either
##     an integer or a vector of length equal to
##     the number of factors.
## RETURN
##   An object of class 'designfactors'
## EXAMPLES
##    planor.factors(c("A","B","C","P"),c(2,3,6,3))
##     planor.factors(LETTERS[1:12],2)
##     planor.factors(12,2)
##     planor.factors( c("A","B","Block"), 3, block=~Block )
##     zz <- planor.factors( c("A","B","Block"), c(2,3,5))
##     zz@@levels$A <- c("plus","moins")
## planor.factors(factors=list(A=c("plus","moins"), B=1:3, Block=1:5))
## AB <- data.frame( A=c(rep(c("a","b"),3)), B=rep(c("z","zz","zzz"),rep(2,3)), C=1:6  )
## planor.factors(factors=AB)
## -----------------------------------------------------

planor.factors <- function(factors=NULL, nlevels=NULL,
                           block=NULL,
                           ordered=NULL,
                           hierarchy=NULL,
                           dummy=FALSE){
  ## SPECIAL CASES FOR THE ARGUMENTS factors AND nlevels
  ## First case: 'factors' is a scalar or a vector of names
  if( is.numeric(factors) | is.character(factors) ){
    if( is.numeric(factors) && (length(factors)==1) ){
      FACT.names <- LETTERS[seq_len(factors)] }
    else{
      FACT.names <- factors }
    if(length(nlevels)==1){ nlevels <- rep(nlevels, length(FACT.names)) }
    FACT.levels <- lapply( as.list(nlevels), seq_len )
    names(FACT.levels) <- FACT.names
  }
  ## Second case: 'factors' is a dataframe
  else if( is.data.frame(factors) ){
    select <- unlist(lapply(factors, is.factor))
    FACT.levels <- lapply( factors[select], levels )
    FACT.names <- names(FACT.levels)
  }
  ## Third case: 'factors' is a list of factor levels
  else if( is.list(factors) ){
    FACT.levels <- factors
    FACT.names <- names(FACT.levels)
  }
  ## Fourth case: trouble
  else{ stop("Problem with argument 'factors', 'nlevels' or both") }
  
  ## information on the factors:
  FACT.N <- length(FACT.names)
  FACT.nlev <- unlist(lapply(FACT.levels, length))
  ## - block factors


  if(is.null(block)){
    block <- rep(FALSE, FACT.N) }
  else{
    block <- FACT.names %in% all.vars(block) }
  ## - ordered factors
  if(is.null(ordered)){
    ordered <- rep(FALSE, FACT.N) }
  else{
    ordered <- FACT.names %in% all.vars(ordered) }
  ## storage of the 'factors' information
  fact.info <- data.frame(nlev = FACT.nlev,
                          block = block,
                          ordered = ordered,
                          model = NA,
                          basic = NA,
                          dummy = dummy,
                          stringsAsFactors=FALSE)
  rownames(fact.info) <- FACT.names
  FACTORS <- new("designfactors",
                 fact.info = fact.info,
                 levels = FACT.levels)
  ## hierarchy
  if(!is.null(hierarchy)){
    FACTORS <- planor.hierarchy(FACTORS, hierarchy)
  }
  ## information on the pseudofactors' decomposition
  FACTORS <- planor.pseudofactors(FACTORS)
  ##
  return(FACTORS)
} ## end planor.factors
##---------------------------------------------------------------------------
## "planor.model"
##    A function to declare the factorial terms that must be considered
## as non-negligible and the factorial terms that must be estimable
## when the experiment will be analysed.
## ARGUMENTS
##   - model:   formula of the main model
##   - estimate:  optional formula specifying the factorial terms to estimate;
##               if missing, it is considered that all model terms have to be
##               estimated
##   - listofmodels: list of c(model, estimate) pairs, where model and estimate are
##              formulae; using several pairs allows more flexibility in the design
##              constraints  (see Kobilinsky, 2005); estimate is optional
##   - resolution: integer, larger than or equal to 3, equal to the design resolution. When set, the  'model' argument is ignored.
## See  Note.
##   - factors: a 'designfactors' object,
##         typically an output from 'planor.factors'.
##         Should be set, when the 'resolution' argument is set.
## NOTE
##  The user must specify one or the other set of arguments: 
## 1/ either, 'model' or 'listofmodels' or both 
## 2/ or, 'resolution' and 'factors', and possibly 'listofmodels'. 
## When 'model' and 'resolution' are all set,
## 'model' is ignored. 
## The second case, ---  'resolution' and 'factors' are set ---,
## causes the automatic generation of the (model, estimate) pairs: 
## Assuming 'S' denotes the additive formula including all factors, 
## - if 'resolution' is odd, the model formula is '~(S)^(resolution-1)/2',
## - if 'resolution' is even, the model  formula is '~(S)^(resolution/2)' and the estimate formula is  '~(S)^(resolution/2)-1'
## RETURN
##     A list of c(model, estimate) pairs, where model and estimate are
##   formulae
## EXAMPLES
##  planor.model(~A+B+C+D+A:B,~A+B+C+D, listofmodels=list(c(~E+F,~E)))
##  planor.model(~A+B+C+D+A:B,~A+B+C+D, listofmodels=list(c(~E+F,~E), ~G, ~H, c(~M+N,~N)))
## planor.model(resolution=4, factors=planor.factors( factors=c(LETTERS[1:4]),  nlevels=rep(2,4)))
## ------------------------------------------------------

planor.model <- function(model, estimate, listofmodels,
                         resolution, factors){
  ##
  if (!missing(resolution) && !is.null(resolution)) {
    ## Automatic generation of the model, estimate
    if (missing(factors)) {
        stop("Argument factors missing")
      }

     res <- generate.model(resolution, factors)
     model <- res$model
     estimate <- res$estimate
  } ## end (!missing(resolution))

  if(!missing(model) && !is.null(model)){
    if(missing(estimate) || is.null(estimate)) estimate <- model
    if(missing(listofmodels)|| is.null(listofmodels)) listofmodels <- list( c(model,estimate) )
    else listofmodels <- c( list(c(model,estimate)), listofmodels )
  }
  for(i in seq_along(listofmodels)){
    if(class(listofmodels[[i]]) == "formula"){
      listofmodels[[i]] <- c(model=listofmodels[[i]],estimate=listofmodels[[i]])
    }
    if(is.null(names(listofmodels[[i]]))){
      names(listofmodels[[i]]) <- c("Model","Estimate")
    }
  }
  return(listofmodels)
} ## end planor.model

##---------------------------------------------------------------------------
## "generate.model" 
## Internal function
##    Automatic generation of the model from the design resolution and the factors
##
## TITLE Model generation from the design resolution and the factors
## ARGUMENTS
##   - resolution: integer, larger than or equal to 3, equal to the design resolution.
##   - factors: a 'designfactors' object,
##         typically an output from 'planor.factors'
##   - listofmodels: optional list of c(model, estimate) pairs, where model and estimate are
##              formulae; see  'planor.model'.
## RETURN
##     A list of c(model, estimate) pairs, where model and estimate are
##   formulae. 
## NOTE
## Assuming 'S' denotes the additive formula including all factors, 
## - if 'resolution' is odd, the model formula is '~(S)^(resolution-1)/2',
## - if 'resolution' is even, the model  formula is '~(S)^(resolution/2)' and the estimate formula is  '~(S)^(resolution/2)-1'
## EXAMPLES
##  F <- planor.factors(c(LETTERS[1:4], "bloc"),nlevels=rep(2,5))
##  M <- generate.model(3,F)


generate.model <- function(resolution, factors) {
  resolution <- as.integer(resolution)
  if (resolution <3) {
    stop("Resolution must be larger than or equal to 3")
  }
  
  estimate <- NULL
  model <- NULL
  ## Build the model expression in a character string
  factmod <- paste(rownames(factors@fact.info),collapse="+")
  factmod <- paste("~(", factmod,")", sep="")
  mod <- paste("model <- ", factmod)
  if ( (resolution%%2) ==0) {
    ## resolution is even
    degremod <-  resolution/2
    degreestimate <- (resolution/2)-1
    es <- paste("estimate <-", factmod)
    if (degreestimate <= 0) {
      stop("Negative power in the estimate formula")
    }
    if (degreestimate > 1) {
      es <- paste(es,"^",degreestimate,  sep="")
    }
    eval(parse(text=es))
  } ## end  resolution is even
  else {
    ## resolution is odd
    degremod <- (resolution-1)/2
  }
  if (degremod  <= 0) {
      stop("Negative power in the model  formula")
    }
  if (degremod >1) {
    mod <- paste(mod,"^",degremod, sep="")
  }
  eval(parse(text=mod))
  return( list(model=model, estimate=estimate))
} ## end generate.model

##---------------------------------------------------------------------------
## "planor.designkey" 
##  Search for a design key in a possibly mixed factorial context
##  ARGUMENTS
##   - factors:  'factors' can be of 2 types:
## either,  an object of class 'designfactors' (typically an output
##               from 'planor.factors'),
## or a character vector of factor names.
## In this last case, the arguments 'nlevels', 'ordered', 'hierarchy' can also be set.
##   - nlevels:  a vector of level numbers for each factor name.
## Ignored if 'factors' is an object of class 'designfactors'.
##   - block :  an additive model formula to indicate the block factors
##   - ordered:   an additive model formula to indicate the ordered factors
## Ignored if 'factors' is an object of class 'designfactors'.
##   - hierarchy:   a formula or a list of formulae to indicate hierarchy relationships between factors
## Ignored if 'factors' is an object of class 'designfactors'.
## - model: either a list of model-estimate pairs of formulae, typically an output
##               from 'planor.model', or
##  the   formula of the main model. In this last case, the arguments
## 'estimate' and 'listofmodels' can also be set.
##   - estimate:  if 'model' is a formula,
## optional formula specifying the factorial terms to estimate;
##               if missing, it is considered that all model terms have to be
##               estimated.
## Ignored if 'model' is a list.
##   - listofmodels:  if 'model' is a formula,
## list of c(model, estimate) pairs, where model and estimate are
##              formulae; using several pairs allows more flexibility in the design
##              constraints  (see Kobilinsky, 2005); estimate is optional.
## Ignored if 'model' is a list.
##  - resolution: optional integer equal to the design resolution.
## When set, the model is generated from the factors and
## 'model' and 'estimate' are ignored. See Note.
## - nunits: a scalar ; the total number of units in the design
## - base: an optional additive formula to specify the basic factors
## - max.sol: maximum number of solutions before exit
## - randomsearch: a 'logical';
##       if TRUE, the order of admissible elements is randomised
##                   at each new visit forward to 'j'
## - verbose: a 'logical';
## 		  if TRUE, verbose display
## NOTE
## The 'base' formula must be given as an additive model on the basic factors,
##     and the basic factors must be specified in the 'factors' argument. 
## When  'resolution' is set,
## the (model, estimate) pairs are automatically generated: 
## Stating 'S' is the addition of all the factors,
## the model formula is '~(S)^(resolution-1)/2',
## when 'resolution' is odd. 
## When 'resolution' is even, the model  formula is '~(S)^(resolution/2)' and the estimate formula is  '~(S)^(resolution/2)-1'
## RETURN
## When the research is not recursive, an object of class 'listofkeyrings'. Each component is the solutions for a value of the prime.
##  When the research is recursive, an object  of class 'listofdesignkeys'. Each component is a whole solution list across the different primes.
## EXAMPLES
## K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D, nunits=3^3, base=~A+B+C, max.sol=2)
## ## With automatic model generation
## Km <- planor.designkey(factors=c(LETTERS[1:4], "block"),nlevels=rep(2,5),resolution=3,  nunits=2^4)
## -------------------------------------------------------
planor.designkey <- function(
                             ## arguments for planor.factors
                             factors=NULL, nlevels=NULL,
                             block=NULL,
                             ordered=NULL,
                             hierarchy=NULL,
                             ## arguments for planor.model
                             model=NULL, estimate=NULL,
                             listofmodels=NULL,
                             ## or resolution for automatic generation of the model
                             resolution=NULL,
                             ## other arguments
                             nunits=NULL, base=NULL,
                             max.sol=1, randomsearch=FALSE,
			     verbose=TRUE){
  ## Build factors if required
  if (class(factors) != "designfactors") {
    factors <- planor.factors(factors, nlevels,
                              block, ordered, hierarchy)
  }
  ## Automatic generation of the model
  if (!is.null(resolution)) {
    res <- generate.model(resolution, factors )
    model <- planor.model(res$model, res$estimate, listofmodels)
  } else {
    ## Build the model  if not yet a list
    if (!is.list(model)) {
      model <- planor.model(model, estimate, listofmodels)
    }
  }
  ## case when a factor is in base but not in the models
  if( !is.null(base) ){
    model <- c( list(list(Model=base,Estimate=base)), model)
  }

  ## A. TREATMENT AND BLOCK FACTORS + PSEUDOFACTORS
  ## Initial step on 'factors', to guarantee coherence between the
  ## arguments. Also stores information on factors that will be required
  ## in the sequel.
  factors <- planor.harmonize(factors=factors, model=model, base=base)
  ## PATCH 07/07/2012: The following is necessary to handle the case when a prime
  ## is in the units decomposition but not represented in the base or model formulae.
  ## In that case, it adds DUMMY factors.
  FACT.primes <- unique(factors@pseudo.info$nlev)
  UNITS.primes <- unique(factorize(nunits))
  pb <- seq(UNITS.primes)[! (UNITS.primes %in% FACT.primes)]
  for(k in pb){
    newfactor <- paste("DUMMY", UNITS.primes[k], sep="")
    factors <- bind(factors, planor.factors(factors=newfactor,
                                            nlevels=UNITS.primes[k],
                                            dummy=TRUE))
    newform <- as.formula(paste("~",newfactor))
    model <- c( list(list(Model=newform, Estimate=newform)), model)
  }
  factors <- planor.harmonize(factors=factors, model=model, base=base)
  ## END OF PATCH

  ## information on the factors and on the pseudofactors.
  F.info <- factors@fact.info
  P.info <- factors@pseudo.info
  FACT.primes <- unique(P.info$nlev)
  BASIC.info <- P.info[P.info$basic,,drop=FALSE]
  MODEL.info <- P.info[P.info$model,,drop=FALSE]
  ## Tests relative to the hierarchy relationships (BEGIN)
  ## Test 1: the following lines stop if the nesting factor of a hierarchy
  ## formula belongs to the basic factors. It should be made possible
  ## to do that in the future but at present it is not compatible with
  ## the algorithm.  Indeed, in the backtrack search, nested factors
  ## must absolutely precede the factor they are nested in. In
  ## addition, the construction of the predefined columns of PhiStar
  ## does not take account of the hierarchies at all. This is ok only
  ## if there is no nesting factor among the basic ones.
  if( 2 %in% BASIC.info$hiera ){
    stop("Sorry, but no factor can be both a basic factor and the nesting factor in a hierarchy formula. Advice: use the nested factor(s) in the base argument instead of the nesting one") }
  ## Test 2: the following lines test if, in each hierarchy, the number of
  ## levels of the nested factor is a multiple of the number of levels
  ## of the nesting one
  if( !is.null(MODEL.info$hiera) ){
    for( hh in seq(ncol(MODEL.info$hiera)) ){
      nlev.nested <- prod( MODEL.info$nlev[ MODEL.info$hiera[,hh]==1 ] )
      nlev.nesting <- prod( MODEL.info$nlev[ MODEL.info$hiera[,hh]==2 ] )
      if( nlev.nested/nlev.nesting != floor(nlev.nested/nlev.nesting) ){
        stop("In one hierarchy, the nested and nesting factors do not have coherent numbers of levels. This error may happen when one the declared factors is not in the model nor among the base factors.")
      }}}
  ## Tests relative to the hierarchy relationships (END)

  ## B. UNITS FACTORS AND THEIR DECOMPOSITIONS
  UNITS.decomp <- factorize(nunits)
  UNITS.primes <- unique(UNITS.decomp)
  UNITS.primes <- UNITS.primes[order(UNITS.primes)]
  ## check for coherence between primes for units and primes for factors
  if( !all(FACT.primes %in% UNITS.primes) ){
    stop("At least one prime number involved in the numbers of levels of the factors is missing from the number of units")
  }
  ## units quasifactors (<=> distinct primes)
  ## there are 'UQUASI.N' distinct primes whose values are in 'UNITS.primes'
  UQUASI.N <- length(UNITS.primes)
  UQUASI.npseudo <- apply( outer(UNITS.primes, UNITS.decomp, "=="), 1, sum )
  UPSEUDO.N <- sum(UQUASI.npseudo) # (= length(UNITS.decomp))
  ## some work needed on basic factors
  UQUASI.nbasic <- apply( outer(UNITS.primes, BASIC.info$nlev, "=="), 1, sum )
  
  ## C. INELIGIBLE FACTORIAL AND PSEUDOFACTORIAL TERMS
  modelterm.matrices <- planor.modelterms(model)
  if(verbose){
    cat("Preliminary step 1 : processing the model specifications\n")}
  b.F.ineligible <- planor.ineligibleterms(modelterm.matrices)
  if(verbose){
    cat("Preliminary step 2 : performing prime decompositions on the factors\n")}
  b.P.ineligible <- planor.ineligibleset(factors, b.F.ineligible)
  ## Free memory
  if (FREE) {
    rm(modelterm.matrices); gc()
  }

  ## ineligible pseudofactorial terms associated with each prime
   b.PTERMS.q <- matrix(NA, nrow=UQUASI.N,
                           ncol=ncol( b.P.ineligible))

  zzz <- seq_len(UQUASI.N)
  for(k in zzz){
    b.PTERMS.q[k,] <-
      0 < apply(b.P.ineligible[MODEL.info$nlev == UNITS.primes[k], ,
                               drop=FALSE], 2, sum)
  }

  ## D. CALCULATION OF THE KEY MATRICES FOR EACH PRIME
  ##    A quite basic loop except when at least one ineligible element
  ##    has several non zerow Sylow components (Kobilinsky 2000, Sections 7 and 8)
  ##    In that case, the option here is as follows:
  ##     - the search order is on increasing primes
  ##     - for each prime p_k, all ineligible terms with a p_k Sylow component are
  ##       included in the ineligible terms for prime p_k
  ##     - if the search fails, a new one is started without
  ## E0. Two possible cases:
  ## Check if some ineligible elements relate to more than one Sylow subgroup,
  ## that is, if some ineligible term is associated with more than one prime

  indpdt.searches <- all(apply(b.PTERMS.q, 2, sum) == 1)

  ## E1. Preparation: a first loop to prepare needed objects
  ## indexed by the p_k. Essentially, construction
  ## of the predefined columns associated with the basic
  ## factors. This part is quite rough and does not take account at
  ## all of the hierarchy formulae. This is why nesting factors are
  ## presently not allowed among the basic factors
  PREDEF <- vector("list", length=UQUASI.N)
  names(PREDEF) <- UNITS.primes
  for(k in zzz){
    ## preparation of the predefined part (only basic factors at the moment)
    r.k <- UQUASI.npseudo[k]            ## r.k units factors at p.k levels ...
    nbase.k <- max(1, UQUASI.nbasic[k]) ## ... including nbase.k basic factors
    PREDEF[[k]] <- diag(x=1, nrow=r.k, ncol=r.k)[, seq(nbase.k),drop=FALSE]
  }
  ## E2. Second loop to search for the design key matrices
  ## E2a. First case: independent searches over the primes
  PHISTAR <- vector("list", length=UQUASI.N)
  names(PHISTAR) <- UNITS.primes

  if(indpdt.searches){
    ## Loop on the primes
    for(k in zzz){
      if(verbose){ cat( paste("Main step for prime p =", UNITS.primes[k],
                              ": key-matrix search\n") ) }
      p.k <- UNITS.primes[k]
      r.k <- UQUASI.npseudo[k]
      pseudos.k <- MODEL.info$nlev == p.k
      ## test for the absence of a factor or term with a multiple of p.k levels
      if (all(pseudos.k==0)) {
        stop("\n the number of units is a multiple of ", p.k, ". ",
             "So, for technical reasons, there should be at least one factor with",
             " a multiple of ", p.k," levels in the model or base formula.",
             "This problem can be solved by adding a replication factor at e.g. ",
             p.k^r.k,
             " levels in the factors argument and in the formula of the base argument.\n")
      }
      else if(all(b.PTERMS.q[k,]==0)){
        if(verbose){ cat("  no factorial term involved with this prime\n") }
      }
      ## extract from hiera only what relates to prime p.k (BEGIN)
      if("hiera" %in% names(P.info)){
        hiera.k <- P.info$hiera[pseudos.k,,drop=FALSE]
        hiera.ok <- apply(hiera.k, 2, sum) > 0
        if(any(hiera.ok)){
          hiera.k <- hiera.k[,hiera.ok,drop=FALSE]
        }
        else{ hiera.k <- NULL }
      }
      else{ hiera.k <- NULL }
      ## extract (END)
      PHISTAR[[k]] <-
        planor.designkey.basep(p=p.k,
                               r=r.k,
                               b.ineligible=
                                 b.P.ineligible[pseudos.k,
                                                as.logical( b.PTERMS.q[k,]),
                                                drop=FALSE ],
                               hierarchy=hiera.k,
                               predefined=PREDEF[[k]],
                               max.sol=max.sol,
                               randomsearch=randomsearch,
                               verbose=verbose)
    } ## end k
  } ## end indpdt.searches
  ## E2b. Second case: recursive searches required
  else{
    PHISTAR <- planor.designkey.recursive(k = 1,
                                          nb.sol=0,
                                          PVuqf = UNITS.primes,
                                          NPuqf = UQUASI.npseudo,
                                          b.INELIGtpf = b.P.ineligible,
                                          NIVtpf = MODEL.info$nlev,
                                          b.INELIGuqf = b.PTERMS.q,
                                          PREDEF = PREDEF,
                                          max.sol = max.sol)
  } ## end recursive search
  ## FINALISATION OF THE SOLUTION OBJECTS
  ## PHISTAR is empty when there is no solution for any k
  if (length(PHISTAR) > 0) {
    ## PROPER NAMING OF THE COLUMNS OF THE KEY MATRICES
    ## loop on the distinct primes
    for(k in zzz){
      p.k <- UNITS.primes[k]
      r.k <- UQUASI.npseudo[k]
      LIBtpf.k <- rownames(P.info)[(P.info$nlev == p.k) & (P.info$model)]
      LIBupf.k <- rownames(BASIC.info)[(BASIC.info$nlev == p.k) & (BASIC.info$model)]
      LIBupf.k <- c(LIBupf.k, rep("*U*", r.k - length(LIBupf.k)))
      ## FIRST CASE: independent searches
      if(indpdt.searches){
        ## test on existence of solutions (PHISTAR[[k]] is null when none)
        if(!is.null(PHISTAR[[k]])){
          PHISTAR[[k]] <- lapply(PHISTAR[[k]], function(X,LIBtpf.k,LIBupf.k ){
              # when X has one single element, it should be transformed into a matrix
        X <- as.matrix(X)
            colnames(X) <- LIBtpf.k
            rownames(X) <- LIBupf.k
            X <- new("keymatrix", .Data=X, p=p.k)
            return(X)},
                                 LIBtpf.k,LIBupf.k)
          PHISTAR[[k]] <- new("keyring",
                              .Data=PHISTAR[[k]] ,
                              p=p.k ,
                              pseudo.info=P.info,
                              LIB=list(LIBtpf.k, LIBupf.k))
        }
      } ## end of indpdt.searches (step FINALISATION)
      ## SECOND CASE: recursive search
      else{
        ## loop on the solution designkeys
        for(ll in seq(length(PHISTAR))){
          names(PHISTAR[[ll]])[k] <- p.k
          colnames(PHISTAR[[ll]][[k]]) <- LIBtpf.k
          rownames(PHISTAR[[ll]][[k]]) <- LIBupf.k
          PHISTAR[[ll]][[k]] <- new("keymatrix", .Data=PHISTAR[[ll]][[k]],
                                    p=p.k)
        } ## end of the loop on the solution designkeys
      } ## end of the recursive case (step FINALISATION)
    } ## end of the loop on the distinct primes

    ## CLASS DECLARATION
    ## FIRST CASE: independent searches
    if(indpdt.searches){

      
       PHISTAR <- new("listofkeyrings",
                       .Data=PHISTAR,
                       factors=factors,
                       nunits=nunits,
                       model=model)
      ##}
     } # fin (indpdt.searches)
    ## SECOND CASE: recursive search
    else{
      for(ll in seq(length(PHISTAR))){
        PHISTAR[[ll]] = new("designkey", .Data=PHISTAR[[ll]] ,
                 factors=factors,
                 nunits=nunits,
                 model=model,
                 recursive=TRUE)
      }
      PHISTAR <- new("listofdesignkeys",
                     .Data= PHISTAR,
                     factors=factors,
                     nunits=nunits,
                     model=model)
    }
  }  ## end of  "if (length(PHISTAR) > 0)"
  
  else {
    return (NULL)
  }
  

  
  return(PHISTAR)
}  ## end planor.designkey
##---------------------------------------------------------------------------
## THE MOST INTEGRATIVE FUNCTION 05/07/2012
##---------------------------------------------------------------------------
## An integrated function to construct a regular design directly
## See the comments of planor.designkey and planor.design for more information
## EXAMPLE
## mydesign <- regular.design(factors=c("block", LETTERS[1:4]),
##  nlevels=rep(3,5), model=~block + (A+B+C+D)^2, estimate=~A+B+C+D,
##  nunits=3^3, randomize=~block/UNITS)
## -------------------------------------------------------------
regular.design <- function(
                    ## arguments for planor.factors
                    factors=NULL, nlevels=NULL,
                    block=NULL, ordered=NULL, hierarchy=NULL,
                    ## arguments for planor.model
                    model=NULL, estimate=NULL, listofmodels=NULL,
                    ## or resolution for automatic generation of the model
                    resolution=NULL,
                    ## other arguments
                    nunits=NULL, base=NULL,
                    ## design stage
                    randomize=NULL,
                    randomsearch=FALSE,
                    ## information
                    output="planordesign",
                    verbose=FALSE,
                    ...){
  ## SEARCH FOR THE KEY MATRIX
  temp <- planor.designkey(factors=factors, nlevels=nlevels, block=block,
                           ordered=ordered, hierarchy=hierarchy,
                           model=model, estimate=estimate, listofmodels= listofmodels,
                           resolution=resolution,
                           nunits=nunits, base=base,
                           max.sol=1,
                           randomsearch=randomsearch, verbose=verbose)

   
  ## DESIGN GENERATION
  
  
  if (is.null(temp) || any( unlist(lapply(temp, is.null)) )){
    cat( paste("No solution\n") )
    return(NULL)
  }
  final <- planor.design(temp, ,randomize=randomize, ...)
  ## FINALISATION
  ## take out the dummy variables
  final <- final[ , !final@factors@fact.info$dummy ]
  ## keep only the data.frame if required
  if(output=="data.frame"){
    final <- getDesign(final)
  }

  return(final)
} # end regular.design
##---------------------------------------------------------------------------
## 2. MAJOR CALCULUS FUNCTIONS
##---------------------------------------------------------------------------
planor.hierarchy <- function(factors, hierarchy=NULL){
  ## Calculates and stores information on the pseudofactor decomposition
  ## of the factors within the 'fact.info' slot of a 'designfactors' object
  ## Called by planor.factors (internal function)
  ## ARGUMENTS
  ##  factors : an object of class 'designfactors'
  ##  hierarchy : a formula or a list of formulae to indicate hierarchy
  ##              relationships between factors
  ## RETURN
  ##  an object of class 'designfactors'
  ## DETAILS:
  ##  This is mainly a function to be used within other functions. It can
  ##  be used to create or renew hierarchy relationships in a 'designfactors'
  ##  object.
  ## EXAMPLES
  ##  F4v6 <- planor.factors( factors=LETTERS[1:4], nlevels=rep(6,4))
  ##  planor.hierarchy(F4v6, hierarchy=list(C~A, D~B )) # internal
## -----------------------------------------------------------
  FACT.names <- rownames(factors@fact.info)
  FACT.N <- length(FACT.names)
  if(is.null(hierarchy)){
    hiera.mat <- matrix(0, FACT.N, 0)   }
  else{
    if(!is.list(hierarchy)){ hierarchy <- list(hierarchy) }
    hiera.mat <- matrix(0, FACT.N, length(hierarchy))
    for(i in seq_along(hierarchy)){
      nesting.i <- all.vars(hierarchy[[i]])[1]
      nested.i <- all.vars(hierarchy[[i]])[-1]
      hiera.mat[ FACT.names == nesting.i, i ] <- 2
      hiera.mat[ FACT.names %in% nested.i, i ] <- 1
    }
  }

  factors@fact.info <- data.frame(factors@fact.info, hiera=I(hiera.mat))

  ## storage of 'pseudofactors' information
  if(nrow(factors@pseudo.info) != 0){
    factors@pseudo.info <- planor.pseudofactors(factors)@pseudo.info
  }
  ##
  return(factors)
}
##---------------------------------------------------------------------------
planor.pseudofactors <- function(factors){
    ## Calculates and stores information on the pseudofactor decomposition
    ## of the factors within a 'designfactors' object
    ## ARGUMENTS
    ##  factors : an object of class 'designfactors'
    ## RETURN
    ##  the 'designfactors' object with the pseudofactors slot completed or updated
    ## DETAILS
    ##  This is mainly a function to be used within other functions. It should
    ##  be ran for updating each time a 'designfactors' object is modified.
    ## EXAMPLES
    ##  F2 <- planor.factors( factors=c(LETTERS[1:4], "block"), nlevels=c(6,6,4,2,6) )
    ##  M2 <- planor.model( model=~block+(A+B+C+D)^2, estimate=~A+B+C+D )
    ##  F2h <- planor.harmonize(factors=F2[,1:5], model=M2, base=~A+B+D)
    ##  attributes(planor.pseudofactors(F2h)) # internal

    ## ---------------------------------------------------------------
    fact.info <- factors@fact.info
    FACT.N <- nrow(fact.info)
    ## prime decompositions
    FACT.npseudo <- rep(NA, FACT.N) ## nbrs of pseudofactor 'children' per factor
    PSEUDO.nlev <- NULL             ## prime nbrs of levels of each pseudofactor

    for(i in seq_len(FACT.N)){
        split.s <- factorize(fact.info$nlev[i])
        FACT.npseudo[i] <- length(split.s)
        PSEUDO.nlev <- c(PSEUDO.nlev, split.s)
    }
    ## pseudofactor names
    PSEUDO.parent <- rep(seq_len(FACT.N), FACT.npseudo) ## 'parents' of the pseudofactors
    numero <- unlist( sapply(FACT.npseudo , seq_len) )
    PSEUDO.num <- paste("_", as.character(numero), sep="")
    PSEUDO.num[ (FACT.npseudo == 1)[PSEUDO.parent] ] <- ""
    PSEUDO.names <- paste( names(factors)[PSEUDO.parent], PSEUDO.num, sep="" )

    ## storage of 'pseudofactors' information
    pseudo.info <- data.frame(parent = PSEUDO.parent,
                              nlev = PSEUDO.nlev,
                              stringsAsFactors=FALSE)
    pseudo.info <- cbind(pseudo.info,
                         fact.info[PSEUDO.parent,
                                   !(names(fact.info) %in% c("names","nlev"))])
    rownames(pseudo.info) <- PSEUDO.names
    factors@pseudo.info <- pseudo.info
    ##
    return(factors)
}
##---------------------------------------------------------------------------
## "planor.harmonize"  
## Harmonize the factors originating from a list of factors, a list of models, and a list of basic factors
## ARGUMENTS
##  - factors: can be of 2 types:
## either,  an object of class 'designfactors' (typically an output
##               from 'planor.factors'),
## or a character vector of factor names.
## In this last case, the arguments 'nlevels', 'ordered', 'hierarchy' can also be set.
##   - nlevels:  a vector of level numbers for each factor name.
## Ignored if 'factors' is an object of class 'designfactors'.
##   - ordered:   an additive model formula to indicate the ordered factors
## Ignored if 'factors' is an object of class 'designfactors'.
##   - hierarchy:   a formula or a list of formulae to indicate hierarchy relationships between factors
## Ignored if 'factors' is an object of class 'designfactors'.
## - model: either a list of model-estimate pairs of formulae, typically an output
##               from 'planor.model', or
##  the   formula of the main model. In this last case, the arguments
## 'estimate' and 'listofmodels' can also be set.
##   - estimate:  if 'model' is a formula,
## optional formula specifying the factorial terms to estimate;
##               if missing, it is considered that all model terms have to be
##               estimated.
## Ignored if 'model' is a list.
##   - listofmodels:  if 'model' is a formula,
## list of c(model, estimate) pairs, where model and estimate are
##              formulae; using several pairs allows more flexibility in the design
##              constraints  (see Kobilinsky, 2005); estimate is optional.
## Ignored if 'model' is a list.
##   -  base: an optional formula to specify the basic factors. These factors  must belong to the factors argument
## RETURN
##   An object of class 'designfactors' very similar to 'factors', but with two additional columns in slots  'fact.info' and 'pseudo.info':
##
##  - 'model' (logical, TRUE for factors present in some formula)
##
##  - 'basic' (logical, TRUE for basic factors)
## NOTE
##     This function is essentially a check that the factors in all three arguments
##     are coherent, even though it performs some additional work.
##     The function stops if it detects a model or basic factor that is absent from
##     'factors'. This is because the number of levels of such a
##     factor is unknown and so the design search cannot proceed.
##     Besides, the function eliminates the factors that do appear neither in
##     'model' nor in 'base' and it reorders the factors by putting first the basic  ones.
## EXAMPLES
##     F2 <- planor.factors( factors=c(LETTERS[1:4], "block"), nlevels=c(6,6,4,2,6) )
##     M2 <- planor.model( model=~block+(A+B+C+D)^2, estimate=~A+B+C+D )
##     planor.harmonize(factors=F2[,1:5], model=M2,base=~A+B+D)
## -----------------------------------------------------
planor.harmonize <- function(
                             ## arguments for planor.factors
                             factors=NULL, nlevels=NULL,
                             ordered=NULL,
                             hierarchy=NULL,
                             ## arguments for planor.model
                             model=NULL, estimate=NULL,
                             listofmodels=NULL,
                             ## other arguments
                             base=NULL){
  ## Build factors if required
  if (class(factors) != "designfactors") {
    factors <- planor.factors(factors, nlevels,
                              ordered,hierarchy)
  }
  ## Build the model  if required
  if (!is.list(model)) {
    model <- planor.model(model, estimate, listofmodels)
  }

  ## Names of the factors in the 'factors' argument
  factors.names <- rownames(factors@fact.info)
  ## Names of the factors in the 'base' argument
  
  if(is.null(base)) base.names <- vector(mode="character",length=0)
  else base.names <- attr(terms(base),"term.labels")
  
  ## Names of the factors in the model argument
  ## remark: the names are searched for only in the right side of the formulae
  model.names <- vector("character",0)
  for(i in seq_along(model)){
    model.names <-
      unique(c(model.names,
               all.vars( rev(as.list(model[[i]][[1]]))[[1]] ),  # model part
               all.vars( rev(as.list(model[[i]][[2]]))[[1]] ))) # estimate part
  }
  ## Continue only if all factors in the models and in the base list are declared
  ## in the factors argument (otherwise we cannot know their numbers of levels)
  missing.check <- ! (unique(c(model.names, base.names)) %in% factors.names)
  if( any(missing.check) ){
    missing.names <- unique(c(model.names, base.names))[missing.check]
    stop(paste("\n*****\n",
               "Model or basic factor(s) <", missing.names,
               "> not specified in the factors argument.\n",
               "Consequence: unknown number(s) of levels.\n",
               "*****", sep=""))
  }

  ## Classification of the factors
  factors.type <-
    3 - 2*(factors.names %in% model.names) - 1*(factors.names %in% base.names)
  ## factors.type equals:
  ## 0 for a factor in the model and basic factors
  ## 1 for a factor in the model only
  ## 2 for a factor in the basic factors only
  ## 3 for a factor outside the model and basic factors

  ## Store the information



  factors@fact.info$model <- factors.type <= 1
  factors@fact.info$basic <- (factors.type == 0) | (factors.type == 2)
  ## Keep only the factors of type 0, 1, or 2
    factors <- factors[ factors.type < 3 ]
  
  ## Re-order to put basic factors first
  if(! is.null(base)){
    is.basic <- factors@fact.info$basic
    ## Application of the method "bind" of designfactors
    if(!all(is.basic)){
      factors <- bind( factors[is.basic], factors[!is.basic] )
    }
  }
  

  return(factors)
}
##---------------------------------------------------------------------------
planor.modelterms <- function(modlist){
    ## Creates the 0-1 matrix pairs of model and estimate terms.
    ## Internal function
    ## ARGUMENTS
    ##  modlist: list of c(model, estimate) pairs, where model and estimate are
    ##           formulae; typically, an output from planor.model
    ## RETURN
    ##  a list of c(model, estimate) pairs,
    ##  where model and estimate are 0-1 matrices indexed by factors in rows and
    ##  factorial terms in columns
    ## DETAILS
    ##  the factors indexing rows are those appearing in the model and estimate formulae
    ## EXAMPLES
    ##  modelexample <- planor.model(~A+B+C+D+A:B,~A+B+C+D)
    ##  planor.modelterms(modelexample)
    ##  M2 <- planor.model( model=~block+(A+B+C+D)^2, estimate=~A+B+C+D )
    ##  planor.modelterms(M2)
  ## ------------------------------------------------------

    modmat <- vector("list", length=length(modlist))
    for(i in seq_along(modlist)){
        ## step i, matrix of model terms
        Mterms.i <- attributes(terms(modlist[[i]][[1]]))
        ModelTerms.i <- Mterms.i$factors
        NullModel.i <- length(ModelTerms.i) == 0 # case when model=~1
        if( !NullModel.i ){
          if(Mterms.i$response > 0){ ModelTerms.i <- ModelTerms.i[-1,,drop=FALSE] }
          ## check for marginality completeness
          if(max(ModelTerms.i)>1){ stop("Sorry, nesting in the model formulae not possible") }
          ModelTerms.i <- cbind(MU=0, ModelTerms.i)
        }
        ## step i, matrix of estimate terms
        Eterms.i <- attributes(terms(modlist[[i]][[2]]))
        EstimTerms.i <- Eterms.i$factors
        NullEstim.i <- length(EstimTerms.i) == 0 # case when estim=~1
        if( !NullEstim.i ){
          if(Eterms.i$response > 0){ EstimTerms.i <- EstimTerms.i[-1,,drop=FALSE] }
          ## check on marginality completeness
          if(max(EstimTerms.i) > 1){
            warning(paste("The estimate formula ", i,
                          " contains at least one interaction without all its marginal terms.\n",
                          " Such marginal terms may be unestimable in the final design.\n", sep=""))
            EstimTerms.i[EstimTerms.i == 2] <- 1
          }
        }
        ## step i, four different cases for the model and estimate parts
        ## first case: model=~1 and estimate=~1
        if(NullModel.i & NullEstim.i){
          ModelFine.i <- NULL
          EstimFine.i <- NULL
        }
        ## second case: model=~1 and estimate NOT=~1
        else if(NullModel.i){
          ModelFine.i <- NULL
          EstimFine.i <- EstimTerms.i
        }
        ## third case: model NOT=~1 and estimate=~1
        else if(NullEstim.i){
          ModelFine.i <- ModelTerms.i
          EstimFine.i <- NULL
          }
        ## fourth case: model NOT=~1 and estimate NOT=~1
        else{
          ## making of coherent rows between model and estimate matrices
          mLIBtf.i <- rownames(ModelTerms.i)
          eLIBtf.i <- rownames(EstimTerms.i)
          LIBtf.i <- unique( c(mLIBtf.i, eLIBtf.i) )
          Ntf.i <- length(LIBtf.i)
          Nmterms.i <- ncol(ModelTerms.i)
          ModelFine.i <- matrix(0, Ntf.i, Nmterms.i,
                                dimnames=list(LIBtf.i, colnames(ModelTerms.i)))
          ModelFine.i[mLIBtf.i,] <- ModelTerms.i[mLIBtf.i,]
          Neterms.i <- ncol(EstimTerms.i)
          EstimFine.i <- matrix(0, Ntf.i, Neterms.i,
                              dimnames=list(LIBtf.i, colnames(EstimTerms.i)))
          EstimFine.i[eLIBtf.i,] <- EstimTerms.i[eLIBtf.i,]
          } ## end of the fourth case

        modmat[[i]] <- list(model=ModelFine.i, estimate=EstimFine.i)
    } ## end of the loop on (model,estimate) pair i
    return(modmat)
}
##---------------------------------------------------------------------------
planor.ineligibleterms <- function(modmat){
    ## Calculates the ineligible factorial terms of a planor-type model
    ## Internal function
    ## ARGUMENTS
    ##  modmat: a list of model-estimate pairs of 0-1 matrices, typically an output
    ##         from planor.modelterms
    ## RETURN
    ##  A 0-1  matrix with each row associated with a factor and each column with
    ##  an ineligible factorial term;
    ## EXAMPLE
    ##  M2 <- planor.model( model=~block+(A+B+C+D)^2, estimate=~A+B+C+D )
    ##  M2terms <- planor.modelterms(M2)
    ##  print(planor.ineligibleterms(M2terms)) # internal
  ## --------------------------------------------------------

    LIBtf <- c()
    for(i in seq_along(modmat)){
        LIBtf <- unique( c(LIBtf, rownames(modmat[[i]]$model), rownames(modmat[[i]]$estimate)) )
    }
    Ntf <- max(1,length(LIBtf))

    ## Calculation of ineligible terms -
    ## Initialisation of b.ineligible to the identity matrix of order
    ## equal to the total number of factors in the formulae. This imposes that
    ## all levels of each factor will be represented in the designs,
    ## whereas it was optional in the APL PLANOR
  
  
  
    b.ineligible <- as.matrix( diag(Ntf))

## End initialisation

    for(i in seq_along(modmat)){
        ModelTerms.i <- modmat[[i]]$model
        EstimTerms.i <- modmat[[i]]$estimate

        ## first case: model=~1 and estimate=~1
        if(is.null(ModelTerms.i) & is.null(EstimTerms.i)){
          stop("The pair (model=~1,estimate=~1) is not allowed.") }
        ## second case: model=~1 and estimate NOT=~1
        else if(is.null(ModelTerms.i)){
          b.ineligible.i <- EstimTerms.i }
        ## third case: model NOT=~1 and estimate=~1
        else if(is.null(EstimTerms.i)){
          b.ineligible.i <- ModelTerms.i[,-1] }
        ## fourth case: model NOT=~1 and estimate NOT=~1

        else{
          b.ineligible.i <- symmdiff(ModelTerms.i, EstimTerms.i) }
        ## adaptation to the whole factor list
        iLIBtf.i <- rownames(b.ineligible.i)
        colnames(b.ineligible.i) <- NULL
        b.ineligibleF.i <-matrix(0, nrow= Ntf, ncol=ncol(b.ineligible.i),
                                      dimnames=list(LIBtf, colnames(b.ineligible.i)))
        b.ineligibleF.i[iLIBtf.i,] <- b.ineligible.i[iLIBtf.i,]
          b.ineligible <- cbind(b.ineligible, b.ineligibleF.i)

        ## elimination of redundant columns
        ## a rough way to avoid last rows to disappear
          b.ineligible <-  rbind(b.ineligible,1)

         # the second argument of convertinto.basep is 2 to get a 0-1 matrix
          ineligible.codes <- unique( convertfrom.basep(t(b.ineligible),2) )
          b.ineligible <- t(convertinto.basep(ineligible.codes,2))
        b.ineligible <- as.matrix(b.ineligible[-(Ntf+1),,drop=FALSE]) ## rough way again
    }
    rownames(b.ineligible) <- LIBtf

    return(b.ineligible)
} ## end planor.ineligibleterms
##---------------------------------------------------------------------------
planor.ineligibleset <- function(factors, b.ineligible){
    ## Construction of the set of ineligible pseudofactorial effects from a
    ## set of ineligible model terms.
    ## Internal function
    ## ARGUMENTS
    ##  - factors: a 'designfactors' object
    ##  - ineligible: an output from planor.ineligibleterms, that is, a 0-1 matrix
    ##              with one row per factor and one column per ineligible term
    ## OUTPUT
    ##  A 0-1 matrix with s rows associated with pseudofactors and each column
    ##  an ineligible pseudofactorial term
    ## DETAILS
    ##  The ineligible effects derive from the ineligible terms through the
    ##  decomposition of factors into pseudofactors with a prime number of levels.
    ##  This function includes three main steps:
    ##  1) splitting of the factors and factorial terms into "quasifactors" and
    ##     "quasi-factorial terms", based on decompositions of level numbers into
    ##     powers of distinct primes
    ##  2) first reduction of the ineligible terms (see Kobilinsky 2000)
    ##  3) splitting of the quasifactors and quasifactorial terms into "pseudofactors"
    ##     and "pseudofactorial terms", based on decompositions of level numbers
    ##     that are powers of a prime into primes
    ## EXAMPLE
    ##  Fset <- planor.factors( factors=c(LETTERS[1:3]), nlevels=c(6,6,4) )
    ##  Iset <- cbind(diag(3), c(1,1,1))
    ##  rownames(Iset) <- c("A","B","C")
    ##  print(planor.ineligibleset(Fset, Iset)) # internal
    ## --------------------------------------------------------------
    
    ## restriction to the factors present in the rows of 'ineligible'
    F.select <- names(factors) %in% rownames(b.ineligible)
    factors <- factors[ F.select ]
    FACT.info <- factors@fact.info
    PSEUDO.info <- factors@pseudo.info
    ## reordering of the rows of 'ineligible', in case
        b.ineligible <- as.matrix(b.ineligible[rownames(FACT.info),,drop=FALSE])

    
    ## NOW, FACT.info, PSEUDO.info and ineligible correspond to the same factors
    ## in the same order
    
    ## Number of factors
    FACT.N <- nrow(FACT.info)
    ## Number of factorial terms in 'ineligible'
    FTERMS.N <- ncol(b.ineligible)
    ## Info on factors and pseudofactors
    FACT.nlev <- FACT.info$nlev         ## level nbrs
    PSEUDO.parent <- PSEUDO.info$parent    ## parents
    PSEUDO.nlev <- PSEUDO.info$nlev    ## level nbrs
    ## distinct primes
    PRIMES <- unique(PSEUDO.nlev)
    ## nbrs of pseudofactors per factor
    FACT.npseudo <- FACT.info$npseudo
    
    ## QUASIFACTORS:
    ## As defined here, there is one quasifactor per factor and per
    ## distinct prime in its decomposition. For example, a factor at 12
    ## levels has 2 quasifactors at 4 and 3 levels and 3 pseudofactors
    ## at 2, 2, 3 levels. So quasifactors are an intermediate level of
    ## decomposition between factors and pseudofactors. This is useful
    ## for the 'first split' described below.
    ##
    ## nbr of distinct primes per factor = nbr of quasifactors per factor
    FACT.nquasi <- tapply( PSEUDO.nlev, PSEUDO.parent, function(x){length(unique(x))} )
    FACT.cquasi <- c(0, cumsum(FACT.nquasi))
    QUASI.N <- sum(FACT.nquasi)
    ## *LIST* of the decomposed numbers of levels of each factor
    NLEV.decomp <- tapply( PSEUDO.nlev, PSEUDO.parent, list )
    ## *UNLIST* of the distinct primes of each factor
    QUASI.prime <- unlist( lapply( NLEV.decomp, unique ) )
    ## *UNLIST*  of the associated exponents = nbr of pseudofactors per factor
    QUASI.npseudo <- unlist( lapply( NLEV.decomp,
                                    function(x){apply(outer(unique(x),x,"=="),1,sum)}) )
    
    ## 1) First split: factor -> one (temporary) "quasi"-factor per prime component
    ## that step is not in Kobilinsky 2000. It fulfills part of the reduction
    ## algorithm p.10 directly by avoiding the creation of interactions between
    ## quasifactors associated with the same initial factor.
    ## Matrix "ineligible" is turned into matrix "iqtI" with one row per
    ## quasifactor and one column per "quasi" ineligible term, that is, a term
    ## obtained by splitting the component factors of an initial factorial term
    ## into their associated quasifactors
    
    ## FTERMS.nquasi : numbers of quasifactorial terms per ineligible factorial term
    ##               (vector of length FTERMS.N)
        FTERMS.nquasi <- apply( diag(FACT.nquasi) %*% b.ineligible, 2, function(x) prod(x[x>0])  )

    
    FTERMS.cquasi <- c(0, cumsum(FTERMS.nquasi))
    QTERMS.N <- sum(FTERMS.nquasi)
    ## initialisation of iqtI and loop on the ineligible terms that must be split
    ## iqtI: ineligible quasifactorial terms
    iqtI <-  b.ineligible[ rep(seq_len(FACT.N), FACT.nquasi),
                          rep(seq_len(FTERMS.N), FTERMS.nquasi) ]
    
    for(j in seq_len(FTERMS.N)[FTERMS.nquasi > 1]){
        tfactors.j <- seq_len(FACT.N)[ as.logical(b.ineligible[,j]) ]
        Ntf.j <- length(tfactors.j)
        qrows.j <- as.logical( rep( b.ineligible[,j], FACT.nquasi ) )
        qcols.j <- FTERMS.cquasi[j] + seq_len(FTERMS.nquasi[j])
        iqtI[ qrows.j, qcols.j ] <- 0
        
        qterms.i <- t(crossing(FACT.nquasi[tfactors.j])) +
            matrix( FACT.cquasi[tfactors.j], Ntf.j, FTERMS.nquasi[j] )
        qterms.j <- t( matrix(seq_len(FTERMS.nquasi[j]), FTERMS.nquasi[j], Ntf.j) ) +
            FTERMS.cquasi[j]
        iqtI[cbind(c(qterms.i), c(qterms.j))] <- 1
    }
    
    ## 2) First reduction: application of the algorithm on page 10 of Kobilinsky 2000
    ## iqtR: reduced set of ineligible quasifactorial terms
    
    ## QTERMS.nprime : numbers of distinct non-zero primes in each ineligible quasifactorial
    ##                 term (vector of length QTERMS.N)
    

    QTERMS.nprime <- apply( diag(x=QUASI.prime, nrow=length(QUASI.prime)) %*% iqtI, 2,
                           function(x){ length(unique(x[x!=0])) } )
    b.iqtR <- NULL
    mdiff <- diff(range(QTERMS.nprime))

    iqtI <- as.matrix(iqtI)
    
    while( ncol(iqtI)>0 ){
        m <- min(QTERMS.nprime)
        selectm <- QTERMS.nprime == m
            b.iqttq <- as.matrix(iqtI[, selectm, drop=FALSE ])
            b.iqtR <- cbind(b.iqtR, b.iqttq)


        iqtI <- iqtI[, !selectm, drop=FALSE ]
        QTERMS.nprime <- QTERMS.nprime[ !selectm, drop=FALSE ]
        qft.j <- 0
        while( (qft.j < ncol(b.iqttq)) && (ncol(iqtI) > 0) ){
            qft.j <- qft.j + 1
            qrows.j <- rep(FALSE, QUASI.N)
            for(k in unique( QUASI.prime[ as.logical(b.iqttq[,qft.j]) ] )){
                qrows.j <- qrows.j | (QUASI.prime == k)
            }
            iqtK.j <- rep(FALSE, ncol(iqtI))
            for(l in seq_len(ncol(iqtI))){
                iqtK.j[l] <- all( b.iqttq[qrows.j,qft.j] == iqtI[qrows.j,l] )
            }
            iqtI <- as.matrix(iqtI[, !iqtK.j ])
            QTERMS.nprime <- QTERMS.nprime[ !iqtK.j ]
        } # end of the loop on qft.j
        if( (ncol(iqtI)>0) ){
            mdiff <- diff(range(QTERMS.nprime))
            if( mdiff==0 ){
                    b.iqtR <- cbind(b.iqtR, iqtI)
                iqtI <-  iqtI[,0]
            }
        }
    }   ## end of the loop on m
    
    ## 3) Second split: quasifactor -> exponent-many pseudo-factors
    ## Matrix "iqtR" is turned into matrix "iptR" with one row per
    ## pseudofactor and one column per pseudo ineligible term, that is, a term
    ## obtained by splitting the component quasifactors of a quasi factorial term
    ## into their associated pseudofactors
    ## Npterms : total number of pseudofactorial ineligible terms
    QTERMS.newN <- ncol(b.iqtR)
        QTERMS.npseudo <- apply( diag(QUASI.npseudo) %*% b.iqtR, 2, function(x) prod(2^(x[x>1]) - 1)  )

    QTERMS.cpseudo <- c(0, cumsum(QTERMS.npseudo))
    ##  iptR <- iqtR[ rep(seq(QUASI.N), QUASI.npseudo),
    ##                rep(seq(QTERMS.newN), QTERMS.npseudo) ]
        b.iptR <- b.iqtR[ rep(seq(QUASI.N), QUASI.npseudo),
                         rep(seq(QTERMS.newN), QTERMS.npseudo) ]

    for(j in seq_len(QTERMS.newN)[QTERMS.npseudo > 1]){
        prows.j <- as.logical( rep( b.iqtR[,j], QUASI.npseudo ) )
        pcols.j <- QTERMS.cpseudo[j] + seq_len(QTERMS.npseudo[j])
        b.iptR[ prows.j, pcols.j ] <- 0
        tqfactors.j <- seq_len(QUASI.N)[ as.logical(b.iqtR[,j]) ]
        combis.j <- crossing(2^QUASI.npseudo[tqfactors.j]-1)
        submat.j <- NULL
        for(f in seq_len(length(tqfactors.j))){
            submat.j <- cbind(submat.j, convertinto.basep(combis.j[,f],2))
        }
        b.iptR[ prows.j, pcols.j ] <- t( submat.j )
    }
    rownames(b.iptR) <- rownames(PSEUDO.info)
    
    return(b.iptR)
} ## end planor.ineligibleset
##---------------------------------------------------------------------------
planor.designkey.basep <- function(p, r, b.ineligible,
                                   hierarchy, predefined=NULL,
                                   max.sol=1, randomsearch=FALSE,
				   verbose=FALSE){
  ## Searches for a design key matrix when all s treatment factors are at p levels
  ## Internal function
  ## ARGUMENTS
  ##  - p: a prime number
  ##  - r: integer (the power of p that gives the number of units)
  ##  ineligible: a 0-1 matrix of ineligible (pseudo-)factorial terms
  ##  - hierarchy: NULL or a matrix 0-p
  ##  - predefined: if specified, a (r x f) matrix of prespecified defining
  ##              relationships (f smaller than s)
  ##  - max.sol: maximum number of solutions before exit
  ##  - randomsearch: if TRUE, the order of admissible elements is randomised
  ##                at each new visit forward to j
  ##  - verbose: logical TRUE to trace the execution
  ## RETURN
  ##  a list of design key matrices of size r x s with elements in Zp
## -----------------------------------------------------------------
  ## PRELIMINARIES: (a) ineligible characters
  s <- nrow(b.ineligible)
  ## turn the ineligible factorial terms into ineligible p-group characters
  if(p!=2){ b.ineligible <- representative.basep(b.ineligible,p) }

  ## Indices of the last two non-zero coefficients in each ineligible trt character:
  ## - first row of ineligible.lnz: indices of the last but one non-zero coefficients;
  ## - second row of ineligible.lnz: indices of the last non-zero coefficients.
  ## The columns of ineligible and ineligible.lnz are then reordered
  ineligible.lnz <- apply(b.ineligible, 2,
                             function(x){
                               nz <- seq_along(x)[x!=0]
                               if(length(nz) < 2) maxis <- c(0,max(nz))
                               else maxis <- nz[length(nz)+c(-1,0)]
                               return(maxis)})

  ineligible.reorder <- order(ineligible.lnz[2,],ineligible.lnz[1,])
  ##  ineligible <- ineligible[ , ineligible.reorder]
     b.ineligible <-  b.ineligible[ , ineligible.reorder]

  ineligible.lnz <- ineligible.lnz[ , ineligible.reorder]
  if (FREE) {
    rm(ineligible.reorder); gc()
  }

  ## PRELIMINARIES: (b) hierarchies
  ## Uses the convention that nesting factors are associated with a 2 in 'hierarchy'
  ## and nested ones with a 1
  ## Hmaxj : 0 if there is no hierarchy constraint on factor A_j
  ##         else, index of the higher factor that has to be nested in A_j
  ## Hset : list of index vectors of the factors that have to be nested in A_j
  Hset <- vector("list", length=s)
  Hmaxj <- rep(0, length=s)
  if(!is.null(hierarchy)){
    for(j in seq_len(ncol(hierarchy))){
      nesting <- seq_len(s)[hierarchy[,j] == 2]
      for(k in nesting){
        Hset[[k]] <- seq_len(s)[hierarchy[,j] == 1]
        Hmaxj[k] <- max(Hset[[k]])
      }
    }
  }

  ## PRELIMINARIES: (c) predefined part of PhiStar:
  ## if no predefined matrix, the first column of PhiStar is forced
  ## to be equal to (1 0 ... 0)
  if(is.null(predefined)) {
    b.PhiStar <- matrix(c(1,rep(0,r-1)), r, 1)
  }  else {
     b.PhiStar <- as.matrix(predefined)
   }

  f <- ncol(b.PhiStar)

  if(s < f)
    stop("There is likely to be too many base factors for this prime. ")
  if(s == f){
        check <- !any(apply(((b.PhiStar %*% b.ineligible)%%p)==0, 2, all))

    if(check){
      ## check is 1 when a column is entirely null
      if(verbose){ cat("  no need (all columns are predefined)\n") }
      return(list(matrix(b.PhiStar[,], ncol=ncol(b.PhiStar))))
    }
    else stop("Design key completely predefined but inadequate. Check the base factors.")
  } ## end (s == f)
  if(verbose){
    if((s-f)==1){ cat("  => search for column",s,".\n") }
    else{ cat("  => search for columns",f+1,"to",s,"\n") }
  }

  if(s > f){
    Hset <- Hset[f+seq(1:(s-f))]
    Hmaxj <- Hmaxj[f+seq(1:(s-f))]
  }


  ## INITIALISATIONS FOR THE BACKTRACK SEARCH
  ## Ustar = non-zero elements of U* = initial candidate columns for PhiStar
  U.nb <- (p^r)-1
  b.Ustar <- t(convertinto.basep(seq_len(U.nb),p))

  ## initially admissible elements for j [(r x n_j) p-matrices]
  b.iaA <- vector("list", length=s-f)
  ## status of initially admissible elements for j [n_j-length vectors]
  siaA <- vector("list", length=s-f)
  ## currently admissible elements for j [vectors of selected indices between 1 and n_j]
  wiaA <- vector("list", length=s-f)
  ## smallest visited index k since the previous visit in j
  ktr <- rep(0,s-f)
  ## nbs of initially admissible elements for j
  liaA <- rep(0,s-f)

  ## BACKTRACK SEARCH
  ## Remarks:
  ## 1. When j is reached forward, the admissible elements of U*, for
  ##    column f+j  of the key matrix, must be identified. We assume that the ineligible
  ##    elements are all in canonical form, that is, their last non-zero coordinate is
  ##    equal to 1. So we do not need it hence the selection of their first (f+j-1)
  ##    coordinates at step j (forward) below (see Kobilinsky, 2000, end of page 16)
  ## 2. The management of initially admissible elements and hierarchies is different
  ##    from that proposed in Kobilinsky, 2000. Here, new 'initially admissible'
  ##    elements are calculated as soon as hierarchy constraints have changed for
  ##    a given factor A_j.
  PhiStar.solution <- NULL
  jprev <- 0  ;  j <- 1
  while(j > 0){
    b.PhiStar <- b.PhiStar[,seq(f+j-1), drop=FALSE]

    ## STEP 1: IF column j is reached forward
    ## -> identification and updating of the admissible elements
    if(jprev < j){
      ## STEP 1.A: generation of the INITIALLY admissible elements (iae) for j
      ##  -> performed ONLY IF this is the first visit to j with the present
      ##     set of Hset[[j]] columns
      ##  -> this part will remain valid during all the backtrack search except
      ##     for factors A_j subject to changing hierarchy constraints
      ##     (Hmaxj[j] > f)
      ## ---
      ## The condition just below is satisfied if this is the first visit to j
      ## or if A_j is subject to hierarchy constraints that have changed since
      ## the previous visit to j
      if(ktr[j] <= Hmaxj[j]){
        if(verbose && (ktr[j] == 0)){ cat("      first visit to column", f+j,"\n") }
        ## 1.A.a: elements of U* that satisfy the hierarchy constraints
        ## ---
        ## The IF condition just below is satisfied only if A_j is subject to
        ## hierarchy constraints and they have changed since the previous visit
        if(0 < Hmaxj[j]){
          ## b.aux = b.PhiStar[,Hset[[j]],drop=FALSE]
          b.aux = b.PhiStar[,Hset[[j]],drop=FALSE]
          b.candidates <- subgroup.basep(b.aux, p, all=TRUE)
        }
        ## ELSE -> no hierarchy constraint
        else{ b.candidates <- b.Ustar  }
        ## 1.A.b: selection based on the relevant ineligible subset
        k <- max(f,Hmaxj[j])
        select.kj <- (ineligible.lnz[1,] <= k) & (ineligible.lnz[2,] == (f+j))

        if (all(select.kj==FALSE)){


          iae.j <- rep(TRUE, ncol(b.candidates))
        }
        else{
            b.inelig.kj <- as.matrix(b.ineligible[ seq_len(f+j-1),
                                                    select.kj, drop=FALSE ])
          iae.j <- planor.kernelcheck.basep(b.PhiStar, b.candidates,
                                            b.inelig.kj,p)
        }

        ##        iaA[[j]] <- b.candidates[, iae.j, drop=FALSE]
        if (all(!iae.j)) {
          cat(paste("No solution for column ", f+j, " of the design key\n"))
          return(PhiStar.solution)
        }
          b.iaA[[j]] <- b.candidates[, iae.j, drop=FALSE]
        ## 1.A.c: Information storage after the visit to j
        liaA[j] <- sum(iae.j)
        siaA[[j]] <- rep(s, liaA[j])
        ktr[j] <- k+1
      }

      ## STEP 1.B: UPDATING of the STATUS of the initially admissible elements
      ##           for factor A_j
      ## -> based on the ineligible elements with last-but-one non-zero element
      ##    between ktr[j] and j-1 and with last non-zero element equal to j
      ## ---
      ## The IF condition just below is satisfied when kernel checks are necessary
      if(ktr[j] < f+j){
        ## Make elements eligible if their status may have changed
        admiss.j <- siaA[[j]]
        admiss.j[ admiss.j >= ktr[j] ] <- s
        ## Loop over the last-but-one non-zero indices
        for(k in seq.int(ktr[j], f+j-1)){
          ## relevant ineligible elements
          select.kj <- (ineligible.lnz[1,] == k) & (ineligible.lnz[2,] == (f+j))

          ## relevant initially admissible elements
          tocheck <- seq_len(liaA[j])[admiss.j == s]

          if ( (length(tocheck > 0)) & any(select.kj==TRUE) ){
            b.inelig.kj <- as.matrix(b.ineligible[ seq_len(f+j-1),
                                                      select.kj, drop=FALSE ])
            ae.j <- planor.kernelcheck.basep(b.PhiStar,
                                             b.iaA[[j]][, tocheck,drop=FALSE],
                                             b.inelig.kj,
                                             p)
            admiss.j[ tocheck[!ae.j] ] <- k
          }

        }
        siaA[[j]] <- admiss.j
      }

      ## STEP 1.C: UPDATED admissible elements for factor A_j
      admis.j <- seq_len(liaA[j])[siaA[[j]]==s]
      if(randomsearch){
        wiaA[[j]] <- admis.j[sample(length(admis.j))] }
      else{
        wiaA[[j]] <- admis.j }
      ktr[j] <- f+j
    }
    ## END OF STEP 1
    ##
    ## STEP 2: INCLUSION OF THE NEXT ADMISSIBLE ELEMENT IN PhiStar
    wj <- wiaA[[j]]
    ## CASE 2.A: If there is a next admissible element, add it to PhiStar
    if(length(wj) > 0){
      ##         newcolj <- (iaA[[j]])[,wj[1]]
      b.newcolj <- (b.iaA[[j]])[,wj[1]]
      wiaA[[j]] <- wj[-1]
      b.PhiStar <- cbind(b.PhiStar, b.newcolj, deparse.level=0)

      ## CASE 2.A.a: If this was the last column to fill in PhiStar
      if(j == (s-f)){
        ## Save the solution
        PhiStar.solution <- c(PhiStar.solution, list(b.PhiStar[,,drop=FALSE]))
        ## Return if enough solutions have been found
        if(length(PhiStar.solution) == max.sol){
          cat("The search is closed: max.sol = ",
              length(PhiStar.solution), "solution(s) found \n")
          return(PhiStar.solution)
        }
        ## Otherwise look at the next possibility
        jprev <- j ; j <- j
      }
      ## CASE 2.A.b: If this was not the last column, go forward to next j
      else{
        jprev <- j ; j <- j+1
      }
    }
    ## CASE 2.B: If there is no next admissible element
    else{
      ## Update information
      
      
        ktr.change <- (seq(s-f) >= j) & (ktr > j-1)
        ktr[ ktr.change ] <- j-1

      ## Then go backward
      jprev <- j ; j <- j-1
    }
  }
  cat("The search is closed: ",
      length(PhiStar.solution), "solutions found \n")
  return(PhiStar.solution)
} ## end planor.designkey.basep
##---------------------------------------------------------------------------
planor.designkey.recursive <- function(k,nb.sol,PVuqf,NPuqf,
                                       b.INELIGtpf,NIVtpf,
                                       b.INELIGuqf,PREDEF,
                                       max.sol=1){
  ## Organises the loop over primes (Sylow subgroups) to calculate
  ## the design key submatrices
  ## Internal function
  
  ## ARGUMENTS
  ##  k: index of the prime number to be treated
  ##  nb.sol: current number of solutions found
  ##  PVuqf: vector of the (units) prime numbers
  ##  NPuqf: powers of the prime numbers to get the number of units
  ##  b.INELIGtpf: ineligible pseudofactorial terms
  ##  NIVtpf: vector of the pseudofactors' numbers of levels
  ##  b.INELIGuqf: correspondence matrix between ineligible pseudofactorial terms and primes

  ##  PREDEF: predefined columns of the key matrices (listed by primes)
  ##  max.sol: maximum number of solutions before exit
  ## RETURN
  ##  a list of design keys, where "design key" means a list with one key matrix per prime
## -------------------------------------------------------------

  Nuqf <- length(PVuqf)
  p <- PVuqf[k]
  r <- NPuqf[k]

  ## PRELIMINARIES: (a) ineligible characters
  ## elimination of the ineligible elements not relevant any more
  m <- matrix(as.logical(b.INELIGuqf[seq.int(k,Nuqf),]), ncol=ncol(b.INELIGuqf))
  keepcol <- apply( m , 2, any)
  if (all(keepcol==FALSE))
    stop("INTERNAL ERROR: all the keepcol are FALSE in the function planor.designkey.recursive")
  b.INELIGtpf <- b.INELIGtpf[ , keepcol,drop=FALSE]
  b.INELIGuqf <- b.INELIGuqf[ , keepcol,drop=FALSE]

  ## selection of the ineligible elements whose 'last' Sylow subgroup is p
  m <- matrix(as.logical(b.INELIGuqf[seq_len(Nuqf)>k,]), ncol=ncol(b.INELIGuqf))
  select <- as.logical(b.INELIGuqf[k,]) & (!apply( m , 2, any))
    b.ineligible <- as.matrix(b.INELIGtpf[ NIVtpf==p, select ,drop=FALSE])

  ## turn the ineligible factorial terms into ineligible p-group characters
  b.ineligible <- representative.basep(b.ineligible,p)
  s <- nrow(b.ineligible)
  ## Identification of the last non-zero coefficients in each ineligible trt character
    ineligible.lnz <- apply(b.ineligible, 2,
                             function(x){max(seq_along(x)[x!=0])}
                           )


  ## PRELIMINARIES: (b) predefined part of PhiStar:
    b.PhiStar <- PREDEF[[k]]

    f <- ncol(b.PhiStar)
  if(s < f)
    stop("Fewer factors than predefined factors")

  ## INITIALISATIONS FOR THE RECURSIVE SEARCH
  ## Calculation of the set of initially admissible elements of U*
  ## admissible = at first all non-zero elements of U*, then a subset of those
  ## admissible.keep = indicator of the U* elements to keep
    b.admissible <- t(convertinto.basep(seq_len((p^r)-1),p))

  nb.admissible <- ncol(b.admissible)
  ## Backtrack search - preliminaries
  iaA <- list(length=s-f)                 ## indices des admissibles pour j
  liaA <- rep(NA,s-f)                     ## nbs d'admissibles pour j
  niaA <- rep(0,s-f)                      ## positions en cours dans les indices

  ## RECURSIVE BACKTRACK SEARCH
  ## 
  PhiStar.solutions <- NULL
  ## 
  jprev <- 0  ;  j <- 1
  while(j > 0){
      b.PhiStar <- b.PhiStar[,seq(f+j-1), drop=FALSE]

    if(jprev < j){
      ## when j increases, the admissible candidates in U*, for
      ## column f+j of the key matrix, must be identified
      ## attention, we assume that the ineligible elements are all in canonical form,
      ## that is, their last non-zero coordinate is equal to 1 - so we do not need it
      ## hence the selection of their first (f+j-1) coordinates at step j below
      ## see Kobilinsky, 2000, end of page 16
         b.ineligible.j <-  as.matrix(b.ineligible[ seq_len(f+j-1), ineligible.lnz==(f+j), drop=FALSE ])


      admissible.keep <- planor.kernelcheck.basep(b.PhiStar, b.admissible, b.ineligible.j, p)
      iaA[[j]] <- seq_len(nb.admissible)[admissible.keep]
      liaA[j] <- length(iaA[[j]])
      niaA[j] <- 0
    }
    if(niaA[j] < liaA[j]){
      niaA[j] <- niaA[j]+1
      newcolj <- (iaA[[j]])[niaA[j]]
        b.PhiStar <- cbind(b.PhiStar,b.admissible[,newcolj])


      ## 
      if(j == (s-f)){ ## WE HOLD A SOLUTION FOR THE CURRENT PRIME p
        ## 
        cat(paste("Solution for p =",p,"\n"))
        ## FIRST CASE: k is maximum, then we look for other solutions or we return
        if(k == Nuqf){
          PhiStar.solutions <- c(PhiStar.solutions,
                                 list(list(matrix(b.PhiStar[,], ncol=ncol(b.PhiStar)))))
          ## start getting out from the recursive search if the number of solutions is ok
          nb.sol <- nb.sol + 1
          if(nb.sol == max.sol){
            ## print success message
            cat("The search is closed: max.sol = ",
                length(PhiStar.solutions), "solutions found \n")
            ## and return
            return(PhiStar.solutions) }
        }
        ## SECOND CASE: k is not maximum, then we have to go to the next prime
        if(k < Nuqf){
          ## SECOND CASE: (a) preparation before the recursive call
          ## identification of the ineligible terms belonging to the current
          ## p-Sylow group and to at least another one not yet accounted for
          m <- matrix(as.logical(b.INELIGuqf[seq_len(Nuqf)>k,]), ncol=ncol(b.INELIGuqf))
          selcol <- as.logical(b.INELIGuqf[k,]) & apply(m , 2, any)
          testcol <- seq_len(ncol(b.INELIGuqf))[ selcol ]
          elimcol <- rep(NA,length(testcol))
          ## identification of the ineligible terms that, in addition,
          ## are not in the kernel of the current PhiStar
          for(ip in seq_along(testcol)){
            ##ineligible.ip <-
            ##   representative.basep(INELIGtpf[NIVtpf==p,testcol[ip],
            ##                                            drop=FALSE],p)
            b.ineligible.ip <-
              representative.basep(as.matrix(b.INELIGtpf[NIVtpf==p,testcol[ip],drop=FALSE]),p)
            b.PhiKer.ip <- (b.PhiStar %*% b.ineligible.ip)%%p
            elimcol[ip] <- all( apply(b.PhiKer.ip, 2, sum) != 0 )
          }
          elimcol <- testcol[elimcol]
          ## SECOND CASE: (b) recursive call
          next.primes.solutions <- Recall(k = k+1,
                                          nb.sol = nb.sol,
                                          PVuqf = PVuqf,
                                          NPuqf = NPuqf,
                                          b.INELIGtpf = b.INELIGtpf[,-elimcol],
                                          NIVtpf = NIVtpf,
                                          b.INELIGuqf = b.INELIGuqf[,-elimcol],
                                          PREDEF = PREDEF,
                                          max.sol=max.sol)
          ## SECOND CASE: (c) solution management if the recursive call has been successful
          if(!is.null(next.primes.solutions)){## we hold at least one new global solution
            nb.sol <- nb.sol + length(next.primes.solutions)
            ## Complete and save the solution
            for(ll in seq(length(next.primes.solutions))){
              next.primes.solutions[[ll]] <- c(list(matrix(b.PhiStar[,], ncol=ncol(b.PhiStar))),
                                             next.primes.solutions[[ll]])
            }
            PhiStar.solutions <- c(PhiStar.solutions, next.primes.solutions)
            if( FREE ){ rm(next.primes.solutions); gc() }
            ## Return if the number of solutions is ok
            if( nb.sol == max.sol ){ return(PhiStar.solutions) }
          } ## end of SECOND CASE: (c)
        } ## end of SECOND CASE (k < Nuqf)
        ## if not returned in the FIRST or SECOND CASES, then look for other solutions
        jprev <- j ; j <- j
      } ## end of "WE HOLD A SOLUTION FOR THE CURRENT PRIME p" (j == s-f)
      else{
        jprev <- j ; j <- j+1
      }
    } ## end of "if(niaA[j] < liaA[j])"
    ## go one column backward
    else{
      jprev <- j ; j <- j-1
    }
  } ## end of "while( j > 0 )"

  ## if we get here, then the number of solutions is not ok yet
  if(k == 1){
    if(is.null(PhiStar.solutions)) cat("No solution found\n")
    else if( nb.sol < max.sol) cat( paste("The search is closed: ",
                                          length(PhiStar.solutions),
                                          " solutions have been found\n") )
  }
  return(PhiStar.solutions)
} ## end planor.designkey.recursive
##---------------------------------------------------------------------------
## 3. SECONDARY CALCULUS FUNCTIONS
##---------------------------------------------------------------------------
planor.kernelcheck.basep <- function(b.PhiStar, b.admissible, b.IneligibleSet, p){
    ## Checks whether any of the N elements in an ineligible set belongs to
    ## the kernel of cbind(PhiStar, admissible[,k]), with PhiStar a p-morphism
    ## matrix and 'admissible' a set of n candidates for making the next
    ## column of PhiStar.
  ## Internal function
    ## ARGUMENTS
    ##  - b.PhiStar: a (r x s) matrix of integers modulo p
    ##  - b.admissible: a (r x n) matrix of integers modulo p
    ##  - b.IneligibleSet: a (s x N) matrix of integers modulo p
    ##  - p: a prime number
    ## RETURN
    ##  a logical vector of length n
    ## DETAILS
    ##  for each column A_k of the matrix 'b.admissible', the function evaluates
    ##  the (r x (s+1)) matrix (PhiStar|A_k). It returns TRUE if no vector
    ##  (I_j'|1)' belongs to the kernel of (PhiStar|A_k), where I_j is the jth column
    ##  of 'IneligibleSet'
## --------------------------------------------------------
    r <- nrow(b.admissible)
    nb.admissible <- ncol(b.admissible)
    nb.ineligible <- ncol(b.IneligibleSet)
    
    
    
    
    
    b.ImagesIS <- (- b.PhiStar %*% b.IneligibleSet)%%p
    test <- rep(2, nb.admissible) ## init par n'importe quelle valeur non NA
    test <- as.logical(.Call("PLANORloopkernelcheck",
                             as.integer(r),
                             as.integer(nb.admissible),
                             as.integer(nb.ineligible),
                             b.ImagesIS,
                             b.admissible,
                             test=as.integer(test)))


    return(test)
} ## end planor.kernelcheck.basep
##---------------------------------------------------------------------------

weightorder.basep <- function(b.mat,p,factnum,blocklog){
  ## Reordering of matrix columns taking account of their weights and trt/block type
  ## Internal function
  ## ARGUMENTS
  ##  - b.mat:  matrix of pseudofactorial effects
  ##  - p: a prime
  ##  - factnum: a numeric factor to identify rows of mat associated
  ##           with the same factor
  ##  - blocklog: a logical vector to identify rows of mat associated
  ##           with a block factor
  ##
  ## RETURN
  ##  the same matrix with columns reordered, and weights given as attributes.
  ##  For each factorial term T confounded with the mean (columns of b.mat),
  ##  the attributes are:
  ##   $trt.weight : number of distinct treatment factors in T
  ##   $trt.pseudoweight : number of distinct treatment pseudofactors in T
  ##   $blc.weight : number of distinct block factors in T
  ##   $blc.pseudoweight : number of distinct block pseudofactors in T
  ## -----------------------------------------------------------------
  nr <- nrow(b.mat)
  nc <- ncol(b.mat)
  labels <- rownames(b.mat)

  ## 1st PART = TRT WEIGHTS
  retour1 <- list( weight= as.double(rep(0, nc)),
                      pseudoweight= as.double(rep(0, nc)),
                      binrank = as.double(rep(0, nc)),
                      modrank = as.double(rep(0, nc)))
  if(sum(!blocklog) > 0){
    
    
    b.mat1 <-  as.matrix(b.mat[!blocklog, , drop=FALSE])

    factnum1 <- factnum[!blocklog]
    labels1 <- labels[!blocklog]
    nr1 <- sum(!blocklog)
     if ( sum(!blocklog) == 1)
        {
      retour1 <- list( weight= (b.mat1[,] > 0)*1,
                      pseudoweight= (b.mat1[,] > 0)*1,
                      binrank = rep(0,nc),
                      modrank = rep(0,nc) )
    }
    else{
retour1 <- .Call("PLANORweightorder",
                       as.integer(nr1), as.integer(nc),
                       as.integer(p),
                       as.integer(factnum1),
                       b.mat1, retour1)
    }
  }

  ## 2nd PART = BLOCK WEIGHTS
  retour2 <- list( weight= as.double(rep(0, nc)),
                      pseudoweight= as.double(rep(0, nc)),
                      binrank = as.double(rep(0, nc)),
                      modrank = as.double(rep(0, nc)))
  if(sum(blocklog) > 0){
       b.mat2 <- matrix(b.mat[blocklog,], ncol=ncol(b.mat))

    factnum2 <- factnum[blocklog]
    labels2 <- labels[blocklog]
    nr2 <- sum(blocklog)
    if(sum(blocklog) == 1){
      retour2 <- list( weight= (b.mat2[,] > 0)*1,
                      pseudoweight= (b.mat2[,] > 0)*1,
                      binrank = rep(0,nc),
                      modrank = rep(0,nc) )
    }
    else{
         retour2 <- .Call("PLANORweightorder",
                       as.integer(nr2), as.integer(nc),
                       as.integer(p),
                       as.integer(factnum2),
                       b.mat2, retour2)
    }
  }
  trt.only <- retour2$weight == 0
  reorder <- order(trt.only,
                   retour1$weight, retour2$weight,
                   retour1$binrank, retour1$modrank)

  
  

  b.mat <- b.mat[,reorder,drop=FALSE]


  attributes(b.mat)$trt.weight <- retour1$weight[reorder]
  attributes(b.mat)$trt.pseudoweight <- retour1$pseudoweight[reorder]
  attributes(b.mat)$blc.weight <- retour2$weight[reorder]
  attributes(b.mat)$blc.pseudoweight <- retour2$pseudoweight[reorder]
  rownames(b.mat) <- labels

  return(b.mat)
} ## end weightorder.basep
##---------------------------------------------------------------------------
## Preliminary printing conventions

printpower <- function(p.k){
  cat(paste("\n********** Prime ", p.k), " design **********\n\n")
}

wprofile <- function(x){
  ## calculates the "profile" of the values in a vector of integers
  w <- unique(x)
  nb <- apply( outer(w, x, "==")*1, 1, sum )
  profile <- paste(w,nb,sep="^")
  return(profile)
}
##---------------------------------------------------------------------------
## 4. OUTPUT FUNCTIONS
##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
## "planor.design.levels"
## Generate a full factorial n1 x n2 x ... x ns design with columns
## considered as factors
## ARGUMENTS
## - key a vector of integers of length s
## - start  an integer from where to start the series of symbols
## RETURN an integer matrix with prod(n) rows and s columns giving all
##  combinations along the rows, in lexicographic order
## Examples
## planor.design.levels(rep(2,3))
## --------------------------------------

planor.design.levels <- function(key,start=1){

    OUT <- crossing(key,start=start)
    for(j in ncol(OUT)){
        OUT[[j]] <- factor(OUT[[j]])
    }
     return(OUT)
}

## --------------------------------------
## "planor.design" method for "numeric"
## --------------------------------------
setMethod("planor.design", signature(key="numeric"),
          definition=planor.design.levels)


##--------------------------------------------------------------------------
## printgmat : utility function
## Display a matrix.
## To bound the display amount,  we print the first maxprint first
## lines and columns, only. 
##--------------------------------------------------------------------------
printgmat <- function(mat, maxprint=getOption("planor.max.print", default=20)) {
  prnrow <- min(nrow(mat), maxprint)
  prncol <- min(ncol(mat), maxprint)
  prmat <- matrix(c(mat[1:prnrow, 1:prncol]),
                  nrow=prnrow,
                  dimnames=list(dimnames(mat)[[1]][1:prnrow],
                    dimnames(mat)[[2]][1:prncol]))
  print( prmat )
  if ( prnrow < nrow(mat))
    cat(paste("\nThe first",  prnrow, "rows on a total of", nrow(mat)))
  if ( prncol < ncol(mat))
    cat(paste("\nThe first",  prncol, "columns on a total of", ncol(mat),"\n"))
  cat("\n")
} ## end printgmat
##n ------------------------------------
