#
# These routines are designed to automate the process of creating the rather complicated
# landscape object.
#
#
#
#

#
# initialize the basic landscape object
#
landscape.new.empty <- function()
{
#    list(list(RndChooseProb=NULL,StartGen=NULL,Extinct=NULL,Carry=NULL,Localprob=NULL,S=NULL,R=NULL,M=NULL))
# KKM 5.20.05................................................................
  #old version
  #tmpdemo <- list(localdem=NULL,epochs=NULL)
  #new version
  tmpdemo <- list(localdem=NULL,localdemK=NULL,epochs=NULL)
#............................................................................
#    list(localdem=list(list(LocalS=NULL,LocalR=NULL,LocalM=NULL)),epochs=NULL)
  rland <- list(intparam=NULL,switchparam=NULL,floatparam=NULL,demography=tmpdemo,loci=NULL,individuals=NULL)
  rland
}

landscape.new.default <- function()
{
  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland)
  rland <- landscape.new.floatparam(rland)
  rland <- landscape.new.switchparam(rland)
  rland
}

landscape.new.example <- function()
{
  rland <- NULL
  rland <- landscape.new.empty()
  
  rland <- landscape.new.intparam(rland, h=2, s=2)
  rland <- landscape.new.switchparam(rland,mp=0)
  rland <- landscape.new.floatparam(rland)


  S <- matrix(c(0, 0,
                1, 0), byrow=TRUE, nrow = 2)
  R <- matrix(c(0, 1.1,
                0, 0), byrow=TRUE, nrow = 2)
  M <- matrix(c(0, 0,
                0, 1), byrow=TRUE, nrow = 2)
  rland <- landscape.new.local.demo(rland,S,R,M)
  
  S <- matrix(rep(0,16), nrow = 4)
  R <- matrix(rep(0,16), nrow = 4)
  M <- matrix(rep(0,16), nrow = 4)
  
  rland <- landscape.new.epoch(rland,S=S,R=R,M=M,carry=c(1000,1000))
  
  rland <- landscape.new.locus(rland,type=0,ploidy=2,mutationrate=0.001,transmission=0,numalleles=5)
  rland <- landscape.new.locus(rland,type=1,ploidy=1,mutationrate=0.005,numalleles=3,frequencies=c(.2,.2,.6))
  rland <- landscape.new.locus(rland,type=2,ploidy=2,mutationrate=0.007,transmission=0,numalleles=6,allelesize=75)
  
  rland <- landscape.new.individuals(rland,c(50,0,50,0))
  rland
}

#
# these routines set up the list of intparams with some sort of reasonable defaults
# the first within a landscape, the second independantly

landscape.new.intparam <- function(rland,h=1,s=1,cg=0,ce=0,totgen=1000,maxland=200000)
{
  rland$intparam <- list(h,s,0,0,cg,ce,totgen,0,maxland)
  names(rland$intparam) <- c("habitats","stages","locusnum","numepochs","currentgen","currentepoch","totalgens","numdemos","maxlandsize")
  rland
}
  
#new.intparam <- function(h=1,s=2,l=1,ne=1,cg=0,ce=0,totgen=1,nd=1,maxland=200000)
#{
#  rl <- list(h,s,l,ne,cg,ce,totgen,nd,maxland)
#  names(rl) <- c("habitats","stages","locusnum","numepochs","currentgen","currentepoch","totalgens","numdemos","maxlandsize")
#  rl
#}


#
# these routines set up the list of floatparams with some sort of reasonable defaults
# the first within a landscape, the second independantly

landscape.new.floatparam <- function(rland, s=0)
{
  rland$floatparam <- list(s)
  names(rland$floatparam) <- c("selfing")
  rland
}

#new.floatparam <- function(s=0)
#{
#  rl <- list(s)
#  names(rl) <- c("selfing")
#  rl
#}


#
# these routines set up the list of switchparams with some sort of reasonable defaults
# the first within a landscape, the second independantly

# KKM 5.20.05 new versions..................................................
landscape.new.switchparam <- function(rland, re=0,rd=0,mp=1,dd=0)
{
  rland$switchparam <- list(re,rd,mp,dd)
  names(rland$switchparam) <- c( "randepoch","randdemo","multp","densdepdemo")
  rland
}

#new.switchparam <- function(re=1,rd=1,mp=1,dd=0)
#{
#  rl <- list(re,rd,mp,dd)
#  names(rl) <- c( "randepoch","randdemo","multp","densdepdemo")
#  rl
#}

#
# is.nsquare
#
# is a square matrix with size n x n

is.nsquare <- function(M,n)
{
  ((dim(M)[1] == dim(M)[2]) && (dim(M)[1] == n))
}

#
# landscape.new.local.demo
#
# Initializes local demographies.  This function is going to require the number of stages 
# in each demography (from "intparam").  This number could also be calculated from user input.
# Will also require three actual matrices (S,R,M for survival, reproduction,
# and male function, respecively) input by user.  Each matrix also has to be the same size.

landscape.new.local.demo <- function(rland,S,R,M,k=0)
{
 if (k==0){  #matrix at ZPD or there is no density dependence
  if (is.null(rland$demography$localdem))
    {
      rland$demography$localdem <- list(NULL)
      demonum <- 1
    }
  else
    {
      demonum <- length(rland$demography$localdem) + 1
      rland$demography$localdem[[demonum]] <- list(LocalS=NULL,LocalR=NULL,LocalM=NULL)
    }

  if (is.nsquare(S,rland$intparam$stages) &&
      is.nsquare(R,rland$intparam$stages) &&
      is.nsquare(M,rland$intparam$stages))
    {
      rland$demography$localdem[[demonum]]$LocalS <- S
      rland$demography$localdem[[demonum]]$LocalR <- R
      rland$demography$localdem[[demonum]]$LocalM <- M
      rland$intparam$numdemos <- length(rland$demography$localdem)
    }
  else
    {
      stop("Matricies do not conform to stages set in intparam!")
    }
} #end if k==0
 else{ #k==1; matrix at carrying capacity
   if (is.null(rland$demography$localdemK))
    {
      rland$demography$localdemK <- list(NULL)
      demonumK <- 1
    }
  else
    {
      demonumK <- length(rland$demography$localdemK) + 1
      rland$demography$localdemK[[demonumK]] <- list(LocalS=NULL,LocalR=NULL,LocalM=NULL)
#      if (demonumK > rland$intparam$numdemos) NEED AN ERROR TRAP HERE
    }

  if (is.nsquare(S,rland$intparam$stages) &&
      is.nsquare(R,rland$intparam$stages) &&
      is.nsquare(M,rland$intparam$stages))
    {
      rland$demography$localdemK[[demonumK]]$LocalS <- S
      rland$demography$localdemK[[demonumK]]$LocalR <- R
      rland$demography$localdemK[[demonumK]]$LocalM <- M
    }
  else
    {
      stop("Matrices do not conform to stages set in intparam!")
    }
  } #end else
  rland
}

# 
# landscape.new.epoch
#
# initializes an epoch.  This includes creating landscape matrices that describe survival,
# reproduction and male function. These matrices are square and have numbers of cols and rows
# equal to the number of pops times number of demographic stages within pops
# (intparam$habitats * intparam$stages)
#
# It also needs to specify the extinction rates and carrying
# capacities of each population (intparam$habitats)
#
# finally, the function needs to define the probability of selecting particular local
# demographies for each population.  (this occurs if switchparam$randdemo==1)
#

landscape.new.epoch <- function(rland,S=NULL,R=NULL,M=NULL,epochprob=1,startgen=0,extinct=NULL,carry=NULL,localprob=NULL)
{
  if (is.null(rland$demography$epochs))
    {
      rland$demography$epochs <- list(NULL)
      epochnum <- 1
    }
  else
    {
      epochnum <- length(rland$demography$epochs) + 1
      rland$demography$epochs[[epochnum]] <- list(RndChooseProb=NULL,StartGen=NULL,Extinct=NULL,
                                                  Carry=NULL,Localprob=NULL,S=NULL,R=NULL,M=NULL)
    }
  
  rland$demography$epochs[[epochnum]]$RndChooseProb <- epochprob
  rland$demography$epochs[[epochnum]]$StartGen <- startgen
  if (is.null(extinct))
    {
      extinct <- rep(0, rland$intparam$habitats)
    }  
  if (length(extinct) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$Extinct <- extinct
    }
  else
    {
      stop("Wrong size for extinct vector", dim((extinct)))
    }
  
  if (is.null(carry)) 
    {
      carry <- rep(1000, rland$intparam$habitats)
    }
  if (length(carry) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$Carry <- carry
    }
  else
    {
      stop("Wrong size for carry vector")
    }

  numdem <- length(rland$demography$localdem)
  if (is.null(localprob)) 
    {
      localprob <- rep(1/numdem, numdem)
    }
  if (length(localprob == numdem))
    {
      rland$demography$epochs[[epochnum]]$Localprob <- localprob
    }
  else
    {
      stop("Wrong size for localprob vector")
    }

  matrixsize <- rland$intparam$habitats * rland$intparam$stages
  
  if (is.null(S))
    {
      S <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.null(R))
    {
      R <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.null(M))
    {
      M <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.nsquare(S,matrixsize) && 
      is.nsquare(R,matrixsize) && 
      is.nsquare(M,matrixsize))
    {
      rland$demography$epochs[[epochnum]]$S <- S
      rland$demography$epochs[[epochnum]]$R <- R
      rland$demography$epochs[[epochnum]]$M <- M
      rland$intparam$numepochs <- length(rland$demography$epochs)
    }
  else
    {
      stop("S, R, and M matricies are not the correct size")
    }
  rland
}



#
# landscape.new.locus
#
# This function should create a locus that is populated with alleles.
# it is passed the: type, ploidy, mutation rate, transmission, number of alleles, allele size
# (needed only for alleles of type 2) and vector representing the frequency of each allele
# (if empty, allelic distribution is uniform)
#
#
# the function should then create a list of alleles: assign 0 to allele birth date,
# assign proportion based upon the frequency vector, assign state depending upon allele type
#
#

#
# new locus that takes into account allele states (needed to import coalescence sims)
#

landscape.new.locus <- function (rland, type = 0, ploidy = 1, mutationrate = 0, transmission = 1, 
    numalleles = 2, allelesize = 50, frequencies = NULL, states = NULL) 
{
    if (!(is.list(rland$loci))) {
        rland$loci <- list(list(type = 0, ploidy = 0, trans = 0, 
            rate = 0, alleles = 0))
        locusnum <- 1
    }
    else {
        locusnum <- length(rland$loci) + 1
        rland$loci[[locusnum]] <- list(type = 0, ploidy = 0, 
            rate = 0, trans = 0, alleles = 0)
    }
    rland$intparam$locusnum <- locusnum
    if (type >= 0 && type <= 2) {
        rland$loci[[locusnum]]$type <- as.integer(typelookup(type))
    }
    else {
        stop("Invalid type of locus")
    }
    if (ploidy == 1 || ploidy == 2) {
        rland$loci[[locusnum]]$ploidy <- as.integer(ploidy)
    }
    else {
        stop("Invalid ploidy count")
    }
    rland$loci[[locusnum]]$rate <- mutationrate
    if (transmission == 0 || transmission == 1) {
        rland$loci[[locusnum]]$trans <- as.integer(transmission)
    }
    else {
        stop("Invalid transmission number")
    }
    if (numalleles >= 0) {
        rland$loci[[locusnum]]$alleles <- makealleles(type, numalleles, 
            allelesize, frequencies, states)
    }
    else {
        stop("Need non-negative numbers of alleles")
    }
    rland
}

typelookup <- function(type)
  {
    type <- type + 251
    type
  }



makealleles <- function(type,numalleles,allelesize,frequencies,states)
{
  retval <- 0

  if(is.null(frequencies))
    {
      frequencies <- rep(1.0/numalleles, numalleles)
    }

  if(length(frequencies) != numalleles)
    {
      stop("Frequency list is not the right size")
    }
  
  if(type == 0 || type == 1)
    {
      retval <- vector("list", numalleles)
      for (x in 1:numalleles)
        {
          retval[[x]]$aindex <- as.integer(x)
          retval[[x]]$birth <- as.integer(0)
          retval[[x]]$prop <- frequencies[x]
          if (is.null(states))
            {
              retval[[x]]$state <- as.integer(x)
            } else
          {
            retval[[x]]$state <- as.integer(states[x])
          }
        }
    }
  else if(type == 2)
    {
      retval <- vector("list", numalleles)
      for (x in 1:numalleles)
        {
          retval[[x]]$aindex <- as.integer(x)
          retval[[x]]$birth <- as.integer(0)
          retval[[x]]$prop <- frequencies[x]
          if (is.null(states))
            {
              retval[[x]]$state <- geneseq(allelesize)
            } else {
              retval[[x]]$state <- states[x]
            }
        }
    }
  retval          
}

geneseq <- function(size)
  {
    if(size <= 0)
      {
        stop("Allelesize must be positive")
      }
    retval <- floor(runif(size)*4)
    retval[retval==0] <- 'A'
    retval[retval==1] <- 'T'
    retval[retval==2] <- 'G'
    retval[retval==3] <- 'C'
    paste(retval, sep="", collapse="")
  }
    

#
# landscape.new.individuals
#
# should take the landscape as it stands and use the c++ method Landscape::popsizeset to populate
# the landscape with individuals, this function does expect a distribution of individuals for
# each habitat*stage  combination in the landscape
#

landscape.new.individuals <- function(rland, PopulationSizes)
  {
    rland <- landscape.coerce(rland,noind=T)
    rland <- .Call("populate_Rland",rland,PopulationSizes,PACKAGE="rmetasim")
    rland
  }

#
# a convenience function that provides the highest column number for demographic information in
# the individuals matrix
#
landscape.democol <- function()
  {
    as.integer(.Call("num_demo_cols",PACKAGE="rmetasim"))
  }
