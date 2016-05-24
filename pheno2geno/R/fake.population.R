#
# fake.population.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: fake.population, fakePheno.internal, fakeFounders.internal, simBC.internal
#           fakePhysicalMap.internal, convertMap.internal, fakeMixUps
#

# fake.population
#
# DESCRIPTION:
#  Simulating object of class population 
# OUTPUT:
#  An object of class population
#
fake.population <- function(n.founders = 4, n.offspring = 100, n.markers=100,n.chromosomes=10, type = c("riself", "f2", "bc", "risib"), n.mixups=0, verbose=FALSE,...){
  type <- match.arg(type)
  if(!(is.numeric(n.founders))) stop("n.founders should be numeric\n")
  if(!(is.numeric(n.offspring))) stop("n.offspring should be numeric\n")
  if(!(is.numeric(n.markers))) stop("n.markers should be numeric\n")
  if(!(is.numeric(n.chromosomes))) stop("n.chromosomes should be numeric\n")
  if(n.founders<4){
    warning("too small n.founders, changing to 4\n")
    n.founders <- 4
  }
  
  if(!(n.founders%%2==0)){
    warning("n.founders should be even, changing to",n.founders+1,"\n")
    n.founders <- n.founders+1
   }
  if(length(type)>1) type <- type[1]
  if(n.offspring<10){
    warning("too small n.offspring, changing to 10\n")
    n.offspring <- 10
  }
  if(n.markers<10){
    warning("too small n.markers, changing to 10\n")
    n.markers <- 10
  }
  if(n.markers<n.chromosomes){
    warning("n.markers cannot be smaller than n.chromosomes, changing n.markers to 10*n.chromosomes\n")
    n.markers <- 10*n.chromosomes
  }

  ### fuction itself
  map <- sim.map(rep(100,n.chromosomes),n.mar=(n.markers/n.chromosomes), include.x=FALSE,)
  fake <- sim.cross(map,type=type, n.ind=n.offspring, ...)
  geno <- t(pull.geno(fake))
  if(type=="bc"){
    geno <- t(apply(geno,1,simBC.internal))
  }
  map <- convertMap.internal(map)
  physicalMap <- fakePhysicalMap.internal(map)
  colnames(geno) <- paste("RIL",1:ncol(geno),sep="_")
  pheno <- t(apply(geno,1,fakePheno.internal))
  rownames(pheno) <- rownames(geno)
  colnames(pheno) <- colnames(geno)
  founders <- t(apply(pheno,1,fakeFounders.internal,n.founders))
  rownames(founders) <- rownames(geno)
  colnames(founders) <- 1:n.founders
  colnames(founders)[1:(n.founders/2)] <- paste("Founder",1,1:(n.founders/2),sep="_")
  colnames(founders)[(n.founders/2+1):n.founders] <- paste("Founder",2,(n.founders/2+1):n.founders,sep="_")
  #geno[which(geno==2)] <- 0
  foundersGroups <- c(rep(0,(n.founders/2)),rep(1,(n.founders/2)))
  if(n.mixups>0){
    pheno <- fakeMixUps.internal(pheno,n.mixups)
  }
  population <- create.population(pheno, founders, foundersGroups, geno, map, physicalMap, verbose=verbose)
  class(population)[2]<- type
  check.population(population)
  invisible(population)
}

############################################################################################################
#                                  *** fakePheno.internal ***
#
# DESCRIPTION:
#  subfunction of fake.population - simulating phenotype data using genotype data simulated by sim.cross
# OUTPUT:
#  row of offspring phenotype matrix
############################################################################################################
fakePheno.internal <- function(genoRow,maxScale=10,maxError=3){
  scalingF <- runif(1,1,maxScale)
  errorF <- runif(length(genoRow),0,maxError)
  genoRow <- (genoRow*scalingF) + errorF
  invisible(genoRow)
}

############################################################################################################
#                                  *** fakeFounders.internal ***
#
# DESCRIPTION:
#  subfunction of fake.population - simulating founders phenotype data using offspring phenotype data
# OUTPUT:
#  row of parental phenotype matrix
############################################################################################################
fakeFounders.internal <- function(phenoRow,n.founders){
  errorF <- runif(n.founders,0,2)
  up <- runif(1,-1,1)
  diffExprRate <- runif(1,0.1,1)
  cur_mean <- mean(phenoRow)
  if(up>=0){
    foundersRow <- c(rep((cur_mean-diffExprRate*cur_mean),(n.founders/2)),rep(cur_mean+diffExprRate*cur_mean,(n.founders/2))) + errorF
  }else{
    foundersRow <- c(rep((cur_mean+diffExprRate*cur_mean),(n.founders/2)),rep(cur_mean-diffExprRate*cur_mean,(n.founders/2))) - errorF
  }
  invisible(foundersRow)
}

############################################################################################################
#                                  *** fakePhysicalMap.internal ***
#
# DESCRIPTION:
#  simulating physical map using genetic one
# OUTPUT:
#  map of the same type
############################################################################################################
fakePhysicalMap.internal <- function(map){
  for(i in 1:nrow(map)){
    errorF <- runif(1,0,100)
    if(errorF>90){
      newChrom <- round(runif(1,1,(max(map[,1]))))
      map[i,1] <- newChrom
    }
  }
  map <- cbind(map,map[,2]+300)
  invisible(map)
}

############################################################################################################
#                                  *** convertMap.internal ***
#
# DESCRIPTION:
#  convert rqtl type map into population type one
# OUTPUT:
#  map of population class type - rownames - names of the markers, first column - numbers of chromosomes,
#    second - position of the marker on chromosome
############################################################################################################
convertMap.internal <- function(map){
  map_ <- NULL
  for(i in 1:length(map)){
    cur_chr <- cbind(rep(i,length(map[[i]])),map[[i]])
    map_ <- rbind(map_,cur_chr)
  }
  invisible(map_)
}

############################################################################################################
#                                  *** simBC.internal ***
#
# DESCRIPTION:
#  convert rqtl type map into population type one
# OUTPUT:
#  map of population class type - rownames - names of the markers, first column - numbers of chromosomes,
#    second - position of the marker on chromosome
############################################################################################################
simBC.internal <- function(genoRow){
  n.ind <- length(genoRow)
  toCount <- round(runif(1,1,2))
  toChange <- (3-toCount)
  numberOfToCount <- sum(genoRow==toCount)
  errorRate <- round(n.ind*(runif(1,0,5)/100))
  expected <- round(0.75*n.ind)
  cat(numberOfToCount,expected,errorRate,"\n")
  while(numberOfToCount<(expected-errorRate)){
    genoRow[round(runif(1,1,length(genoRow)))] <- toCount
    numberOfToCount <- sum(genoRow==toCount)
  }
  while(numberOfToCount>(expected+errorRate)){
    genoRow[round(runif(1,1,length(genoRow)))] <- toChange
    numberOfToCount <- sum(genoRow==toCount)
  }
  other <- sum(genoRow==toChange)
  cat(numberOfToCount,other ,"\n")
  invisible(genoRow)
}

############################################################################################################
#                                  *** mixUpPheno.internal ***
#
# DESCRIPTION:
#  convert rqtl type map into population type one
# OUTPUT:
#  map of population class type - rownames - names of the markers, first column - numbers of chromosomes,
#    second - position of the marker on chromosome
############################################################################################################
fakeMixUps.internal<- function(pheno, n.mixups){
  for(i in 1:n.mixups){
    toMix <- sample(1:nrow(pheno),2)
    print(toMix)
    temp_ <- pheno[toMix[1],]
    pheno[toMix[1],] <- pheno[toMix[2],]
    pheno[toMix[2],] <- temp_
  }
  invisible(pheno)
}
