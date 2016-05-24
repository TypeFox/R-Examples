# James Niehaus
# assign.R
# Tutorial
# Spring 2003
# calculate the assignment index for individuals within
# a landscape


##
# Returns a list with 3 elements
# the first is an array stating whether an individual (row) is misassigned
# the second is an array that gives the probabilites of an individual (row) belonging to a
#   popluation (col)
# the third is a list that gives the probabilites of an individuals (list index) loci (row)
#   belonging to a specific population (col)
assignmentTest <- function(Rland, verbose=F) 
{  
  if( is.landscape(Rland) )
  {
    nPops <- Rland$intparam$habitats
    nLoci <- Rland$intparam$locusnum
    nInd  <- countPopulation(0,Rland)

    sampleset <- sampleFromSubPopulations(Rland, 4)
    
    indProb <- array(rep(0,(nPops * nLoci)), dim=c(nInd,nPops))    
    lociProb <- list()
    misassigned <- data.frame(misassigned=rep(0, length(sampleset)),
                              row.names = sampleset )
    lociIndex <- matrix(nrow=nLoci, ncol=nPops, byrow=TRUE)    
    
    # Constuct the Individual and Loci probabilities tables
    for (curInd in sampleset)
    {
      if (verbose)
        {
          print(c("*** curInd:",curInd))
        }
      aIndex <- rep(1, nPops)      
      for (curLoci in 1:nLoci)
        {
          if (verbose) 
            {
              print(c("curLoci:",curLoci))
            }
          freq <- indxfreqNormal(curInd, curLoci, Rland)
          lociIndex[curLoci,] <- lociAIndex(freq, curInd, curLoci, Rland) 
          aIndex <- aIndex * lociIndex[curLoci,]
        }
      lociProb[[curInd]] <- lociIndex
      indProb[curInd,] <- aIndex
    }
    indProb <- -log(indProb)
    indProb[(indProb == Inf)] <- 0

    # construct the missassigned table
    for (curInd in 1:nInd)
      {
        myProb <- indProb[curInd,]
        myPop <- landscape.populations(Rland)[curInd]
        misassigned[[1]][which(sampleset == curInd)] <-
          !(myProb[myPop] == max(myProb))
      }
    
    list(misassigned, indProb, lociProb)
  }
  else
  {
    print("Parameter is not a landscape")
  }
}

sampleFromSubPopulations <- function (Rland, sampleSize)
{
  sample <- c()
  nPops <- Rland$intparam$habitats
  pops <- landscape.populations(Rland)
  for (x in 1:nPops) {
    sample <- c(sample, sample(which(pops == x), sampleSize))
  }
  sort(sample)
}
countPopulation <- function (nPop, Rland)
{
  if(nPop == 0)
  {
    dim(Rland$individuals)[1]    
  }
  else
  {
    sum(landscape.populations(Rland) == nPop)
  }
}

indxfreqWithout <-function(IndNum, lnum=1, Rland)
  {
    lv<-landscape.locus(Rland,lnum)[,c(FALSE,FALSE,FALSE,rep(TRUE,(ncol(landscape.locus(lnum,Rland))-3)))];
    if (landscape.ploidy(Rland)[lnum]==1)
      {
        lv<-lv[-IndNum]
        table(landscape.populations(Rland)[-IndNum],lv)
      }
    else
      {        
        lv<-lv[-IndNum,]
        lv <- c(lv[,1],lv[,2])
        table(rep(landscape.populations(Rland)[-IndNum],2),lv)
      }
  }

indxfreqNormal <- function(IndNum, lnum=1, Rland)
  {
    nPops <- Rland$intparam$habitats
    myPop <- landscape.populations(Rland)[IndNum]

    freq <- indxfreqWithout(IndNum, lnum, Rland)
    for(x in 1:nPops)
      {
        if (x == myPop)
          {
            geneCount <- countPopulation(x,Rland) - 1
          }
        else
          {
            geneCount <- countPopulation(x,Rland)
          }
        geneCount <- geneCount * landscape.ploidy(Rland)[lnum]

	freq[x,] <- freq[x,] / geneCount
      }
    freq[freq == NA] <- 0
    freq
  }

lociAIndex <- function(freq, curInd, curLoci, Rland) 
  { 
    myAlleleIndex <- getMyAlleleIndex(freq, curInd, curLoci, Rland)
    if (length(myAlleleIndex) == 1)
      {
        retval <- freq[,myAlleleIndex]
      }
    else if (length(myAlleleIndex) == 2)
      {
        if (myAlleleIndex[1] == myAlleleIndex[2])
          {
            retval <-(freq[,myAlleleIndex[1]])^2
          }
        else
          {
            retval <- 2 * freq[,myAlleleIndex[1]] * freq[,myAlleleIndex[2]]
          }
      }
    else
      {
        retval <- 1
      }
    retval
  }

getMyAlleleIndex <- function(freq, curInd, curLoci, Rland)
  {
    myAllele <- landscape.locus( Rland,lnum=curLoci)[curInd,c(-1,-2,-3)]
    if (landscape.ploidy(Rland)[curLoci] == 1)
    {
      rv <- which(as.numeric(colnames(freq))==myAllele)
    }
    else
    {
      rv <- c(which(as.numeric(colnames(freq))==myAllele[1]),which(as.numeric(colnames(freq))==myAllele[2]))      
    }
    rv
  }
 
