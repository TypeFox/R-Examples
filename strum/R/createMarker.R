#==============================================================================
# File: createMarker.R
#
# Author: Nathan Morris
#
# Notes: Function to create strumMarker object using HapMap data.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Create strumMarker object given hapmap data
#------------------------------------------------------------------------------
createStrumMarker = function(
                      hapMapData,
                      populationRecombRate = 50,
                      errorRate = 0,
                      mutationRate = 0,
                      missingRate = 0,
                      coding = c(0,1,2),
                      returnIBD = FALSE,
                      intervalIBD = 10
                    )
{
  rownames(hapMapData) = hapMapData$rsID

  nonMarkerCols = which(names(hapMapData) %in% c("rsID","phys_position","chr"))

  # calculate allele frequency and determine the minor Allele for each marker
  #---------------------------------------------------------------------------
  minorMajorAllele = apply(hapMapData[,-nonMarkerCols],1,
                           function(alleles)
                           {
                             ret = table(alleles)
                             ret = ret/sum(ret)
                             return(names(sort(ret)))
                           })

  # recode the genotypes; replace minors with 1 and majors with 0
  #---------------------------------------------------------------
  hapD = as.matrix(hapMapData[,-nonMarkerCols])
  hapD = as.data.frame(t(hapD))

  hapM = mapply(function(D, M)
                {
                  if( length(M) == 1 )
                    ret = rep(0, length(D))
                  else
                    ret = ifelse(D==M[1], 1, 0)
                  return(ret)
                }, hapD, minorMajorAllele)

  # Create object
  #---------------
  retObj <- new("strumMarker",
                markerFacts = data.frame(markerName = hapMapData$rsID,
                                         chrom = hapMapData$chr,
                                         mapPos = hapMapData$phys_position/1000000,
                                         minorAllele = sapply(minorMajorAllele,function(a) return(ifelse(length(a)==1, NA, a[1]))), 
                                         majorAllele = sapply(minorMajorAllele,function(a) return(a[length(a)]))),
                haplotypes = t(hapM),
                populationRecombRate = populationRecombRate,
                errorRate = errorRate,
                mutationRate = mutationRate,
                missingRate = missingRate,
                coding = coding,
                returnIBD = returnIBD,
                intervalIBD = intervalIBD)

  .printInfoLine("Creating strumMarker", "Done", 52, 0)

  return(retObj)
}
