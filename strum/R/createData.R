#==============================================================================
# File: createStrumData.R
#
# Author: Nathan Morris
#
# Notes: Make strumData given data & data type.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# A function to create strumData object
#------------------------------------------------------------------------------
createStrumData = function(inData, dType, ibdFileName=NULL, fileType="SAGE")
{
  # Check if a dataframe. If not, stop().
  #---------------------------------------
  if( !is.data.frame(inData) )
      stop("Data must be a data.frame!")

  # If empty data, stop().
  #------------------------
  if( !length(inData) )
    stop("Data must be non-empty data.frame!")

  if( dType == "Pedigree" )
  {
    # Check for the required fields. If not, stop().
    #------------------------------------------------
    if( .checkPed(inData) == 1 )
      stop("If dataType is Pedigree, data must contain: family, id, father, and mother.")

    # Build pedigree structure & get phi matrices
    #---------------------------------------------
    pedData = .orderPedigrees(inData)
    #pos = .getParentPos(pedData)
    peds = split(pedData, pedData$family)
    phi = .getPhi(peds)

    retObj = new("strumData", dataType = dType, dataVals = pedData, phi = phi)

    # Import IBD information
    #------------------------
    if( !is.null(ibdFileName) )
      retObj = importIBD(retObj, ibdFileName, fileType)

  } else if( dType == "RawData" )
  {
    # Check for the required fields. If not, make them.
    #---------------------------------------------------
    if( .checkPed(inData) == 1 )
    {
      N = nrow(inData)
      #family = rep(1, N)
      id     = seq(1:N)
      fm     = rep(0, N)

      if( !("family" %in% names(inData)) )
      {
        orgNames = names(inData)
        inData = cbind(inData, id)
        names(inData) = c(orgNames, "family")
      }
      if( !("id" %in% names(inData)) )
      {
        orgNames = names(inData)
        inData = cbind(inData, id)
        names(inData) = c(orgNames, "id")
      }
      if( !("mother" %in% names(inData)) )
      {
        orgNames = names(inData)
        inData = cbind(inData, fm)
        names(inData) = c(orgNames, "mother")
      }
      if( !("father" %in% names(inData)) )
      {
        orgNames = names(inData)
        inData = cbind(inData, fm)
        names(inData) = c(orgNames, "father")
      }
    }

    retObj = new("strumData", dataType = dType, dataVals = inData)

  } else
  {
    stop("Data type must be Pedigree or RawData!")
  }

  .printInfoLine("Creating strumData", "Done", 52, 0)

  return(retObj)
}

#------------------------------------------------------------------------------
# Get kinship coefficient matrix using pedigree package
#------------------------------------------------------------------------------
.getPhi = function(pos)
{
  phi = lapply(pos,
               function(pedi) 
               {
                 simplePed = pedi[, c("id", "father", "mother")]
                 
                 myn = nrow(simplePed)

                 makeA(simplePed,which = rep(TRUE,myn))

                 A = read.table("A.txt")

                 Amat = matrix(0, nrow=myn, ncol=myn)

                 for( i in 1:nrow(A) )
                 {
                   x = A[i,1]
                   y = A[i,2]
                   v = A[i,3]
                   Amat[x,y] = v
                   Amat[y,x] = v
                 }

                 rownames(Amat) = simplePed$id
                 colnames(Amat) = simplePed$id

                 return(Amat)
               })

  file.remove("A.txt")

  return(phi)
}
