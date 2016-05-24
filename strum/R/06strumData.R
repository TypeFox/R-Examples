#==============================================================================
# File: strumData.R
#
# Author: Nathan Morris
#
# Notes: strumData class & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumData class
#------------------------------------------------------------------------------
setClass("strumData",
         representation(
           dataType = "character", 
           dataVals = "data.frame", 
           phi      = "list",
           ibd      = "strumIBD"),	
         prototype = list(
           dataType = "",
           dataVals = data.frame(),
           phi      = list(),
           ibd      = new("strumIBD"))
)

#------------------------------------------------------------------------------
# 'dataType' accessor functions:
#------------------------------------------------------------------------------
setGeneric('dataType', function(object) standardGeneric('dataType'))
setMethod('dataType', signature(object = 'strumData'),
          function(object)
          {
            return(object@dataType)
          }
)

setGeneric('dataType<-', function(object,value) standardGeneric('dataType<-'))
setMethod('dataType<-', signature(object = 'strumData'),
          function(object, value)
          {
            if( value %in% c("Pedigree", "RawData", "MeanCov") )
            {
              object@dataType <- value
              return(object)
            } else
            {
              stop("Incorrect dataType. Must be 'Pedigree', 'RawData', or 'MeanCov'")
            }
          }
)

#------------------------------------------------------------------------------
# 'dataVals' accessor functions:
#------------------------------------------------------------------------------
setGeneric('dataVals', function(object) standardGeneric('dataVals'))
setMethod('dataVals', signature(object = 'strumData'),
          function(object)
          {
            return(object@dataVals)
          }
)

setGeneric('dataVals<-', function(object,value) standardGeneric('dataVals<-'))
setMethod('dataVals<-', signature(object = 'strumData'),
          function(object, value)
          {
            object@dataVals <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'phi' accessor functions:
#------------------------------------------------------------------------------
setGeneric('phi', function(object) standardGeneric('phi'))
setMethod('phi', signature(object = 'strumData'),
          function(object)
          {
            return(object@phi)
          }
)

setGeneric('phi<-', function(object,value) standardGeneric('phi<-'))
setMethod('phi<-', signature(object = 'strumData'),
          function(object, value)
          {
            object@phi <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'ibd' accessor functions:
#------------------------------------------------------------------------------
setGeneric('ibd', function(object) standardGeneric('ibd'))
setMethod('ibd', signature(object = 'strumData'),
          function(object)
          {
            return(object@ibd)
          }
)

#------------------------------------------------------------------------------
# "[" generic functions
#------------------------------------------------------------------------------
setMethod("[", "strumData",
          function(x, i, j, ... )
          {
            if( x@dataType == "RawData" )
            {
              x@dataVals = x@dataVals[i,j,...]

            } else if( x@dataType == "Pedigree" )
            {
              #check for unordered pedigrees
              tmp = x@dataVals[i,] #Get ith row

              #is it maybe pedigrees that don't stick together that is the problem
              if( !all(tmp$family == sort(tmp$family)) )
                stop("Index operation cannot result in unordered pedigrees.") 

              if( !missing(j) )
              {
                tmp2 = tmp[,j,drop=FALSE] #Get jth element from ith row
                missingCols = c("family","id","father","mother")[!(c("family","id","father","mother") %in% names(tmp2))]

                if( length(missingCols) > 0 )
                  tmp = cbind(tmp2, tmp[,missingCols])
              }
              
              idsByFamily = split(tmp$id,tmp$family)
              familyIDs = names(idsByFamily)
              
              #index the phis
              tmpPhi = lapply(familyIDs,
                              function(fam)
                              {
                                idsForThisFam = idsByFamily[[fam]]
                                return(x@phi[[fam]][idsForThisFam,idsForThisFam])
                              })

              names(tmpPhi) = familyIDs
              
              #loop through markers and families
              #index the markers
              tmpIBD = lapply(names(x@ibd@ibdMatrix),
                              function(marker)
                              {
                                tmpIBDbyMarker = lapply(familyIDs,
                                                        function(fam)
                                                        {
                                                          idsForThisFam = idsByFamily[[fam]]
                                                          thisIbdMat = x@ibd@ibdMatrix[[marker]][[fam]]
                                                          return(thisIbdMat[idsForThisFam,idsForThisFam])
                                                        })
                                names(tmpIBDbyMarker)=familyIDs
                                return(tmpIBDbyMarker)
                              })

              names(tmpIBD) = names(x@ibd@ibdMatrix)

              x@dataVals      = tmp
              x@phi           = tmpPhi
              x@ibd@ibdMatrix = tmpIBD

            }else
            {
              stop("strumData object has invalid dataType slot.")
            }

            return(x)
          }
)

setReplaceMethod("[", "strumData",
                 function(x, i, j, value)
                 {
                   if( x@dataType == "RawData" )
                   {
                     x[i,j] <- value
                   } else if( x@dataType == "Pedigree" )
                   {
                     x[i,j] <- value
                   } else
                   {
                     stop("strumData object has invalid dataType slot.")
                   }
                 }
)

#------------------------------------------------------------------------------
# "$" generic functions
#------------------------------------------------------------------------------
setMethod("$", "strumData",
          function(x, name)
          {
            if( x@dataType == "RawData" )
            {
              x = x[[name]]
            } else if( x@dataType == "Pedigree" )
            {
              x = x[[name]]
            }else
            {
              stop("strumData object has invalid dataType slot.")
            }
          }
)

setReplaceMethod("$", "strumData",
                 function(x, name, value)
                 {
                   if( x@dataType == "RawData" )
                   {
                     x[[name]] <- value
                   } else if( x@dataType == "Pedigree" )
                   {
                     x[[name]] <- value
                   } else
                   {
                     stop("strumData object has invalid dataType slot.")
                   }
                 }
)

#------------------------------------------------------------------------------
# "[[" generic functions
#------------------------------------------------------------------------------
setMethod("[[", "strumData",
          function(x, name)
          {
            if( x@dataType == "RawData" )
            {
              x = x[[name]]
            } else if( x@dataType == "Pedigree" )
            {
              x = x[[name]]
            } else
            {
              stop("strumData object has invalid dataType slot.")
            }
          }
)

setReplaceMethod("[[", "strumData",
                 function(x, name, value)
                 {
                   if( x@dataType == "RawData" )
                   {
                     x[[name]] <- value
                   } else if( x@dataType == "Pedigree" )
                   {
                     x[[name]] <- value
                   } else
                   {
                     stop("strumData object has invalid dataType slot.")
                   }
                 }
)

#------------------------------------------------------------------------------
# "show" generic function
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumData"),
          function(object) 
          {
            if( length(object@dataVals) == 0 )
            {
              cat("Empty strumData object.\n")

            } else
            {
              value = as.character(object@dataType)
              cat("\nData type: ",value,"\n",sep="")

              if( object@dataType == "Pedigree" | object@dataType == "RawData" )
              {
                rval = nrow(object@dataVals)
                cval = ncol(object@dataVals)
                cat("Data size: ",rval," entries, ",cval," variables\n",sep="")

                if( nrow(object@dataVals) > 5 )
                  cat("\nFirst 5 rows of data values:\n")
                else
                  cat("\nData values:\n")

                print(head(object@dataVals, 5))

                if( object@dataType == "Pedigree" )
                {
                  cat("\nphi object contains ", length(object@phi)," matrices:\n")
                  cat("First matrix: \n")
                  print(head(object@phi, 1))

                  print(object@ibd)
                }

              }
            }
          }
)			

#------------------------------------------------------------------------------
# "importIBD" generic function
#------------------------------------------------------------------------------
setGeneric('importIBD', function(object,fileName,fileType="SAGE") standardGeneric('importIBD'))
setMethod("importIBD", signature(object = "strumData"),
          function(object, fileName, fileType="SAGE") 
          {
            if( file.access(fileName, mode = 0) != 0  )
            {
              cat("Importing IBD failed.  File ", fileName, " does not exist!\n")
              return(object)
            }

            if( fileType == "SAGE" )
            {
              object@ibd = .importIBDSAGE(object, fileName)
              return(object)
            } else
            {
              cat("Importing IBD failed.  Other type of IBD file will be supported in the future.")
              return(object)
            }
          }
)			

#------------------------------------------------------------------------------
# Import data from SAGE ibd file and create an IBD object.
#------------------------------------------------------------------------------
.importIBDSAGE = function(sDataObj, fileName)
{
  ibdFile = scan(fileName, what = "", sep="\n", quiet = TRUE)

  if( substr(ibdFile[1],10,12) != "2.0" )
    stop("SAGE IBD file is not in version 2.0.")

  i = 1
  while( substr(ibdFile[i],1,8) != "ANALYSIS" )
  {
    i = i+1
  }

  analysisLineNumber = i
  #cat("analysisLineNumber = ", analysisLineNumber, "\n")

  analysisBlockStartLineNumber = i + 2
  #cat("analysisBlockStartLineNumber = ", analysisBlockStartLineNumber, "\n")

  while( substr(ibdFile[i],1,7) != "MARKERS" )
  {
    i = i+1
  }

  analysisBlockEndLineNumber = i - 2
  #cat("analysisBlockEndLineNumber = ", analysisBlockEndLineNumber, "\n")

  markersLineNumber = i
  #cat("markersLineNumber = ", markersLineNumber, "\n")

  markersBlockStartLineNumber = i + 2
  #cat("markersBlockStartLineNumber = ", markersBlockStartLineNumber, "\n")

  while( substr(ibdFile[i],1,2) != "#=" )
  {
    i = i+1
  }

  markersBlockEndLineNumber = i - 1
  #cat("markersBlockEndLineNumber = ", markersBlockEndLineNumber, "\n")

  numbMarkers = markersBlockEndLineNumber - markersBlockStartLineNumber + 1
  #cat("numbMarkers = ", numbMarkers, "\n")

  ibdDataBlockStartLineNumber = markersBlockEndLineNumber + 4
  #cat("ibdDataBlockStartLineNumber = ", ibdDataBlockStartLineNumber, "\n")

  ibdDataBlockEndLineNumber = length(ibdFile)
  #cat("ibdDataBlockEndLineNumber = ", ibdDataBlockEndLineNumber, "\n")

  ibdInfo = ibdFile[(analysisBlockStartLineNumber:analysisBlockEndLineNumber)]
  #print(ibdInfo)

  markerInfo = t(sapply(ibdFile[(markersBlockStartLineNumber:markersBlockEndLineNumber)], 
                        function(x) strsplit(x," ")[[1]],
                        USE.NAMES=FALSE))
  markerInfo = data.frame(markerInfo)

  names(markerInfo) = c("Marker","Position")
  markerInfo$Marker = as.character(markerInfo$Marker)

  # Extract the pedigree ids
  #--------------------------
  pedID       = as.character(dataVals(sDataObj)$family)
  lengthPedID = length(pedID)
  pedIDs      = unique(pedID)

  # Initialize the ibd matrices to the kinship coeficient
  #-------------------------------------------------------
  ibd = list()
  for( m in markerInfo$Marker )
  {
    ibd[[m]] = sDataObj@phi
  }

  # Now start reading in the actual IBD data
  #------------------------------------------
  f = function(markerNumber,ibdNumber)
  {
    fldNumber = (markerNumber-1)*3 + ibdNumber+1
    return(c( (17+2)*(fldNumber-1)+1, (17+2)*(fldNumber)-2))
  }

  for( thisLine in ibdFile[ibdDataBlockStartLineNumber:ibdDataBlockEndLineNumber] )
  {
    fld = 0
    flds = c("pedigree","ind1","ind2")
    pair = list()
    lastDelim = 0
    for( i in 1:nchar(thisLine) )
    {
      if( substr(thisLine,i,i) == "," )
      {
        fld = fld+1
        pair[[ flds[fld] ]] = .trimSpace(substr(thisLine, lastDelim + 1, i-1))
        lastDelim = i
      }

      if( fld == 3 )
        break
    }

    startIBDCol = i+2

    for( m in 1:numbMarkers )
    {
      mName = markerInfo$Marker[m]

      # Read in f0
      # Start and end positions for marker m f0; startEnd is a vector of 2
      #--------------------------------------------------------------------
      startEnd = f(m,0) + startIBDCol

      # Check for -'s here! if -, don't set value
      #-------------------------------------------
      f0Char = substr(thisLine, startEnd[1], startEnd[2])

      if( f0Char != "-----------------" )
      {
        f0Numeric = as.numeric(f0Char)

        # Read in f2
        # Start and end positions for marker m f2
        #-----------------------------------------
        startEnd = f(m,2) + startIBDCol

        # Check for -'s here!
        #---------------------
        f2Char = substr(thisLine, startEnd[1], startEnd[2])
        if( f2Char != "-----------------" )
        {
          f2Numeric = as.numeric(f2Char)

          # Find f1
          #---------
          f1Numeric = 1 - f0Numeric - f2Numeric

          # Put ibd and f2 into the matrix
          #--------------------------------
          tryCatch({thisIbd = .5*f1Numeric + f2Numeric
                    ibd[[mName]][[pair$pedigree]][pair$ind1,pair$ind2] = thisIbd
                    ibd[[mName]][[pair$pedigree]][pair$ind2,pair$ind1] = thisIbd},
                    error = function(ex)
                            {
                              print(ex);
                              stop("Mismatch between pedigree and ibd files: ped ID = ",pair$pedigree,"\n")
                            }
                  )
        }
      }
    }
  }

  # Create object with data just parsed from file
  #-----------------------------------------------
  sIbd <- new("strumIBD",
              ibdMarker = markerInfo,
              ibdMatrix = ibd)

  .printInfoLine("Importing S.A.G.E. IBD file", "Done", 52, 0)

  return(sIbd)
}
