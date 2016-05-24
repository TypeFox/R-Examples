circularExtension <- function(i,
                               j,
                               layers,
                               rowLimit,
                               colLimit){
  ## validity checks:

  ## 1. if any of the layers is present more than once in the the
  ## argument to layers, "stop" and issue an error message.

  if(any(duplicated(layers)))
    {
      stop(paste("\nWrong actual argument to formal argument 'layers':\n",
                 "Each layer must only appear once!\n\n"))
    }

  ## 2. if layers is not NULL but a "0" appears in the vector, stop
  ## and issue an error message

  if(0%in%layers)
    {
      stop(paste("\nWrong actual argument to formal argument 'layers':\n",
                 "'layers' can either be 'NULL' (no circular layer at all) or a vector\n",
                 " with the layers given as integers, without any '0'!.\n\n")) 
    }

  ## 3. stop if any of layers < 1

  if(any(layers < 1))
    {
      stop(paste("\n\nWrong actual argument to formal argument 'layers':\n",
                 "There must be no negative layers!\n\n"))
    }
  
  ## 4. stop if any of the subscripts i and j exceed the row
  ## resp. col.limit or are smaller than 1

  if((as.integer(i) > as.integer(rowLimit)) |
     (as.integer(i) < 1) |
     (as.integer(j) > as.integer(colLimit)) |
     (as.integer(j) < 1))
    {
      stop(paste("\nArguments to 'i' and 'j' must be within 1:rowLimit resp.\n",
                 "1:colLimit\n\n"))
    }

  diagSubscripts <- matrix(nrow=0,ncol=2)

  ## if layers = NULL, do nothing and return an empty matrix

  if(is.null(layers))
    {
      return(diagSubscripts)
    }

  else
    {
      for(n in layers)
        {
          ## i starts at zero o'clock as i - n and rises in
          ## steps of 1 to i + n at six o'clock (excluding the
          ## i of the center).
          ## 
          zeroToSix <- c((i-n):(i-1),
                         (i+1):(i+n))
          ## Then at six o'clock i falls from i + n in steps
          ## of 1 to i-n at zero o'clock (exluding the i of the
          ## center). This is simply the opposite of zeroToSix.
          sixToZero <- zeroToSix[length(zeroToSix):1]
          ##
          wholeI <- c(zeroToSix,sixToZero)
          ##
          ## j starts at nine o'clock as j - n and rises in
          ## steps of 1 to j + n at three o'clock (excluding
          ## the j of the center).
          ## 
          nineToThree <- c((j-n):(j-1),
                           (j+1):(j+n))
          ## Then at three o'clock j falls from j + n in
          ## steps of 1 to j-n at nine o'clock (exluding the
          ## j of the center). This is simply the opposite
          ## of nineToThree.
          threeToNine <- nineToThree[length(nineToThree):1]
          ##
          wholeJ <- c(nineToThree, threeToNine)
          ##
          ## wholeJ still needs to be tuned with wholeI so
          ## that it to goes from zero to zero o'clock 
          wholeJ <- wholeJ[c(n+1:(length(wholeJ)-n),1:n)]
          ##
          ## Prevention of indexing when subscript out of
          ## bounds: if the subscript in wholeI is < i, it
          ## means that it goes up toward row = 1. Therefore
          ## it must not be smaller than 1. If the subscript
          ## is > i, it means that it goes down toward row =
          ## rowLimit so it must not be larger.
          ##
          wrongI <- which((wholeI < 1)
                          | (wholeI > rowLimit))
          ## if one is TRUE, subscripting is prevented by
          ## setting i and the corresponding j to NULL.
          wholeI[wrongI] <- is.null(wholeI[wrongI])
          wholeJ[wrongI] <- is.null(wholeJ[wrongI])
          ##
          ## The same for wholeJ, just that the subscript
          ## must not be smaller than 1 or larger than
          ## colLimit. Subscripts already set to NULL
          ## are not tested again.
          ##
          wrongJ <- which(!(is.null(wholeJ))
                          & (wholeJ < 1)
                          |(wholeJ > colLimit))               
          ##
          wholeJ[wrongJ] <- is.null(wholeJ[wrongJ])
          wholeI[wrongJ] <- is.null(wholeI[wrongJ])
          ##
          diagSubscripts <- rbind(diagSubscripts,
                                  cbind(wholeI, wholeJ))
        } ## end n
      diagSubscripts <- diagSubscripts[which(diagSubscripts[,1] != 0),,drop=FALSE]
      return(diagSubscripts)
    }
} ## end of function


