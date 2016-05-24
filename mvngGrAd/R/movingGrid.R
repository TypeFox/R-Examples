movingGrid <- function(rows,
                       columns,
                       obsPhe,
                       shapeCross,
                       layers,
                       excludeCenter = TRUE)
  {
    FUNcall <- list(quote(movingGrid),
                    rows = rows,
                    columns = columns,
                    obsPhe = obsPhe,
                    shapeCross = shapeCross,
                    layers = layers,
                    excludeCenter = excludeCenter)
    
    FUNcall <- as.call(FUNcall)
    
    ## maximum possible number of values, given the choosen design
    
    plusCenter <- 0
    
    if(identical(excludeCenter,FALSE))
      {
        plusCenter <- 1
      }
    
    maxValues <- sum(length(unlist(shapeCross)), ## cross
                      sum(layers*4),         ## circular layers
                      plusCenter)

     ## validity checks
##################

     ## 1. stop if the grid contains no values

     if(identical(as.integer(maxValues), as.integer(0)))
       {
         stop(paste("\nThe grid must consist of at least one value (cell)!\n\n"))
       }

     ## 2. give a warning if the grid consists only of the center cell

     if(identical(as.integer(maxValues), as.integer(1)) &&
        identical(excludeCenter,FALSE))
       {
         warning(paste("\nThe grid consisted only of the center cell. This makes\n",
                       "no sense and the adjusted values will give no information at all!\n\n"))
       }
    

     ## calculate the row and column limits of the    
     ## field map:
    
     rowLimit <- max(rows)
     colLimit <- max(columns)

     ## bring the row and col subscripts and the
     ## obsPhe into one matrix
     rowCol <- cbind(rows,columns,obsPhe)
     ##
     ## empty matrix with max(rowCol[,1]) rows and
     ## max(rowCol[,2]), columns
     ##
     fieldMap <- matrix(nrow = rowLimit,
                        ncol = colLimit)
     ##
     ## fill the matrix with rowCol[,3] by
     ## indexing via the rowCol[,1:2] matrix.
     ##
     fieldMap[rowCol[,1:2]] <- rowCol[,3]
     ##
    
     ## to store the moving means
     movingMeanMap <- matrix(nrow = rowLimit,
                             ncol = colLimit)

     ## and the number of included (non NA) values

     nValuesMap <- matrix(nrow = rowLimit,
                          ncol = colLimit)

     
                          
     for(i in 1:rowLimit) ## rows
       {
         for(j in 1:colLimit) ## columns
           {
             ## coordinates of the cells in the "cross"
             
             coordCross <- extendHorVer(i = i,
                                        j = j,
                                        shapeCross = shapeCross,
                                        rowLimit = rowLimit,
                                        colLimit = colLimit)

             ## coordinates of the cells in the circular layers

             coordCirc <- circularExtension(i = i,
                                            j = j,
                                            layers = layers,
                                            rowLimit = rowLimit,
                                            colLimit = colLimit)

             ## coordinates of all the cells in the grid:

             ## excluding the center:

             if(identical(excludeCenter,FALSE))
               {
                 coordGrid <- rbind(cbind(i,j),
                                    coordCross,
                                    coordCirc)
               }

             else
               {
                 coordGrid <- rbind(coordCross,
                                    coordCirc)
               }
             

             ## nValues
             nValuesMap[i,j] <-
               sum(!is.na(fieldMap[coordGrid[,1:2,drop = FALSE]]))
             

             ## moving mean calculation ########
             ## ################################

             ## To prevent the occurence of NaNs, calculate the mean
             ## only if there is at least one non-NA value, else set
             ## the moving Mean to NA
             
             if(sum(!is.na(fieldMap[coordGrid[,1:2,drop = FALSE]])) > 0)
               {
                 movingMeanMap[i,j] <-
                   mean(fieldMap[coordGrid[,1:2,drop = FALSE]],
                        na.rm = TRUE)
               }
             
             else
               {
                 movingMeanMap[i,j] <- NA
               }

           } ## end j
             
       } ## end i
     


     ## adjustment procedure ##
     ## #################### ##

     ## 1. create a vector with the moving mean results that,
     ## corresponds to row, column and obsPhe

     movingMeans <- movingMeanMap[rowCol[,1:2,drop = FALSE]]

     ## and for the nValues

     nValues <- nValuesMap[rowCol[,1:2,drop = FALSE]]

     ## 2. calculate the coefficient of the regression of
     ## the moving means on the phenotypic values
     
     model <- lm(rowCol[,3] ~ movingMeans,
                 na.action = "na.exclude")

     coefAdj <- model$coefficients[2]

     names(coefAdj) <- NULL

     ## 3. calculate the adjusted phe. values
     
     adjPhe <-
       rowCol[,3] -
         (coefAdj*(movingMeans - mean(movingMeans,na.rm = TRUE)))

     adjPhe <- round(adjPhe,3)

     ## 4. calculate the correlation between the moving means
     ## and the phenot. values

     correlation <-
       cor(rowCol[,3],movingMeans, use = "complete.obs")


### give a warning and suggestions when moving means could not be
### calculated for all entries

     if(any(is.na(movingMeans)))
       {
         warning(paste("\n\nMoving means could not be calculated for all entries\n",
                       "because no non-NA values were available.\n",
                       "A solution is to change or enlarge the design of the grid.\n\n"))
       }

    
     result <- new("movG",
                   movingMeanMap = movingMeanMap,
                   row = as.integer(rowCol[,1]),
                   col = as.integer(rowCol[,2]),
                   observedPhe = rowCol[,3],
                   adjustedPhe = adjPhe,
                   movingMean = movingMeans,
                   nValues = nValues,
                   adjModel = model,
                   correlation = correlation,
                   maxValues = as.integer(maxValues),
                   FunCall = FUNcall)
     
     return(result)

  }
