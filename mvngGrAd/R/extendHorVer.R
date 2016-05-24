extendHorVer <- function(i,
                          j,
                          shapeCross,
                          rowLimit,
                          colLimit)
   {
     ## Validity checks:
     ## 1. if shapeCross is not a list of length 4, stop and issue an error message

     if(!(is.list(shapeCross) && identical(length(shapeCross),as.integer(4))))
       {
         stop(paste("Wrong actual argument to formal argument 'shapeCross':\n",
                    "'shapeCross' must be a list of length 4!\n\n"))
       }
     ## 2. if there are duplications in the elements
     ## of the argument to shapeCross stop and issue an error message.

     anyDup <- sapply(shapeCross,duplicated)

     if(any(sapply(anyDup,any)))
       {
         stop(paste("\nWrong actual argument to formal argument 'shapeCross':\n",
                    "There must be no duplications within the elements of 'shapeCross'!\n\n"))
       }

     ## 3. if there are any '0' in any of the elements to shapeCross,
     ## stop and issue a error message.

     if(0%in%unlist(shapeCross))
       {
         stop(paste("\n\nWrong actual argument to formal argument 'shapeCross':\n",
                    "The elements of 'shapeCross' can either be NULL (no cells included\n",
                    "in this direction) or a vector of integer values, without any '0'.\n\n"))

       }
     ## 4. stop if any of shapeCross < 1

     if(any(unlist(shapeCross) < 1))
       {
         stop(paste("\n\nWrong actual argument to formal argument 'shapeCross':\n",
              "The elements of 'shapeCross' must not be negative!\n\n"))
       }
     ## 5. stop if any of the subscripts i and j exceed the row
     ## resp. col.limit or are smaller than 1

     if((as.integer(i) > as.integer(rowLimit)) |
        (as.integer(i) < 1) |
        (as.integer(j) > as.integer(colLimit)) |
        (as.integer(j) < 1))
       {
         stop(paste("\nArguments to 'i' and 'j' must be within 1:rowLimit resp.\n",
                    "1:colLimit\n\n"))
       }


     ## Horizontal: for each value in shapeCross
     ## (l) go i-l up respective i+1 down. j
     ## remains the same since the column
     ## doesn't change. Test whether the row
     ## subscript gets out of bounds, if so,
     ## both the row and column subscript
     ## remain NULL. This makes sure that no
     ## subscripting takes place latter
     ## on. To be sure and to save space,
     ## remove NULL subscripts at the end. 
     ## If the value given for shapeCross is NULL,
     ## don't do anything.
     ## Vertical: principally the same, just
     ## that j changes towards right or left
     ## and i remains constant.
     includedDown <- matrix(0,nrow=length(shapeCross[[1]]),
                            ncol=2)     
     pos <- 1
     if(!is.null(shapeCross[[1]]))
       {
         for(l in shapeCross[[1]])
           {
             if(i+l <= rowLimit)
               {
                 includedDown[pos,] <-
                   c(i+l,j)
               }
             pos <- pos + 1
           } ## end l
       }
     
     includedUp <- matrix(0,nrow=length(shapeCross[[2]]),
                          ncol=2)
     ##
     pos <- 1
     if(!is.null(shapeCross[[2]]))
       {
         for(l in shapeCross[[2]])
           {
             if(i-l >= 1)
               {
                 includedUp[pos,] <-
                   c(i-l,j)
               }
             pos <- pos + 1
           } ## end l
       }
     ## j + l > colLimit, NULL
     ## i - l < 1, NULL

     includedLeft <- matrix(0,nrow= length(shapeCross[[3]]),
                            ncol=2)
     pos <- 1
     if(!is.null(shapeCross[[3]]))
       {
         for(l in shapeCross[[3]])
           {
     
             if(j - l >= 1)
               {
                 includedLeft[pos,] <-
                   c(i, j - l)
               }
             pos <- pos + 1
           } ## end l
       }
     
     includedRight <- matrix(0,nrow=length(shapeCross[[4]]),
                             ncol=2)
     pos <- 1
     if(!is.null(shapeCross[[4]]))
       {
         for(l in shapeCross[[4]])
           {
     
             if(j + l <= colLimit)
               {
                 includedRight[pos,] <- 
                   c(i,j + l)
               }
             pos <- pos + 1
           } ## end
       }
     verHor <- rbind(includedDown,
                     includedUp,
                     includedLeft,
                     includedRight)
     verHor <- verHor[which(verHor[,1] != 0),,drop=FALSE]     
     return(verHor)
   } ## end of function
