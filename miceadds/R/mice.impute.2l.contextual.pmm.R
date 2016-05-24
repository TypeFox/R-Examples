mice.impute.2l.contextual.pmm <- function (y, ry, x, type , imputationWeights = NULL , 
                interactions=NULL , quadratics = NULL , ...){

   vname <- get("vname", pos = parent.frame()) # get variable name            
   newstate <- get( "newstate" , pos = parent.frame() )  
   # data preparation
   xcov <- .a2l.contextual.auxiliary( y = y  , ry=ry , x=x , type=type , ...)
   #------
   # pmm imputation at level 2
#   print( cat("\n"))
#print( head(xcov))
   ximp <- mice.impute.weighted.pmm( y= y , ry=ry, x = xcov , imputationWeights = imputationWeights , 
            interactions= interactions , quadratics = quadratics ,   ... ) 
   return(ximp)
    }
