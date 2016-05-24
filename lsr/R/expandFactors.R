# file:    expandFactors.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# expandFactors() takes a data frame as input, and returns the same data frame,
# but with all factor variables replaced by the corresponding contrasts. It's
# actually just a wrapper to model.matrix()
expandFactors <- function( data, ... ) {
  
  attr(data,"na.action") <- na.pass # don't drop NA
  df <- model.matrix( as.formula( paste("~",names(data),collapse="+")), data, ... )
  print(df)
  df <- df[,-1,drop=FALSE] # remove intercept
  attr(df,"contrasts") <- NULL
  attr(df,"assign") <- NULL
  
  return( as.data.frame(df) )
}
