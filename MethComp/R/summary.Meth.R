summary.Meth <-
function( object, ... )
{
# Make a table of no. replicates for each item
# Table of item by method
tt <- with( object, tapply( !is.na(y), list(item,meth), sum ) )
# What sort of numbers of replicates does actually exist?
xx <- as.numeric( names( table( tt ) ) )
# Add a row of each to tt so that the result of tt is rectangular
tt <- rbind( tt, matrix(xx,length(xx),ncol(tt)) )
# Subtract 1 from the result to compensate for the added rows
XX <- t(apply(tt,2,table)-1)
if( dim(XX)[1]==1 )
  {
  XX <- t(XX)
  colnames(XX) <- xx
  }
# Make column names longer to niceify output
n.spc <- ceiling((11-nchar(paste(colnames(XX),collapse="")))/length(colnames(XX)))-1
n.spc <- max( c(0,n.spc) )
colnames(XX) <- paste( paste(rep(" ",n.spc),collapse=""), colnames(XX), sep="" )
# Add the sum column
XX <- addmargins(XX,2)
# Total no. observations by method
XX <- cbind( XX,
             with( object, tapply( !is.na(y), list(meth), sum ) ),
             with( object, tapply(        y , list(meth), min ) ),
             with( object, tapply(        y , list(meth), median ) ),
             with( object, tapply(        y , list(meth), max ) )
              )
# Niceify the result
colnames( XX )[ncol(XX)-(4:0)] <- c("#Items",
                                    paste("#Obs:",sum(!is.na(object$y))),
                                    "Values:  min",
                                    "med",
                                    "max")

names(dimnames(XX)) <- c("Method"," #Replicates")
XX
}
