# Frequency count and percentage table
# Table 9.40 :  Hilbe, JM (2011) Negative Binomial Regression, 2nd ed, Cambridge Univ Press
myTable <- function(x) { 
myDF <- data.frame( table(x) ) 
myDF$Prop <- prop.table( myDF$Freq ) 
myDF$CumProp <- cumsum( myDF$Prop ) 
myDF 
} 









