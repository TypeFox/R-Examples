nnumattr <-
function(df){
# Number of numeric attributes
 pct<-0
 for (i in 1:dim(df)[2]){
  if (is.numeric(df[,i])) pct<-pct+1
 }
 pct
}

