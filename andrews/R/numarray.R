numarray <-
function(df){
# Numeric attributes only
 dfm<-0
 m<-dim(df)[2]
 for (a in 1:m){
  if (is.numeric(df[,a])){
   if (length(dfm)>1) dfm<-cbind(dfm,df[,a]) else dfm<-df[,a]
  }
 }
 dfm
}

