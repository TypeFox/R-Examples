`compnames` <-
function(nobj){
#################################################################
#  function to generate a vector of labels
#  for comparisons in sinclair ordering
#################################################################
     nobj  <-get("nobj",get("ENV", environment(patt.design)))
     str<-NULL
     for (j in 2:nobj) {
         for (i in 1:(j-1) ){
             str<-c(str, paste("u",i,j,sep=""))
         }
     }
     assign("undecnames",str,get("ENV", environment(patt.design)))
     str
}
