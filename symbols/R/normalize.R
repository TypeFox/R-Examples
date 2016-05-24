normalize <-
function(ar){
# Normalization into <0,1>
 # ar<-(ar-mean(ar))/sd(ar)
 if (is.numeric(ar)) {ar<-(ar-min(ar))/(max(ar)-min(ar))}
 ar
}

