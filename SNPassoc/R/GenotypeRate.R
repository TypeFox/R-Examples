`GenotypeRate` <-
function(x)
{
 temp<-sum(!is.na(x))/length(x)
 ans<-temp*100
 ans
}

