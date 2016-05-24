"UPmultinomial" <-
function(pik) 
{if(any(is.na(pik))) stop("there are missing values in the pik vector")
as.vector(rmultinom(1,sum(pik),pik/sum(pik)))
}

