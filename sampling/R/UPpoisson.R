"UPpoisson" <-
function(pik) 
{if(any(is.na(pik))) stop("there are missing values in the pik vector")
as.numeric(runif(length(pik))<pik)
}
