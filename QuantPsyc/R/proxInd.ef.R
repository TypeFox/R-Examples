"proxInd.ef" <-
function (data,i) 
{ 
d <- data[i,]
ind <- as.numeric(proximal.med(d)[5,1])
return(ind)
}

