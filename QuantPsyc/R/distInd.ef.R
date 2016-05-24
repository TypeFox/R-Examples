"distInd.ef" <-
function (data,i) 
{ 
d <- data[i,]
ind <- as.numeric(distal.med(d)[9,1])
return(ind)
}

