"sset" <-
function(data,nsnps,nids,list) .C("sset",as.raw(data), as.integer(nsnps), as.integer(nids), as.integer(list), as.integer(length(list)), cdata = raw(nsnps*ceiling(length(list)/4)), PACKAGE="GenABEL")$cdata

