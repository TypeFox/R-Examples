"put.snps" <-
function(idata) .C("put_snps",as.integer(idata), as.integer(length(idata)), cdata = raw(ceiling(length(idata)/4)), PACKAGE="GenABEL")$cdata

