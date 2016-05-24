library("pdc")

if (!all(codebook(c(1,2,3,4,5),m=3,t=1)==c(0,0,0,0,0,1))) stop("test failed")

if (!all(codebook(c(1,2,3,4,5,4,3,2,1),m=2)==c(0.5,0.5))) stop("test failed")

