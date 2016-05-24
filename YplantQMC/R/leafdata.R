leafdata <- function(plant,...){

	
	if(plant$inputformat == "P"){
    poslen <- plant$pdata$Lt >0
		leaflen <- plant$pdata$L.3[poslen]
		petiolelen <- plant$pdata$L.2[poslen]
		an <- plant$pdata$An.3[poslen]
		az <- plant$pdata$Az.3[poslen]
		or <- plant$pdata$Or[poslen]
		Lts <- plant$pdata$Lt[poslen]
	}
	
	if(plant$inputformat == "Q"){
		leaflen <- plant$qdata$L.3
		petiolelen <- rep(0, plant$nleaves)
		an <- plant$qdata$An.3
		az <- plant$qdata$Az.3
		or <- plant$qdata$Or
		Lts <- plant$qdata$Lt
	}
  
	# Find leaf shape (may vary by leaftype).
	lsh <- c()
	for(i in 1:length(Lts))lsh[i] <- plant$ldata[[Lts[i]]]$leafshape
	
	leafarea <- leaflen^2 * lsh

	# Z coordinate of leaf base.
	# f <- function(x)unname(x$XYZ[x$midribpoints[1],3])
	# heightbase <- sapply(plant$leaves,f)
	heightbase <- plant$leafbasecoor[,3]
	
return(data.frame(leafnr=1:length(leaflen), heightbase=heightbase,
	len=leaflen, area=leafarea, 
	ang=an, az=az, or=or,
	petilen=petiolelen, 
	leafshape=lsh))

}

