analyzeBits <-
function(value,mode=1, bitpos=0){
	if (is.na(value)){
		return(NA)
	}
	
	if ((bitpos < 0)||(bitpos > 15)){
		stop("bitpos has to be between 0 and 15")
	}
	
	mask <- 0

	if (mode == 0)
		mask <- .C("getBits", number=as.integer(value), pos=as.integer(bitpos), m=as.integer(mask), PACKAGE="phenex")$lw;
	if (mode == 1)
		mask <- .C("getLandWater", number=as.integer(value), lw=as.integer(mask), PACKAGE="phenex")$lw;
	if (mode == 2)
		mask <- .C("getCloudmask", number=as.integer(value), cm=as.integer(mask), PACKAGE="phenex")$cm;
	if (mode == 3)
		mask <- .C("getDayNo", number=as.integer(value), day=as.integer(mask), PACKAGE="phenex")$day;
	
	return(mask)
}
