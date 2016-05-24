dataType_converter <- function(from,fromFormat="raster",toFormat="mmap")
{
	dataType_table <- list(
			raster=c(
					"LOG1S",
					"INT1S",
					"INT1U",
					"INT2S",
					"INT2U",
					"INT4S",
					"INT4U",
					"FLT4S",
					"FLT8S",
					NA,
					NA,
					NA,
					NA,
					NA,
					NA),
			mmap=list(
					as.Ctype(logi8()),
					as.Ctype(int8()),
					as.Ctype(uint8()),
					as.Ctype(int16()),
					as.Ctype(uint16()),
					as.Ctype(int32()),
					NA,
					as.Ctype(real32()),
					as.Ctype(real64()),
					as.Ctype(uchar()),
					as.Ctype(logi32()),
					as.Ctype(int24()),
					as.Ctype(uint24()),
					as.Ctype(int64()),
					as.Ctype(cplx())
			)
	)
	
	if(fromFormat=="raster")
	{
		id <- (dataType_table$raster==from) & (!is.na(dataType_table$raster))
	}
	
	if(toFormat=="mmap")
	{
		return(dataType_table$mmap[id][[1]])
	}
}