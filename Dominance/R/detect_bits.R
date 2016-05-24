detect_bits <-
function(bits,set=TRUE)
{
# ----------------------- detect other char than 0 and 1 -> error -------------
# TRUE if all 0 and 1
all_bits <- (regexpr("^[01]*$", bits) > 0) 
if (all_bits ==FALSE)
	{
	print("errror in detect_bits: variable bits must contain only \"0\" and \"1\"",quote=FALSE)
		{
		 return(all_bits)
		 break;
		}
	}#end if

# end detect other char than 0 and 1 -> error
if (set==TRUE) return(match <- gregexpr("1", bits)[[1]])
if (set==FALSE)return(match <- gregexpr("0", bits)[[1]])

}
