"outlier" <-
function (x, opposite = FALSE, logical = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, outlier, opposite = opposite, logical = logical)
    else if (is.data.frame(x)) 
        sapply(x, outlier, opposite = opposite, logical = logical)
    else {
	if (xor(((max(x,na.rm=TRUE) - mean(x,na.rm=TRUE)) < (mean(x,na.rm=TRUE) - min(x,na.rm=TRUE))),opposite)) 
		{
			if (!logical) min(x,na.rm=TRUE)
			else x == min(x,na.rm=TRUE)
		}
		else 
		{
			if (!logical) max(x,na.rm=TRUE)
			else x == max(x,na.rm=TRUE)
		}
	} 
}

