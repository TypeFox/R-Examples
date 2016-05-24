rm.end.NAs <-
function(dat){
	start <- length(dat)
	for(i in start:1){
		if(is.na(dat[i])){
			end <- i
			}else{
				end <- i
				break
				}
		}
	return(dat[1:end])
	}

