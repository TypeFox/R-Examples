"processOrder" <- function (parorder, mod)
{
   modeldiffs <- mod@modeldiffs
   change <- modeldiffs$change 
   for(diffs in change) {
	     for(i in 1:length(parorder)) {
		   if(parorder[[i]]$name == diffs$what) {
			for(j in 1:length(diffs$dataset)) {
			   r <- which(parorder[[i]]$dataset == diffs$dataset[j])
			   if(length(r) > 0)
				parorder[[i]]$dataset <- parorder[[i]]$dataset[-r]
			    
			}   
                   }
             }
   }
   parorder
}
    