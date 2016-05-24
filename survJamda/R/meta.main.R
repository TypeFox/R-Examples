meta.main <-
function(geno.files,surv.data, method = "none")
{	
	options(warn=-1) 
	curr_set = 1:length(geno.files)
	for (y in curr_set){
		x = setdiff(curr_set, y)
		det.set.meta (x, y, geno.files,surv.data, method)
	}
}

