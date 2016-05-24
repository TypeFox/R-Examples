pause <-
function (x=0) { 
	if(x > 0){
		Sys.sleep(x)
	}else{
		cat("Hit <enter> to continue...")
		readline()
		invisible()
	}
}
