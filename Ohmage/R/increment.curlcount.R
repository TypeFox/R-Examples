increment.curlcount <- function(){
	mycount <- getOption("CURLCOUNT");
	if(mycount < 1000){
		options(CURLCOUNT = mycount + 1);
	} else {
		message("Reached 1000 curl requests. Renewing curl handler!")
		renewCurl();
	}
}
