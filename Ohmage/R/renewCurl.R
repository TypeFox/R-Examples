renewCurl <- function(){
	#unlockBinding("OhCurlHandle", as.environment("package:Ohmage"));
	#unlockBinding("OhCurlReader", as.environment("package:Ohmage"));
	#assign("OhCurlHandle", getCurlHandle(), as.environment("package:Ohmage"));
	#assign("OhCurlReader", dynCurlReader(OhCurlHandle, binary = TRUE), as.environment("package:Ohmage"));
	#lockBinding("OhCurlHandle", as.environment("package:Ohmage"));
	#lockBinding("OhCurlReader", as.environment("package:Ohmage"));
	
	options(OhCurlHandle = getCurlHandle());
	options(OhCurlReader = dynCurlReader(getOption("OhCurlHandle"), binary = TRUE));
	options("CURLCOUNT"= 0);
}
