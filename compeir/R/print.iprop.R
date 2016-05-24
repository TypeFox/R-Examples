print.iprop <-
function( 	x,
							event.code = NULL,
							covar.code = NULL, 
							full.sample = FALSE,
							display.digits = 4,
							... 
							){
								
  if(class(x) != "iprop"){stop("Object needs to be of class iprop")}
  object = x
  given.covar.code = c(object$covar.code, object$full.sample.code)

  if(is.null(event.code)) event.code = object$event.code
  if(is.null(covar.code)) covar.code = object$covar.code
  
  event.code = as.character(event.code)
  covar.code = as.character(covar.code)

  covar.labels = rep(NA, length(covar.code))
  for(l in 1:length(covar.code)){
  	covar.labels[l] = object$covar.lab[which(object$covar.code == covar.code[l])]
    }

  if(full.sample){
  	covar.code = c(covar.code, object$full.sample.code)
  	covar.labels = c(covar.labels, object$full.sample.lab)
  	}  
  				
  if(any(!(covar.code %in% given.covar.code))){
  	stop("covar.code *", covar.code[which(!(covar.code %in% given.covar.code))], "* is not contained in iprop object")
  	} 
  if(any(!(event.code %in% object$event.code))){
  	stop("event.code *", event.code[which(!(event.code %in% object$event.code))], "* is not contained in iprop object")
  	} 
  	  
  for ( k in covar.code ){

    cat("\n")
    cat(covar.labels[which(covar.code == k)])

    cat("\n\n")
    
    text <- matrix(NA, ncol=4, nrow=length(event.code))
    text <- as.data.frame(text)
    
    names(text) <- c(
                     "ip",
                     paste("lower .",object$ci.level*100, sep=""),
                     paste("upper .",object$ci.level*100, sep=""),
                     "var(ip)")
    row.names(text) <- event.code
      
    for ( i in event.code ){
      text[i,1] <- round(object$ip[k,i],display.digits)
      text[i,2] <- round(object$conf.lower[k,i],display.digits)
      text[i,3] <- round(object$conf.upper[k,i],display.digits)
      text[i,4] <- format(object$var[k,i],digits=display.digits)
    }
    #row.names(text) <- object$covar.lab[which(object$covar.code %in% covar.code)]
    
    event.labels = rep(NA, length(event.code))
	for(l in 1:length(event.code)){
		event.labels[l] = object$event.lab[which(object$event.code == event.code[l])]
    }
    row.names(text) = event.labels
 	#row.names(text) = object$event.lab[which(object$event.code %in% event.code)]


    print(text, quote=FALSE, right=FALSE, ...)
    cat("\n")
  }
  cat("\n")
  cat(paste("ci.fun = \"",object$ci.fun,"\"\n",sep=""))
}

