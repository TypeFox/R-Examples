print.irates.ratio <-
function( 	x, 
									event.code = NULL,
									display.digits = 4,
									... ){
										
  if(class(x) != "irates.ratio"){stop("Object needs to be of class irates.ratio")}
  object = x
  if( is.null(event.code) ){
  	event.code = object$event.code
  	event.lab = object$event.lab
  	} else{
  		event.lab = object$event.lab[which(object$event.code %in% event.code)]
  		}  	

  event.code = as.character(event.code)

  if(any(!(event.code %in% object$event.code))){
	stop(paste("event.code *", event.code[which(!(event.code %in% object$event.code))], "* is not contained in irates.ratio object"))
	}
	
  covar.lab <- object$covar.lab

  n.event <- length(event.lab)
  
  cat("\n")
  cat(paste(covar.lab[1],":",covar.lab[2]))

  cat("\n\n")
    
  text <- matrix(NA, ncol=4, nrow = n.event)
  text <- as.data.frame(text)
    
  names(text) <- c(
                   "irr",
                   paste("lower .",object$ci.level*100, sep=""),
                   paste("upper .",object$ci.level*100, sep=""),
                   "var(irr)")
  row.names(text) <- event.code
      
  for ( i in event.code ){
    text[i,1] <- round(object$irr[i],display.digits)
    text[i,2] <- round(object$conf.lower[i],display.digits)
    text[i,3] <- round(object$conf.upper[i],display.digits)
    text[i,4] <- format(object$var[i],digits=display.digits)
  }

  row.names(text) <- event.lab

  print(text, quote=FALSE, right=FALSE, ...)
  cat("\n")

  cat("\n")
  cat(paste("ci.fun = \"",object$ci.fun,"\"\n",sep=""))
}

