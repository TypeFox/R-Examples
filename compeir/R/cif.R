cif <-
function (irates, 
					t, 
					event.code=NULL, 
          covar.code=NULL,
          full.sample = FALSE
          ){
                   	
  object = irates
   
  given.covar.code = c(object$covar.code, object$full.sample.code)
  given.covar.lab = c(object$covar.lab, object$full.sample.lab)

  if(any(t<0)) stop("t needs to be >= 0")

  if(is.null(event.code)) event.code = object$event.code
  if(is.null(covar.code)) covar.code = object$covar.code
  
  event.code = as.character(event.code)
  covar.code = as.character(covar.code)
  
  if(full.sample){
  	covar.code = c(covar.code, object$full.sample.code)
  	}
  
  if(any(!(covar.code %in% given.covar.code))){
  	stop("covar.code", paste(covar.code[which(!(covar.code %in% given.covar.code))], collapse = ", "), "is not contained in irates object")
  	}

  if(any(!(event.code %in% object$event.code))){
  	stop("event.code", paste(event.code[which(!(event.code %in% object$event.code))], collapse = ", "), "is not contained in irates object")
  	}
  	
  temp = vector("list", length(event.code))
  
  for(i in 1:length(event.code)){
  	temp[[i]] = as.data.frame(object$ir[covar.code,event.code[i]]/(rowSums(object$ir[covar.code,]))*(1-exp(-rowSums(object$ir[covar.code,]) %o% t)), row.names = covar.code)
  	names(temp[[i]]) = t
  	}	
  	
  	names(temp) = object$event.code[which(event.code == object$event.code)]  	
  res <- list(
              t=t,
              cif=temp,
              event.code=event.code,
              covar.code=covar.code
              )
              
  class(res) <- "cif"
  return (res)
}

