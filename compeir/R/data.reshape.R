data.reshape <-
function (
                        times,
                        events,
                        covar,
                        no.event.code,
                        event.code=NULL,
                        covar.code=NULL
                        ){

  events = unlist(events)
  covar = unlist(covar)
  times = unlist(times)
  
  no.event.code = as.character(no.event.code)
  
  if(is.null(event.code)){
  	tmp = levels(as.factor(events))
  	event.code = tmp[-which(tmp == no.event.code)]
  	}
  if(is.null(covar.code)) covar.code = levels(as.factor(covar))
  
  event.code =  as.character(event.code)
  covar.code =  as.character(covar.code)
  
  if(any(!(covar.code %in% covar))){stop(paste("covar.code", paste(covar.code[which(!(covar.code %in% covar))], collapse = ", "),"is not contained in covar"))}

  if(any(!(event.code %in% events))){stop(paste("event.code", paste(event.code[which(!(event.code %in% events))], collapse = ", "),"is not contained in events"))}

  time.code = "time"
  names <- c(time.code,no.event.code,event.code)
  row.names <- covar.code

  data <- as.data.frame(matrix(NA, ncol=length(names), nrow=length(row.names)))
  names(data) <- names
  row.names(data) <- row.names
		
  for(i in row.names){
  	data[i, time.code] = sum ( as.numeric(times[which(covar %in% i & events %in% c(no.event.code, event.code))]) )
  	
  for ( j in c(no.event.code, event.code) ){
    data[i,j] = length(events[which(events %in% j & covar %in% i)])
    }
  }

  class(data) <- "data.frame"

  return (data)
}

