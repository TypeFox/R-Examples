

csv_from_xes <- function(xesfile) {

	parsed_xes <- parseXES(xesfile)

	n_case_att <- length(parsed_xes$trace.att)
	n_event_att <- length(parsed_xes$event.att)
	n_cases <- length(parsed_xes$traces$'concept:name')
	caseids <- parsed_xes$traces$'concept:name'
	result <- data.frame(stringsAsFactors = FALSE)

	n<-1
	for(i in 1:n_cases){
		data <- add_NAS_to_event_attributes_per_case(parsed_xes$events[[i]])
		n_events <- length(data[[1]])
		case <- as.data.frame(data, stringsAsFactors = F)
		case <- bind_cols(data.frame(rep(caseids[i], n_events), stringsAsFactors = F), case)
		result <- bind_rows(result, case)

	}
	colnames(result)[1]<-"case_concept.name"
	for(i in 1:n_event_att)
		colnames(result)[i+1] <- paste("event",colnames(result)[i+1], sep = "_")
	return(result)
}

add_NAS_to_event_attributes_per_case <- function(x) {
	for(i in 1:length(x))
		v <- (as.vector(sapply(x, length)))
	if(length(unique(v)) > 1){
		# meer dan 1 waarde
		m <- max(v) # max aantal attributen
		for(i in 1:length(x)){
			if(length(x[[i]]) < m){
				x[[i]] <- c(x[[i]], rep(NA, times = (m - length(x[[i]]))))
			}
		}
	}
	return(x)
}



##Handler function for XES parsing.Used by parseXES
handler <- function(){
	#states: log, trace, event
	state <- "log"
	trace.data <- list()
	event.data <- list()
	trace.counter <- 0
	event.counter <- 0


	trace <- function(x,atts){
		state <<- "trace"
		trace.counter <<- trace.counter + 1
		event.data[[trace.counter]] <<- list()
	}

	event <- function(x,atts){
		state <<- "event"
		event.counter <<- event.counter + 1
	}

	endElement <- function(x,...){
		if(x =="trace"){
			state <<- "log"
			event.counter <<- 0
		}
		else if(x == "event"){
			state <<- "log"
		}
	}

	attributes <- function(x,atts){
		if(state =="trace"){
			trace.data[[atts[["key"]]]][trace.counter] <<- atts[["value"]]
		}
		else if (state == "event"){
			event.data[[trace.counter]][[atts[["key"]]]][event.counter] <<- atts[["value"]]
		}
	}


	return(list(
		trace = trace,
		event = event,
		endElement = endElement,
		date = attributes,
		string = attributes,
		int = attributes,
		float = attributes,
		boolean = attributes,
		trace.data = function(){trace.data},
		event.data = function(){event.data}
	))

}

parseXES <- function(logfile){
	temp <- XML::xmlEventParse(logfile,handler())
	tracedata <- temp$trace.data()
	trace.att <- names(tracedata)
	eventdata <- temp$event.data()
	event.att <- unique(unlist(lapply(eventdata,function(x){names(x)})))
	return(list(
		traces = tracedata,
		events = eventdata,
		trace.att = trace.att,
		event.att = event.att
	))
}
