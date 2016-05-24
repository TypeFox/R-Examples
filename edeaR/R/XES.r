
createXES<- function(file,
                     traces,
                     events,
                     classifiers = NULL,
                     logattributes = NULL,
                     caseid_field = NULL){
  # File: The location of the output file
  # Traces: a dataframe where each row represents a trace and each column represents
  #         a trace attribute.
  # Events: a dataframe where each row represents an event and each column represents
  #         an event attribute. (This dataframe also has a column which refers to the
  #         traceid)
  # Classifiers: A list of classifiers. The key represents the name of the classifier
  #              and the value contains a string vector of the respective event attributes
  # Logatrributes:  A list of log atributes. The key represents the attribute name and the
  #                 value represents the attribute value. The attribute type is derived from
  #                 the attribute value
  # Caseid_field: The columnname which acts as traceid in the events dataframe.
  #               DEFAULT: The first column of the events dataframe is used as traceID.



  ##################HELPER FUNCTIONS############################
  add <- function(text){
    xml_i <<- xml_i + 1
    xml[xml_i] <<- text
    #xml
    #cat(text, file, append=TRUE)
  }

  addAttribute <- function(datatype, key, value){
    if(is.null(value)){
      return()
    }
    if(datatype == "date"){
      value = strftime(value,format="%Y-%m-%dT%H:%M:%S.000+00:00")
    }
    add(paste0('<',datatype,' key="',key,'" value="',value,'"/>'))
  }

  addExtensions <- function(attrs){
    #add concept extension if any of the event attributes start with the "concept:" prefix
    if(any(grepl("^concept:", names(attrs)))){
      add('<extension name="Concept" prefix="concept" uri="http://www.xes-standard.org/concept.xesext" />')
    }
    #add time extension if any of the event attributes start with the "time:" prefix
    if(any(grepl("^time:", names(attrs)))){
      add('<extension name="Time" prefix="concept" uri="http://www.xes-standard.org/time.xesext" />')
    }
  }

  addGlobals <- function(data, attrs, scope){
    temp <- sapply(data[,names(attrs)],function(x){all(!is.na(x))})
    globals <- names(temp[temp==TRUE])

    add(paste0('<global scope="',scope,'">'))
      for(key in globals){
        datatype <- attrs[key]
        addAttribute(datatype, key, defaultvalues[[datatype]])
      }
    add('</global>')
  }

  addClassifiers <- function(){
    if(is.null(classifiers)){
      return()
    }
    for(name in names(classifiers)){
      add(paste0('<classifier name="',name,'" keys="', paste(classifiers[[name]], collapse=" "),'"/>'))
    }
  }

  addLogAttributes <- function(){
    if(is.null(logattributes)){
      return()
    }
    for(name in names(logattributes)){
      value = logattributes[[name]]
      datatype = attribute_types[class(value)[1]]
      if(datatype == "date"){
        value = strftime(value,format="%Y-%m-%dT%H:%M:%S.000+00:00")
      }
      add(paste0('<', datatype,' key="',name,'" value="', value,'"/>'))
    }
  }

  addTraces <- function(){
#    apply(traces, 1, addTrace)
    total = dim(traces)[1]
    pb <- txtProgressBar(min=0, max = total, style = 3)
    for(i in 1:total){
      trace <- traces[i,,drop=FALSE]
      addTrace(trace)
      setTxtProgressBar(pb, i)
    }
  }

  addTrace <- function(trace){
    add('<trace>')
    for(name in names(trace_attrs)){
      addAttribute(trace_attrs[name], name, trace[name])
    }
    trace_id <- as.character(trace[[trace_caseid_field]])
    addEvents(events_per_trace[[trace_id]])
    add('</trace>')
  }

  addEvents <- function(trace_events){
    apply(trace_events, 1, addEvent)
  }

  addEvent <- function(event){
    add('<event>')
    for(name in names(event_attrs)){
      addAttribute(event_attrs[name], name, event[name])
    }
    add('</event>')
  }

  detectAttrType <- function(key, data){
    detected = class(data[[key]])[1]
    attribute_types[[detected]]
  }

  get_attr_info<- function(data){
    sapply(colnames(data), detectAttrType, data)
  }


  ############PRELIMINARIES##################
  defaultvalues <- list("string"="default",
                        "int"="0",
                        "date"="1970-01-01T00:00:00.000+00:00")
  attribute_types <- list("factor"="string",
                          "POSIXct"="date",
                          "integer"="int",
                          "ordered"="string",
                          "character"="string")
  trace_attrs <- get_attr_info(traces)

  if(is.null(caseid_field)){
    event_attrs <- get_attr_info(events[,-1])
  }
  else{
    event_attrs <- get_attr_info(events[,names(events)!=caseid_field])
  }

  if(is.null(caseid_field)){
    events_caseid_field <- colnames(events)[1]
  }
  else { #adjustment Gert
  	events_caseid_field <- caseid_field
  }

  events_per_trace <- split(events,events[events_caseid_field])

  if("concept:name" %in% names(trace_attrs)){
    trace_caseid_field <- "concept:name"
  }
  else if(events_caseid_field %in% names(trace_attrs)){
    trace_caseid_field <- events_caseid_field
  }
  else{
    trace_caseid_field <- colnames(traces)[1]
  }

  n_event_attrs = length(event_attrs)
  n_trace_attrs = length(trace_attrs)
  n_classifiers = length(classifiers)
  n_logattributes = length(logattributes)
  n_traces = dim(traces)[1]
  n_events = dim(events)[1]
  maxsize = 4+ 3*n_event_attrs + 3*n_trace_attrs + n_classifiers + n_logattributes + n_traces*(2+n_trace_attrs)+n_events*(2+n_event_attrs)
  xml <- rep(NA,maxsize)
  xml_i = 1
  ############GENERATE XML###################
  #fileConn <- file(file, open="at")
  add('<?xml version="1.0"?>')
  #cat('<?xml version="1.0"?>', file=file)
  add('<log xmlns="http://www.xes-standard.org/" xes.version="2.0">')
  addExtensions(c(trace_attrs, event_attrs))
  addGlobals(traces, trace_attrs, "trace")
  addGlobals(events, event_attrs, "event")
  addClassifiers()
  addLogAttributes()
  addTraces()
  add('</log>')
  xml <- na.omit(xml)
  writeLines(xml, file)
  #close(fileConn)
}
