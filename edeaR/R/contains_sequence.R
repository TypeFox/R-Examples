

contains_sequence <- function(sequence, trace, logical = TRUE, interleavings_allowed = TRUE) {
  if(interleavings_allowed == TRUE){# there are mistakes in this function when used for counting!! Must be careful
    sequence_reg_expr <- paste("<",paste(sequence, collapse = ">(.*)<"),">",sep ="")
    trace_string <- paste("<",paste(trace, collapse = ">,<"),">",sep ="")
    
    count = 0
    if(grepl(sequence_reg_expr,trace_string)){
        if(match(sequence[length(sequence)], trace) < length(trace)){
        trace_suffix <- trace[(match(sequence[length(sequence)],trace)+1):length(trace)]
        count = 1 + Recall(sequence,trace_suffix, logical = FALSE, interleavings_allowed = interleavings_allowed)
      } else
        count = 1;
    } 
    
    if(logical == TRUE){
      if(count>0)
        return(TRUE)
      else
        return(FALSE)
    } else
      return(count)
  }
  else{ #no interleaving activities allowed
  	sequence_reg_expr <- paste("<",paste(sequence, collapse = ">,<"),">",sep ="")
  	trace_string <- paste("<",paste(trace, collapse = ">,<"),">",sep ="")

    count = 0
    if(grepl(sequence_reg_expr,trace_string)){
    	start <- regexpr(sequence_reg_expr,trace_string)+nchar(sequence_reg_expr)
      if(start < nchar(trace_string)){
        trace_suffix_string <- substring(trace_string,((start)+1),nchar(trace_string))
        trace_suffix <- strsplit(substring(trace_suffix_string,2, nchar(trace_suffix_string) - 1), split = ">,<")[[1]]
        count = 1 + Recall(sequence,trace_suffix, logical = FALSE, interleavings_allowed = interleavings_allowed)
      } else
        count = 1;
    } 
    
    if(logical == TRUE){
      if(count>0)
        return(TRUE)
      else
        return(FALSE)
    } else
      return(count)
  }
}


