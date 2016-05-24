"phylpro" <-
    function(input_file, breaks, winHalfWidth, permReps) {
    input_file<-c("Rphylpro", input_file)
    out<-.Call("stepwise", input_file = as.character(input_file), 
      breaks = as.integer(breaks), winHalfWidth = as.integer(winHalfWidth), 
      permReps = as.integer(permReps), PACKAGE="stepwise")    
    if (length(out$winlocs) > 0) class(out)<-"phylpro"
    out
}
