correctFormat <-
function(arg, square.data) {
  if (!is.null(arg)){ #should be a max length of 1
    if (length(arg) <= 1) {
      if (is.character(arg)) { 
        arg_temp <- type.convert(arg, as.is=T)
        arg_temp <- which(colnames(square.data) == arg_temp) ; return(arg_temp)
      }
      else if (is.numeric(arg)) { 
        return(arg)
      }    
      else {
        sprintf('Input for option %s could not be treated as a character or vector', arg)
      } 
    }
    else if (length(arg > 1)){
      argList <- rep(0,length(arg))
      if (is.numeric(arg)) { 
        argList <- arg ; return(unique(argList))
      } 
      else {
        for (each in 1:length(arg)) {
          arg_temp <- type.convert(arg[each], as.is=T)
          if (is.numeric(arg_temp)) {
            argList[each] <- as.numeric(arg_temp)
          } else if (is.character(arg_temp)) {
            argList[each] <- as.numeric(which(colnames(square.data) == arg_temp))
          } else print('errors')
        } 
        return(unique(argList))
      }
    } 
  }
}
