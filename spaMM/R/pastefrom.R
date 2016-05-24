# From http://stackoverflow.com/questions/7307987/logging-current-function-name

pastefrom <- function(..., callstack=sys.calls(),prefix="From ",full.stack=.spaMM.data$options$MESSAGES.FULL.STACK) {  ## callstack is a list of *calls*
  if (full.stack) {
    cs <- clean_cs(callstack) ## cs is a string made of the names of function called in the call stack
  } else {
    lcs <- length(callstack)
    cs <- paste(callstack[[lcs-1]])[1] ## [[lcs-1]] because [[lcs]] is call to pastefrom.
  }
  paste(prefix,cs,": ", ...,sep="")
}

clean_cs <- function(callstack){ ## callstack is a list of *calls*
  val <- sapply(callstack, function(xt){ ## for each call...
    z <- strsplit(paste(xt, collapse="\t"), "\t")[[1]] ## the function name;
    switch(z[1], ## ...inspects the fn name...
        "lapply" = z[3], ## ... for lapply call, get the name of the called function...
        "sapply" = z[3], ## ... etc...
        "do.call" = z[2], 
        "function" = "FUN",
        "source" = "###",
        "eval.with.vis" = "###",
        z[1]
        )
    }) ## val is a character vector
  val[grepl("\\<function\\>", val)] <- "FUN" ## replaces  any pure fn definition [as called by lapply(,function...)] by "FUN"
  val <- val[!grepl("(###|FUN)", val)] ## then removes "###" and "FUN" from the vector
  val <- head(val, -1) ## removes the last element which is the name of the function calling clean_cs, ie typically 'pastefrom'
  paste(val, collapse="|")
}
