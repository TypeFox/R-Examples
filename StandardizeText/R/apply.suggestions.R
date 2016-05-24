apply.suggestions <-
function(suggest,suggested.values,initial.values,standardized.values,na.rm,print.changes){
  suggestion.count <- dim(suggested.values)[1]
  unapplied.suggestions <- 1:suggestion.count
  if(suggest=="auto" || suggest=="a") {
    for(i in 1:suggestion.count) {
      standardized.values[which(initial.values==suggested.values[i,1])] <- as.character(suggested.values[i,2])
    }
    if(print.changes) {
      cat("\nThe following suggested changes were applied:\n")
      print(suggested.values)
    } else {
      cat(sprintf("\n%d suggested changes were applied.\n",suggestion.count))
    }
  } else {
    cat("\nThe following suggestions were identified:\n")
    print(suggested.values)
    cat("\nWhich suggestions do you want to apply: (a)ll, (n)one, or supply line numbers [e.g. 2, 4:7, 9]\n")
    suggest.action <- readline(">>> ")
    if(suggest.action == "a" || suggest.action == "all") {
      for(i in 1:suggestion.count) {
        standardized.values[which(initial.values==suggested.values[i,1])] <- as.character(suggested.values[i,2])
      }
      unapplied.suggestions <- integer(0)
      cat("\nAll suggestions applied.\n")
    }else if(suggest.action == "n" || suggest.action == "none" || suggest.action == "") {
      cat("\nNo suggestions applied.\n")
    }else {
      tryCatch({
        suggest.action <- gsub("([^,:0-9]+)","",suggest.action)
        suggest.action <- paste("c(",suggest.action,")",sep="")
        selected.lines <- eval(parse(text=suggest.action))
        if(is.vector(selected.lines) && is.numeric(selected.lines) && sum(as.integer(selected.lines)!=selected.lines)==0
           && length(unique(selected.lines))==length(selected.lines) && min(selected.lines)>=1 && max(selected.lines)<=suggestion.count) {
          for(i in selected.lines) {
            standardized.values[which(initial.values==suggested.values[i,1])] <- as.character(suggested.values[i,2])
          }
          unapplied.suggestions <- setdiff(unapplied.suggestions,selected.lines)
          if(print.changes) {
            cat("\nThe following suggestions were applied:\n")
            print(suggested.values[selected.lines,])
          } else {
            cat(sprintf("\nThe selected %d suggestions were applied.\n",length(selected.lines)))
          }
        } else {
          cat("\nNo suggestions applied. Please enter in-range integer values.\n")
        }
      }, warning = function(w) {
        cat("\nUnrecognized input. No suggestions applied.\n")
      }, error = function(e) {
        cat("\nUnrecognized input. No suggestions applied.\n")
      })
    }
    if(na.rm && length(unapplied.suggestions)>0) {
      for(i in unapplied.suggestions) {
        standardized.values[which(initial.values==suggested.values[i,1])] <- NA
      }
      if(print.changes) {
        cat("\nThe following names ignored suggestions and were removed:\n")
        print(suggested.values[unapplied.suggestions,1])
      } else {
        cat(sprintf("\nThe nonselected %d names ignored suggestions and were removed.\n",length(unapplied.suggestions)))
      }
    }
  }
  return(standardized.values)
}
