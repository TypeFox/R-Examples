textfun <-
function(text, speed = 0.05)
{
  text <- strsplit(text,""," ")[[1L]]
  Sys.sleep(speed)
  for(w in text) {
    cat(w)
    Sys.sleep(speed)
  }
  Sys.sleep(speed)
  cat("\n")

  return(invisible(NULL))
}

