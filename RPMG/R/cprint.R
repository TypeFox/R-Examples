`cprint` <-
function(a)
{
##########   dump out an R assignemnt statement to the screen
  nam = deparse(substitute(a))

  if(is.character(a))
    {

      b = as.character(paste(sQuote(a), collapse=","))

    }
  else
    {
      b =  paste(a, collapse=",")
    }

  cat(file="",paste(sep="", nam, "=c(", b , ")") , fill=TRUE)

}

