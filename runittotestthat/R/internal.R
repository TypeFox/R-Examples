fixup <- function(y)
{
  if("info" %in% names(yl <- as.list(y)) && is.null(yl$info)) 
  { 
    yl$info <- NULL
    y <- as.call(yl)
  }
  y
}

get_yes_no_response <- function(...)
{  
  prompt <- if(nargs() == 0) 
  {
    "Please respond (y)es or (n)o > " 
  } else 
  {
    paste0(...)
  }
  repeat
  {
    typed <- tolower(readline(prompt))
    if(typed %in% c("y", "ye", "yes"))
    {
      return(TRUE)      
    }
    if(typed %in% c("n", "no"))
    {
      return(FALSE)    
    }
    message("Input unknown.  Please try again.")
  } 
}