`gsr` <-
function(Source, Search, Replace, char=FALSE)
{
 #global search and replace
#from R help archive.

  if (length(Search) != length(Replace))
    stop("Search and Replace Must Have Equal Number of Items\n")

  Changed <- as.character(Source)

  if (char == FALSE){
    for (i in 1:length(Search))
    {
      #cat("Replacing: ", Search[i], " With: ", Replace[i], "\n")
      Changed <- replace(Changed, Changed == Search[i], Replace[i])
    }
  }
  else if (char == TRUE){
    for (i in 1:length(Search))
    {
      Changed <- replace(Changed, Changed == Search[i], paste(Replace[i]))
    }
  } 

  Changed
}

