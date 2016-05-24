## Helper function for 'setDist'
## Clear row names (age groups)

noGroups <-
function(thisVar){
  for (i in seq(5))
    thisVar[[i, 0]] <- NULL
}