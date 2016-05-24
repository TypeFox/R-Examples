s4hm <-
function(dir, model.name)
{
  if(length(bm <- search.bayesx.models(dir)))
    bm <- unique(grep(paste(model.name, "_hlevel", sep = ""), bm, fixed = TRUE, value = TRUE))
  else
    stop(paste("problems reading model results from directory:", dir))

  return(bm)
}

