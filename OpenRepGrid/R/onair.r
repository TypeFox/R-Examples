# convert grid file to string


onAirQuery <- function(x, domain="onair.openrepgrid.org")
{
  #if (local)
    #str <- paste0("http://localhost:8100/?grid=", get.query)
  file <- paste0(tempfile(), ".txt")
  saveAsTxt(x, file)
  l <- readLines(file)
  get.query <- paste(l[-(1:3)], collapse="\n") 
  paste0("http://", domain, "/?grid=", get.query)
}

# onAirQuery(boeker)

onAir <- function(x)
{
  query <- onAirQuery(x)
  browseURL(query)  
}
