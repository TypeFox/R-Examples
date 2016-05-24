`dlcross` <- function(format=c("rqtl", "dlmap", "other"), genobj, pheobj, mapobj, idname="ID", genfile, phefile, mapfile, type, step=0, fixpos=0, estmap=TRUE, ...)
{
  results <- list()

  if (missing(format)) stop("format is a required argument and is missing") 

  if (format == "rqtl")
    results <- dlcross.cross(genobj, pheobj, idname, step, fixpos, estmap)

  if (format=="dlmap")
	results <- dlcross.dlmap(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile, type, step, fixpos, estmap)

  if (format=="other")
	results <- dlcross.other(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile)

  results
}

