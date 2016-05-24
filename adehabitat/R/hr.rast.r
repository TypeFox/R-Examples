"hr.rast" <- function(mcp, w, border=c("include", "exclude"))
{
    ## Verifications
    if (inherits(w, "asc"))
      w <- as.kasc(list(to=w))
    if (!inherits(w, "kasc"))
      stop("Non convenient data")
    if (!inherits(mcp, "area"))
      stop("mcp should be of class \"area\"")
    bord <- match.arg(border)

    ## a list with one element = one polygon
    lpc<-split(mcp[,2:3], mcp[,1])
    output<-list()

    ## use of the function mcp.rast for each polygon
    for (i in 1:length(lpc))
      output[[names(lpc)[i]]]<-mcp.rast(lpc[[i]], w, bord)

    ## the output:
    output<-as.kasc(output)
    return(output)
  }

