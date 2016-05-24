markup <- function(x)
  gsub("<", "&lt;", gsub(">", "&gt;", gsub("&", "&amp;", x)))

catSDML <- function(...) cat(..., append=TRUE, sep="")

etagsSDML <- function(x, ...)
{
  info <- attr(x, "info")
  for (i in 1:length(x)) {
    inf <- if (!is.null(info) && !is.na(info[i]))
      paste(" info=\"", info[i], "\"", sep="")
    else ""
    
    if (is.na(x[i]) && !is.nan(x[i]))
      catSDML("<na", markup(inf), "/>", ...)
    else
      if (is.logical(x))
        catSDML("<", if (x[i]) "T" else "F", markup(inf), "/>", ...)
      else
        tags(x[i], "e", inf, ...)
    
    if((i!=length(x)) & ((i %% 5)==0))
      catSDML("\n", ...)
  }
}

cetagsSDML <- function(x, ...)
{
  info <- attr(x, "info")
  for (i in 1:length(x)) {
    inf <- if (!is.null(info) && !is.na(info[i]))
      paste(" info=\"", info[i], "\"", sep="")
    else ""
    
    if (is.na(x[i]) && !is.nan(x[i]))
      catSDML("<na", markup(inf), "/>", ...)
    else {
      catSDML("<ce", markup(inf), ">", ...)
      tags(Re(x[i]), "r", ...)
      tags(Im(x[i]), "i", ...)
      catSDML("</ce>", ...)
    }
    
    if((i!=length(x)) & ((i %% 3)==0))
      catSDML("\n", ...)
  }
}

tags <- function(x, s, info = "", ...) {
  if (is.nan(x))
    catSDML("<", s, markup(info), "><nan/></", s, ">", ...)
  else if (x == Inf)
    catSDML("<", s, markup(info), "><posinf/></", s, ">", ...)
  else if (x == -Inf)
    catSDML("<", s, markup(info), "><neginf/></", s, ">", ...)
  else
    catSDML("<", s, markup(info), ">", markup(x), "</", s, ">", ...)
}

getAttrSDML <- function(x)
{
	if (!is.null(x$attributes))
	{
		return(x$attributes)
	}
	return(NULL)
}
