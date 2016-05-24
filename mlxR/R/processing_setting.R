processing_setting <- function(s1)
{
  s2=list()
  data.in  <- FALSE
  id.out   <- FALSE
  if (!is.null(s1)){
    if (!is.null(s1$record.file)) {s2$recordFile=s1$record.file}
    if (!is.null(s1$load.design)) {s2$loadDesign=s1$load.design}
    if (!is.null(s1$seed))        {s2$seed=s1$seed}
    if (!is.null(s1$data.in))     {data.in=s1$data.in}
    if (!is.null(s1$id.out))      {id.out=s1$id.out}
    
  }
  return(list(s2, data.in, id.out))
}



