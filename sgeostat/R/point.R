assign("point",
function (dframe,x='x',y='y') {

  if (!is.data.frame(dframe))
    stop("dframe must be a data frame.")

  if (is.na(match(x,names(dframe))))
    stop ("Could not find the X column.")

  if (is.na(match(y,names(dframe))))
    stop ("Could not find the Y column.")

  names(dframe)[match(x,names(dframe))] <- 'x'
  names(dframe)[match(y,names(dframe))] <- 'y'

  o.point <- data.frame(dframe)

  class(o.point) <- c("point","data.frame")
  return(o.point)
})
