# print.pair prints out INFORMATION ABOUT an object of class "pair"
assign("print.pair",
function(x,...) {

  cat(paste('\nPair object:',deparse(substitute(x)),'\n'))

  if (!is.null(attributes(x)$type))
    cat('\n      Type:            ',attributes(x)$type)
  if (!is.null(attributes(x)$theta))
    cat('\n      Theta:           ',attributes(x)$theta)
  if (!is.null(attributes(x)$dtheta))
    cat('\n      Dtheta:          ',attributes(x)$dtheta)
  cat('\n      Number of pairs: ',length(x$from))
  cat('\n      Number of lags:  ',length(unique(x$lags)))
  cat('\n      Max h:           ',max(x$dist))
  cat('\n\n')


})
