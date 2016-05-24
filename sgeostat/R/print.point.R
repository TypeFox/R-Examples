# print.point prints out INFORMATION ABOUT an object of class "point"
assign("print.point",
function(x,...) {

  cat(paste('\nPoint object:',deparse(substitute(x)),'\n'))

  cat(paste('\n   Locations: ',length(x$x),sep=''))
  cat(paste('\n\n   Attributes:\n      ',paste(names(x),
                    collapse='\n      '),sep=''))
  cat('\n\n')


})
