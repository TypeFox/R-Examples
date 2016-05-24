# this is jsut an illustation of all the parameter passed to the function
# and the request
run <- function(...) {
  ohead("Parameters")
  oprint(str(list(...)))
  ohead("Request")
  oprint(str(request))
  done()
}
