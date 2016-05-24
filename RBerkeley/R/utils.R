check_pointer <- function(dbh) {
  .Call("rberkeley_check_pointer", dbh)
}

as.DB <- function(x) {
  if(is.DB(x)) {
    return(attr(x, "conn_id"))
  } else stop("not a DB connection")
}

is.DB <- function(x) {
  inherits(x, "DB") 
}

print.DB <- function(x, ...) {
  cat("<DB Handle>\n")
}
