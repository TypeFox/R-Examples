allelic.exact.test <- function(d0,d1,d2,h0,h1,h2) {
  return( .C("newallelic",as.integer(d0),as.integer(d1),as.integer(d2),as.integer(h0),as.integer(h1),
          as.integer(h2),pv=as.double(0))$pv);
}
