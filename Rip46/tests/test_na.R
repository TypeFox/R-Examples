local({

  ## TODO: port to testthat
require(Rip46)
require(tools)

f <- function(x) stopifnot(all(is.na(  x  )))

stopifnot( as.character(as.ip4("192.168.1.1")) == "192.168.1.1")

f( as.ip4(as.character(NA)) )
f( as.ip4(as.numeric(NA)) )
f( as.ip4(as.integer(NA)) )
f( as.ip4("NA") )
f( as.ip4("300.300.300.300") )
f( as.character(as.ip4("NA")) )

f( as.ip6("QQ:aa:11::") )
f( as.ip6(as.character(NA))  )
f( as.character(as.ip6(as.character(NA)))  )
})