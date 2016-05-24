### IO control.

SPMD.IO <- function(
  max.read.size = 5.2e6,
  max.test.lines = 500,
  read.method = c("gbd", "common"),
  balance.method = c("block", "block0", "block.cyclic")
){
  list(
    max.read.size = max.read.size,                         # 5 MB
    max.test.lines = max.test.lines,                       # test lines
    read.method = read.method[1],                          # read methods
    balance.method = balance.method[1]                     # for gbd only
  )
} # End of SPMD.IO().

