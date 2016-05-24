### TP control.

SPMD.TP <- function(
  bcast = FALSE,
  barrier = TRUE,
  try = TRUE,
  try.silent = FALSE
){
  list(
    bcast = bcast,                   # if bcast object to all ranks at the end
    barrier = barrier,               # if barrier for all ranks at the end
    try = try,                       # if use try in workers
    try.silent = try.silent          # if silent the try message
  )
} # End of SPMD.TP().

