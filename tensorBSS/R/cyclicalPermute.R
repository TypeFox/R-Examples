cyclicalPermute <-
function(m, r){
  if(m == 1) part1 <- NULL
  else part1 <- 1:(m-1)
  part2 <- (m:(r+1))[c((r+1-m+1), (1:(r+1-m)))]
  c(part1, part2)
}
