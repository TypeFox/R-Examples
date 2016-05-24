PZ2Resp = function(PZ, f, PLOT = TRUE){
  if(length(PZ$Sense) == 0){
    PZ$Sense = 1
  }
  if(length(PZ$Knorm) == 0){
    PZ$Knorm = 1
  }
  
  # take poles/zeros structure, such as element of PreSet.Instr output, and generate transfer function
  s = 2i * pi * f
  zeros = PZ$zeros
  poles = PZ$poles
  A0 = PZ$Knorm * PZ$Sense
  resp = 1
  if(is.null(A0)){
    A0 = 1/max(abs(resp))
  }
  for(i in seq(1, length(zeros), length.out = length(zeros))){
    resp = resp * (s - zeros[i])
  }
  for(i in seq(1, length(poles), length.out = length(poles))){
    resp = resp / (s - poles[i])
  }
  resp = resp * A0

  if(PLOT){
      w = which(f > 0)
    plot(f[w], abs(resp)[w], type = 'l', ylim = max(abs(resp[w])) * c(0.01, 1), log = 'xy')
  }
  invisible(resp)
}
