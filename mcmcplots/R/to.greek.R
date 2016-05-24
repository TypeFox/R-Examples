.to.greek <- function(instr) {
  alphabet <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", 
    "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", 
    "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", 
    "psi", "omega")
  has.greek <- sapply(instr, function(l) 
      any(sapply(alphabet, function(x) length(grep(x, l, ignore.case=TRUE)) > 0)))
  instr[has.greek] <- parse(text=gsub(",", "*','*", instr[has.greek]))
  instr
}