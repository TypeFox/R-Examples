untangle.specials <- function (tt, special, order = 1) 
{
  spc <- attr(tt, "specials")[[special]]
  if (length(spc) == 0) 
    return(list(vars = character(0), terms = numeric(0)))
  facs <- attr(tt, "factors")
  fname <- dimnames(facs)
  ff <- apply(facs[spc, , drop = FALSE], 2, sum)
  list(vars = (fname[[1]])[spc], terms = seq(ff)[ff & match(attr(tt, 
                                                                 "order"), order, nomatch = 0)])
}