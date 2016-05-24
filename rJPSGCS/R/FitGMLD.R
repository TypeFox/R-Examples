FitGMLD <-
function(par="input.par",fped="flipped.ped", out.ld.par="out.ld.par")
{
  # The function will create file out.ld.par in the current working directory
  .jcall("FitGMLD","V","main", .jarray(c(par, fped, out.ld.par)))
}

