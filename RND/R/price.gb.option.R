price.gb.option <-
function(r, te, s0, k, y, a, b, v, w)
{
  prob.1 = pgb(k, a = a, b = b, v = (v + 1/a) , w = (w - 1/a) )
  prob.2 = pgb(k, a = a, b = b, v = v ,         w = w )

  cp  = s0 * exp(-y * te) * ( 1 - prob.1) - k * exp(-r*te) * ( 1 - prob.2)
  pp  = cp - s0 * exp(-y * te) + k * exp(-r*te)

  out = list(prob.1 = prob.1, prob.2 = prob.2, call = cp, put = pp)
  out
}
