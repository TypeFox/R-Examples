price.mln.option <-
function(r, te, y, k, alpha.1, meanlog.1, meanlog.2, sdlog.1, sdlog.2)
{

  discount.factor =  exp(-r * te)
  alpha.2         =  1 - alpha.1
  expected.val.1  =  exp(meanlog.1 + (0.5)*(sdlog.1^2))
  expected.val.2  =  exp(meanlog.2 + (0.5)*(sdlog.2^2))
  s0              =  exp((y-r) * te) * (alpha.1 * expected.val.1 + alpha.2 * expected.val.2)

  u1   = (log(k) - meanlog.1)/sdlog.1
  tmp1 = expected.val.1 * (1 - pnorm(u1 - sdlog.1)) - k * (1 - pnorm(u1))
  c1   = discount.factor * tmp1

  u2   = (log(k) - meanlog.2)/sdlog.2
  tmp2 = expected.val.2 * (1 - pnorm(u2 - sdlog.2)) - k * (1 - pnorm(u2))
  c2   = discount.factor * tmp2

  call.option.price = alpha.1 * c1 + alpha.2 * c2
  put.option.price  = call.option.price - s0 * exp(-y*te) + k * discount.factor
  out = list(call = call.option.price, put = put.option.price, s0 = s0)
  out

}
