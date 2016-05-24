price.am.option <-
function(k, r, te, w, u.1, u.2, u.3, sigma.1, sigma.2, sigma.3, p.1, p.2)
{
  mu.vector     =  c(u.1,u.2,u.3)
  sigma.vector  =  c(sigma.1, sigma.2, sigma.3)
  prop.vector   =  c(p.1, p.2, 1 - p.1  - p.2)
  expected.f0   =  sum(prop.vector * exp(mu.vector + (sigma.vector^2)/2) )
  prob.f0.gr.k  =  sum(prop.vector * (1 - pnorm((log(k) - mu.vector)/sigma.vector)))
  prob.f0.ls.k  =  1 - prob.f0.gr.k

  expected.f0.f0.gr.k  =  sum(prop.vector*exp((sigma.vector^2 + 2*mu.vector)/2)*(1 - pnorm((log(k) - mu.vector - sigma.vector^2)/sigma.vector)))/prob.f0.gr.k
  expected.f0.f0.ls.k  =  (expected.f0 - expected.f0.f0.gr.k * prob.f0.gr.k)/prob.f0.ls.k

  call.value = (w * (expected.f0.f0.gr.k- k)*prob.f0.gr.k) + 
               (1 - w) * approximate.max(expected.f0 - k,   exp(-r*te) * (expected.f0.f0.gr.k - k)*prob.f0.gr.k )

  call.value = max(0,call.value)  # This is to make sure call never goes negative
 
  put.value  = (w * (k - expected.f0.f0.ls.k)*prob.f0.ls.k) +
               (1 - w) * approximate.max(k - expected.f0,   exp(-r*te) * (k - expected.f0.f0.ls.k)*prob.f0.ls.k )

  put.value  = max(0,put.value)  # This is to make sure put never goes negative
  
  out = list(call.value = call.value, put.value = put.value, expected.f0 = expected.f0, prob.f0.gr.k = prob.f0.gr.k, prob.f0.ls.k = prob.f0.ls.k,
             expected.f0.f0.gr.k = expected.f0.f0.gr.k, expected.f0.f0.ls.k = expected.f0.f0.ls.k)
  out
}
