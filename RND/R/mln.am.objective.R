mln.am.objective <-
function(theta, s0, r, te, market.calls, call.weights = NA, market.puts, put.weights = NA, strikes, lambda = 1)
{
   w.1     = theta[1]
   w.2     = theta[2]
   u.1     = theta[3]
   u.2     = theta[4]
   u.3     = theta[5]
   sigma.1 = theta[6]
   sigma.2 = theta[7]
   sigma.3 = theta[8]
   p.1     = theta[9]
   p.2     = theta[10]
   p.3     = 1 - p.1 - p.2

   expected.f0   =  sum(c(p.1, p.2, p.3) * exp(c(u.1,u.2,u.3) + (c(sigma.1, sigma.2, sigma.3)^2)/2) )

   ###
   ### Deal with weights
   ###
 
   if ( sum(is.na(call.weights)) >= 1 )  call.weights = rep(1,length(strikes))
   if ( sum(is.na(put.weights))  >= 1 )  put.weights  = rep(1,length(strikes))

   ###
   ###             call   put
   ###  K < E[f0]   w.1    w.2
   ###  K > E[f0]   w.2    w.1
   ###

   ###
   ###  w.1: Case when k < E[f0] & calls
   ###

   ind.1.calls         = which(strikes < expected.f0)
   len.ind.1.calls     = length(ind.1.calls)
   theoretical.calls.1 = numeric(len.ind.1.calls)

   for (i in 1:len.ind.1.calls)
   {
     theoretical.calls.1[i] = price.am.option(k = strikes[i], r=r, te=te, w = w.1, u.1 = u.1, u.2 = u.2, u.3 = u.3, 
                                                      sigma.1 = sigma.1, sigma.2 = sigma.2, sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$call.value
   }

   obj.1.calls = sum(call.weights[1:len.ind.1.calls] * (theoretical.calls.1 - market.calls[1:len.ind.1.calls])^2)

   ###
   ###  w.2: Case when k < E[f0] & puts 
   ###

   ind.2.puts          = ind.1.calls
   len.ind.2.puts      = len.ind.1.calls
   theoretical.puts.2  = numeric(len.ind.2.puts)

   for (i in 1:len.ind.2.puts)
   {
     theoretical.puts.2[i] = price.am.option(k = strikes[i], r=r, te=te, w = w.2, u.1 = u.1, u.2 = u.2, u.3 = u.3, 
                                                      sigma.1 = sigma.1, sigma.2 = sigma.2, sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$put.value
   }

   obj.2.puts = sum(put.weights[1:len.ind.2.puts] * (theoretical.puts.2 - market.puts[1:len.ind.2.puts])^2)


   ###
   ###  w.2: Case when k > E[f0] & calls
   ###

   ind.2.calls         = which(strikes > expected.f0)
   len.ind.2.calls     = length(ind.2.calls)
   theoretical.calls.2 = numeric(len.ind.2.calls)

   for (i in 1:len.ind.2.calls)
   {
     theoretical.calls.2[i] = price.am.option(k = strikes[i+len.ind.1.calls], r=r, te=te, w = w.2, u.1 = u.1, u.2 = u.2, u.3 = u.3, 
                                                      sigma.1 = sigma.1, sigma.2 = sigma.2, sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$call.value
   }

   obj.2.calls = sum(call.weights[ind.2.calls] * (theoretical.calls.2 - market.calls[ind.2.calls])^2)

   ###
   ###  w.1: Case when k > E[f0] & puts 
   ###

   ind.1.puts          = ind.2.calls
   len.ind.1.puts      = len.ind.2.calls
   theoretical.puts.1  = numeric(length(len.ind.1.puts))

   for (i in 1:len.ind.1.puts)
   {
     theoretical.puts.1[i] = price.am.option(k = strikes[i+len.ind.2.puts], r=r, te=te, w = w.1, u.1 = u.1, u.2 = u.2, u.3 = u.3, 
                                                      sigma.1 = sigma.1, sigma.2 = sigma.2, sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$put.value
   }

   obj.1.puts = sum(put.weights[ind.1.puts] * (theoretical.puts.1 - market.puts[ind.1.puts])^2)

  ###
  ### Get the total objective
  ###
 
  total = obj.1.calls + obj.2.calls + obj.1.puts + obj.2.puts + lambda * (s0 - expected.f0 * exp(-r * te))^2

  ###
  ### Finally ... the objective function
  ###

  if ( (w.1 < 0) | (w.1 > 1) | (w.2 < 0) | (w.2 > 1) | (sigma.1 < 0) | (sigma.2 < 0) | (sigma.3 < 0) |  (p.1 < 0) | (p.1 > 1) | (p.2 < 0) | (p.2 > 1) | ( p.1 + p.2 > 1) ) {obj = 10^7} else { obj = total }
  obj
   
}
