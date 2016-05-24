expn_dd_sp = function(pars1,pars2)
{
    tdmodel = pars2[1]  # type of model
    soc = pars2[2]      # stem (1) or crown (2) age
    res = pars2[3]      # number of points to fit spline to
    ca = pars2[4]       # crown or stem age
    plotit = pars2[5]   # plot the points and the fitted spline (1) or not (0)
    lx = pars2[6]       # number of equations to solve
    cond = pars2[7]     # conditioning on crown/stem age (1) or not (0)
    probs = rep(0,lx)
    probs[1 + (cond == 0) * soc] = 1
    expn = rep(soc,res)
    times = seq(from = -abs(ca),to = 0,length.out = res)
    y = lsoda(probs,times,dd_loglik_rhs,c(pars1[1:min(4,length(pars1))],(cond == 1) * soc,tdmodel - 3),rtol = 1E-10,atol = 1E-16)
    probs = y[1:res,2:(lx + 1)]
    if(cond == 1)
    {
       if(soc == 1) { aux = 1:lx }
       if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
       aux2 = rep(aux,res)
       dim(aux2) = c(lx,res)
       probs = probs/t(aux2)
       sumprobs = probs %*% rep(1,lx)
       sumprobs = rep(sumprobs,lx)
       dim(sumprobs) = c(res,lx)
       probs = probs/sumprobs
    }
    expn = probs %*% (soc:(lx + soc - 1))
    expn2 = probs %*% (soc:(lx + soc - 1))^2
    expn_sp = smooth.spline(times,expn,df = 100)
    expn2_sp = smooth.spline(times,expn2,df = 100)
    if(plotit == 1)
    {
       plot(times,expn,type = 'p')
       lines(expn_sp,lty = 1,col = 2)
    }
    return(list(expn_sp,expn2_sp))
}

predictsp = function(sp,t)
{
   pred = predict(sp,-abs(t))
   return(pred$y)
}

expn_dd = function(pars1,pars2,t)
{
    expnn2_sp = expn_dd_sp(pars1,pars2)
    expn = predictsp(expnn2_sp[[1]],t)
    expn2 = predictsp(expnn2_sp[[2]],t)
    return(list(expn,expn2))
}

lamu_dd_td = function(pars1,pars2,t)
{
   tdmodel = pars2[1]
   ddep = tdmodel - 3
   la = rep(pars1[1],length(t))
   mu = rep(pars1[2],length(t))
   K = pars1[3]
   r = (ddep == 5) * pars1[4]
   n0 = (ddep == 2 | ddep == 4)
   nntd = expn_dd(pars1,pars2,t)[[1]]
   if(ddep == 1)
   {
       latd = pmax(0,la - (la - mu)/K * nntd)
       mutd = mu
   } else {
   if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
   {
       y = -(log(la/mu)/log(K+n0))^(ddep != 2.2)
       latd = pmax(0,la * (nntd + n0)^y)
       mutd = mu
   } else {
   if(ddep == 3)
   {
       latd = la
       mutd = mu + (la - mu) * nntd/K
   } else {
   if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
   {
       y = (log(la/mu)/log(K+n0))^(ddep != 4.2)
       latd = la
       mutd = mu * (nntd + n0)^y
   } else {
   if(ddep == 5)
   {
       latd = pmax(0,la - 1/(r+1)*(la-mu)/K * nntd)
       mutd = mu + r/(r+1)*(la-mu)/K * nntd
   }}}}}
   return(list(la = latd,mu = mutd))
}

#print(lamu_dd_td(pars1 = c(0.8,0.1,40),pars2 = c(4,0,1000,15,1,100), t = c(13.23,2.32) ))