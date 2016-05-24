dd_loglik_bw_rhs = function(t,x,pars)
{
lx = length(x) - 1
la = pars[1]
mu = pars[2]
K = pars[3]
if(length(pars) < 6)
{
   kk = pars[4]
   ddep = pars[5]
} else {
   r = pars[4]
   kk = pars[5]
   ddep = pars[6]
}
n0 = (ddep == 2 | ddep == 4)

nn = -1:(lx+2*kk)
lnn = length(nn)
nn = pmax(rep(0,lnn),nn)

if(ddep == 1)
{
    lavec = pmax(rep(0,lnn),la - (la-mu)/K * nn)
    muvec = mu * rep(1,lnn)
} else {
if(ddep == 1.3)
{
    lavec = pmax(rep(0,lnn),la * (1 - nn/K))
    muvec = mu * rep(1,lnn)
} else {
if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
{
    y = -(log(la/mu)/log(K+n0))^(ddep != 2.2)
    lavec = pmax(rep(0,lnn),la * (nn + n0)^y)
    lavec[which((nn + n0)^y == Inf)] = 0
    muvec = mu * rep(1,lnn)
} else {
if(ddep == 2.3)
{
    y = K
    lavec = pmax(rep(0,lnn),la * (nn + n0)^y)
    muvec = mu * rep(1,lnn)
} else {
if(ddep == 3)
{
    lavec = la * rep(1,lnn)
    muvec = mu + (la - mu) * nn/K
} else {
if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
{
    lavec = la * rep(1,lnn)
    y = (log(la/mu)/log(K+n0))^(ddep != 4.2)
    muvec = mu * (nn + n0)^y
    muvec[which((nn + n0)^y == Inf)] = 0
} else {
if(ddep == 5)
{ 
    lavec = pmax(rep(0,lnn),la - 1/(r+1)*(la-mu)/K * nn)
    muvec = muvec = mu + r/(r+1)*(la-mu)/K * nn
}}}}}}}
xx = c(0,x[1:lx],0)

dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
#dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]

dG = x[1 + (kk == 0)]

return(list(c(dx,dG)))
}