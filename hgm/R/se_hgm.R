##########################################
# se_hgm.R 
# ver. 2014/03/24
# $OpenXM: OpenXM/src/R/r-packages/hgm/R/se_hgm.R,v 1.1 2014/03/25 02:25:26 takayama Exp $
##########################################
#
# Description
#
#  Holonomic gradient method (hgm) for given two points.
#
#
# Usage
#
#  hgm.se.hgm(th0, G0, th1, dG.fun, times=NULL, fn.params=NULL)
#
#
# Arguments
#
#  th0: a vector (d-dim, say); the initial point of th.
#
#  G0: a vector (r-dim, say); the initial value of G.
#
#  th1: a vector (d-dim); the final point of th.
#
#  dG.fun: a function; ``right hand sides'' of the Pfaffian system.
#         The inputs are d-dim and r-dim vectors. The output is a d*r-dim array.
#
#  times: a vector; times in [0,1] at which explicit estimates for G are desired.
#         If time = NULL, the set {0,1} is used, and only the final value is returned.
#
#  fn.params: a list of parameters passed to the function dG.fun.
#         If fn.params = NULL, no parameter is passed to dG.fun.
#
#
# Details
#
#  The function hgm.se.hgm computes the value of a holonomic function
#  at a given point, using HGM.
#  This is a ``Step 3'' function in the HGM framework.
#
#  The Pfaffian system assumed is
#    d G_j / d th_i = (dG.fun(th, G))_{i,j}
#
#  The inputs of hgm.se.hgm are the initial point th0, initial value G0, final point th1,
#  and Pfaffian system dG.fun. The output is the final value G1.
#
#  If the argument `times' is specified, the function returns a matrix,
#  where the first column denotes time, the following d-vector denotes th,
#  and the remaining r-vector denotes G.
#
#
# Usage
#
#  See se_demo.R
#
##########################################


hgm.se.hgm = function(th0, G0, th1, dG.fun, times=NULL, fn.params=NULL){
  if(is.null(times)){
    show.trace = FALSE
    times = c(0,1)
  }else{
    show.trace = TRUE
  }
  params = list(th = th0, dG.fun = dG.fun, v = th1-th0, fn.params = fn.params)
  my.ode = function(tau, G, params){
    th = params$th
    dG.fun = params$dG.fun
    v = params$v
    fn.params = params$fn.params
    th = th + tau * v
    if(is.null(fn.params))  dG = dG.fun(th, G)
    else dG = dG.fun(th, G, fn.params)
    G.rhs = v %*% dG
    list(G.rhs)
  }
  rk.res = rk(G0, times, my.ode, params)
  G1 = as.vector(rk.res[nrow(rk.res), 1+(1:length(G0))])
  if(show.trace) return(rk.res)
  else return(G1)
}
# EOF
