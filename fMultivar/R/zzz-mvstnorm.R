
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  dmvsnorm
#  pmvsnorm
#  rmvsnorm
# FUNCTION:             DESCRIPTION:
#  dmvst
#  pmvst
#  rmvst
################################################################################


################################################################################
# Obsolete Functions:
#
#  dmvsnorm(x, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim))
#  pmvsnorm(q, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim))
#  rmvsnorm(n, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim))
#  
#  dmvst(x, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim), df=4)
#  pmvst(q, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim), df=4)
#  rmvst(n, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim), df=4)


# -----------------------------------------------------------------------------


dmvsnorm <- 
  function(x, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim))
    sn::dmsn(x=x, xi=mu, Omega=Omega, alpha=alpha)

pmvsnorm <- 
  function(q, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim))
    sn::pmsn(x=q, xi=mu, Omega=Omega, alpha=alpha)

rmvsnorm <- 
  function(n, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim))
    sn::rmsn(n=n, xi=mu, Omega=Omega, alpha=alpha)


# -----------------------------------------------------------------------------


dmvst <- 
  function(x, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim), df=4)
    sn::dmst(x=x, xi=mu, Omega=Omega, alpha=alpha, nu=df)

pmvst <- 
  function(q, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim), df=4)
    sn::pmst(x=q, xi=mu, Omega=Omega, alpha=alpha, nu=df)

rmvst <- 
  function(n, dim=2, mu=rep(0, dim), Omega=diag(dim), alpha=rep(0, dim), df=4)
    sn::rmst(n=n, xi=mu, Omega=Omega, alpha=alpha, nu=df)


###############################################################################


mvFit <- 
  function(x, method = c("snorm", "st"), fixed.df = NA, 
           title = NULL, description = NULL, trace = FALSE)
{   
  method <- match.arg(method)

  if (method == "snorm") {
    ans <- msnFit(x, trace=trace)
  }
      
  if (method == "st") {
    if (is.na(fixed.df)) fixed.nu <- NULL else fixed.nu <- fixed.df
    ans <- mstFit(x, fixed.nu=fixed.nu, trace=trace)
  }

  # Return Value:
  ans

}


###############################################################################


