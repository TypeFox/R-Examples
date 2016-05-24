################################################################################
# routines for checking water balances
#-------------------------------------------------------------------------------
#
################################################################################

# Get the current storage excluding naything in the channel
# this can be used initialise water balances
current.storage <- function(groups, stores, ichan=NULL)
{
  if(!is.null(ichan))
  {
    stores <- stores[-ichan,]
    groups <- groups[-ichan,]

  }

  # storage remaining in root and unsaturated zones plus saturated storage in water table
  # subtracting the sd gives the effective overall storage across all zone. -ve values indicate deficit

  return(weighted.average(stores$srz+stores$ex+stores$suz-stores$sd, groups$area))
}

# net input at time step
current.input <- function(groups, rain, ae, qr)
{
  return(as.numeric(rain - weighted.average(ae, groups$area) - qr))
}

# check water balance from output of Dynamic TOPMODEL run
water.balance <- function(groups, stores, dt, storage.in, ichan,
                          qsim,
                          rain, ae=0)
{
  eff.rain <- rain-ae
  eff.rain <- apply(eff.rain[1:nrow(qsim),], MARGIN=1, function(x)sum(x*groups$area))/sum(groups$area)

  # rain distribution
  #areas <- matrix(rep(groups$area, nrow(rain)), ncol=nrow(groups), byrow=T)
#  tot.rain <- rowSums(()*areas)/sum(groups$area)
#  tot.rain<-tot.rain[1:nrow(qsim)]
  # difference between net input and net storage gain
  storage.gain <- current.storage(groups, stores, ichan)-storage.in
  wb <- dt*sum(eff.rain - qsim[,1], na.rm=T)-storage.gain
  return(wb)
}



