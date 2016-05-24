#' Calculates the simulated exposure profile (EE, NEE, PFE, EEE) by use of the Hull-White model. Two sets of results are provided:
#' one after taking into account the marging agreement and one assuming that there is no marging agreement present
#' @title Calculated the Simulated Exposure Profile
#' @param discount_factors The discount curve derived from the spot curve
#' @param spot_curve The curve derived from interpolating the market spot rates
#' @param col The margin agreement
#' @param trades The list of the trade objects
#' @param time_points The timepoints that the analysis is performed on
#' @param sim_data A list containing simulation-related data (model parameters and number of simulation)
#' @return A list containing the exposure profile (both collateralized and uncollateralized)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'

CalcSimulatedExposure = function(discount_factors, time_points, spot_curve, col, trades, sim_data)
{
  num_of_points = length(time_points)
  num_of_trades = length(trades)
  Swap_MtMs         = matrix(0,nrow=sim_data$num_of_sims, ncol = num_of_points)
  Swap_MtMs_coll    = matrix(0,nrow=sim_data$num_of_sims, ncol = num_of_points)

  timesteps_diff = diff(time_points)

  spot_interest_rate = spot_curve[1]
  forward_curve = (discount_factors[1:(num_of_points-1)]/discount_factors[2:num_of_points]-1)/timesteps_diff
  forward_curve = c(spot_interest_rate,forward_curve)
  forward_diff=diff(forward_curve)
  theta=forward_diff+sim_data$mean_reversion_a*forward_curve[2:length(forward_curve)]
  set.seed(30269)
  interest_rates    = rep(0,length(num_of_points))
  interest_rates[1] = spot_interest_rate
  random_numbers     = matrix(runif(sim_data$num_of_sims*(num_of_points-1)),nrow=sim_data$num_of_sims,ncol=num_of_points-1)

for(index in 1:num_of_trades)
{
maturity  = trades[[index]]$Ei
swap_rate = trades[[index]]$swap_rate
BuySell   = ifelse(trades[[index]]$BuySell=='Buy',1,-1)

  time_points_temp   = time_points[time_points<=maturity]
  num_of_points_temp = length(time_points_temp)
  A = rep(0,length(time_points_temp))
  B = (1-exp(-sim_data$mean_reversion_a*(maturity-time_points_temp)))/sim_data$mean_reversion_a

  disc_factors = matrix(0,nrow=num_of_points_temp,ncol=num_of_points_temp)

  dt = maturity/num_of_points_temp

  for(j in 1:sim_data$num_of_sims)
  {
    Floating_leg = rep(0,num_of_points_temp)
    for(i in 2:num_of_points_temp)
      interest_rates[i] = interest_rates[i-1] + (theta[i-1]-sim_data$mean_reversion_a*interest_rates[i-1])*timesteps_diff[i-1]+ sim_data$volatility*qnorm(random_numbers[j,i-1])*sqrt(timesteps_diff[i-1])

    for (i in 1:num_of_points_temp)
      A[i] = discount_factors[num_of_points_temp]/discount_factors[i]*exp(B[i]*forward_curve[i]-(sim_data$volatility^2/(4*sim_data$mean_reversion_a))*(1-exp(-2*sim_data$mean_reversion_a*time_points_temp[i]))*B[i]^2)

    for (i in 1:num_of_points_temp)
    disc_factors[i,1:(num_of_points_temp-i+1)] = A[num_of_points_temp-i+1]*exp(-B[num_of_points_temp-i+1]*interest_rates[1:(num_of_points_temp-i+1)])


    for (i in 0:(num_of_points_temp-1))
      Floating_leg[i+1] = 1-disc_factors[num_of_points_temp-i,i+1]

    Fixed_leg = dt*swap_rate*colSums(disc_factors,na.rm = TRUE)

    Swap_MtM = (Floating_leg - Fixed_leg)*BuySell

    Swap_MtMs[j,1:num_of_points_temp]      = Swap_MtMs[j,1:num_of_points_temp]      + Swap_MtM
    Swap_MtMs_coll[j,1:num_of_points_temp] = Swap_MtMs_coll[j,1:num_of_points_temp] + col$ApplyThres(Swap_MtM)
  }
  trades[[index]]$MtM = Swap_MtMs[1,1]
}
  exposure_profile = list()

  exposure_profile$EE_uncoll  = apply(apply(Swap_MtMs,c(1,2),function(x) ifelse(x>=0,x,NA)),2,mean,na.rm=TRUE)
  exposure_profile$NEE_uncoll = apply(apply(Swap_MtMs,c(1,2),function(x) ifelse(x<0,x,NA)),2,mean,na.rm=TRUE)
  exposure_profile$EE_uncoll[is.na(exposure_profile$EE_uncoll)] = 0
  exposure_profile$NEE_uncoll[is.na(exposure_profile$NEE_uncoll)] = 0
  exposure_profile$PFE_uncoll = apply(apply(Swap_MtMs,c(1,2),function(x) ifelse(x>0,x,NA)),2,quantile,sim_data$PFE_Percentile,na.rm=TRUE)
  exposure_profile$PFE_uncoll[is.na(exposure_profile$PFE_uncoll)] = 0

  exposure_profile$EE  = apply(apply(Swap_MtMs_coll,c(1,2),function(x) ifelse(x>=0,x,NA)),2,mean,na.rm=TRUE)
  exposure_profile$EE[is.na(exposure_profile$EE)] = 0
  exposure_profile$NEE = apply(apply(Swap_MtMs_coll,c(1,2),function(x) ifelse(x<0,x,NA)),2,mean,na.rm=TRUE)
  exposure_profile$NEE[is.na(exposure_profile$NEE)] = 0
  exposure_profile$PFE = apply(apply(Swap_MtMs_coll,c(1,2),function(x) ifelse(x>0,x,NA)),2,quantile,sim_data$PFE_Percentile,na.rm=TRUE)
  exposure_profile$PFE[is.na(exposure_profile$PFE)] = 0

  EEE = rep(0,num_of_points)

  EEE[1] = exposure_profile$EE[1]

  for(i in 2:(num_of_points))
    EEE[i] = max(EEE[i-1],exposure_profile$EE[i])

  exposure_profile$EEE = EEE

  return(exposure_profile)
}
