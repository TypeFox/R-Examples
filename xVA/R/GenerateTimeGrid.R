GenerateTimeGrid = function(col, maturity)
{
remargin_freq_annum = col$remargin_freq/360
mpor_days_annum = col$mpor_days/360
time_points = seq(0,maturity,remargin_freq_annum)
time_points_lookback = seq(remargin_freq_annum-mpor_days_annum,maturity-mpor_days_annum,remargin_freq_annum)

time_points = sort(c(time_points_lookback,time_points))

return(time_points)

}