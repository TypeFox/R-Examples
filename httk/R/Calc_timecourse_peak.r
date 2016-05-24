# This function calculates the peak concentration given a PK timecourse (concentration vs. time):
calc_timecourse_peak <- function(timecourse,conc.col=2,time.col=1)
{
  return(max(timecourse[,conc.col]))
}
