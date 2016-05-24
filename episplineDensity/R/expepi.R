expepi <- function (data, lower = NULL, upper = NULL, N = 10, M = 5, order=2,
          softinfo, opt.args, opt.local.args, postproc.controls = postproc.control()) 
{
#
# Set up for exponential epi-spline fit. Parameters are:
#   lower: numeric, lower end of range for data (default 0)
#   upper: numeric, upper end of range for data (default missing)
#       N: integer, number of segments (so one less than # of mesh points)
#   order: integer, polynomial order (currently only 2 is supported) 
#
if (order != 2) stop ("Only expepi of order 2 is supported for now")

out <- preprocess.data (data, lower, upper)
#
# Set up "epiparameters."
#
epi <- list (m0 = out$m0, mN = out$mN, Ndiscr = N, order = order)
#
# Set up "softinfo". Use the softinfo arguments.
#
soft.list <- list (N = N, M = M, order = order)
if (!missing (softinfo))
    soft.list[names (softinfo)] <- softinfo
soft <- do.call ("setup.softinfo", soft.list)

case <- list (samplesize = length (data))

constraints <- generateLinCon (softinfo = soft, epiparameters = epi, x = out$data)
if (mode(constraints) == "numeric" && constraints[1] == -1)
    stop ("Failure to make constraints")

solve.out <- solveML (data, softinfo = soft, epiparameters = epi, caseinfo = case, constraints = constraints,
             opts = list (maxeval = 10000),
             local.opts = list ("algorithm" = "NLOPT_LD_SLSQP"),
             postproc.controls = postproc.controls)

return (solve.out)
}
