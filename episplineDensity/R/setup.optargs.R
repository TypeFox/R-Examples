setup.optargs <- function (param.length, opts, local.opts) 
{
#
# Set up optimizer global and local arguments.
#
# param.length: integer, length of parameter vector
#         opts: list of named arguments to global optimizer
#   local.opts: list of named arguments to local optimizer
#
default.opts <- list (algorithm = "NLOPT_LD_AUGLAG", print_level = 0, 
                      maxeval = 2500, xtol_rel = 1e-5, 
                      xtol_abs = rep (1e-5, param.length))
#
# If opts is supplied, use its entries in place of the defaults.
#
if (!missing (opts))
    default.opts[names (opts)] <- opts
opts <- default.opts
#
# If local.opts is supplied, use its entries in place of defaults.
#
default.local.opts <- list (algorithm = "NLOPT_LD_SLSQP", 
                      maxeval = 1000, xtol_rel = 1e-5)
if (!missing (local.opts))
    default.local.opts[names (local.opts)] <- local.opts
#
# Local opts get attached to the globals; that's how nloptr() wants it.
#
opts$local_opts <- default.local.opts

return (opts)
}
