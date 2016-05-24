plot.autofitVariogram = function(x, plotit = TRUE, ...)
{ 
    shift = 0.03
    labels = as.character(x$exp_var$np)
    vario = xyplot(gamma ~ dist, data = x$exp_var, panel = autokrige.vgm.panel,
                labels = labels, shift = shift, model = x$var_model,# subscripts = TRUE,
                direction = c(x$exp_var$dir.hor[1], x$exp_var$dir.ver[1]),
                ylim = c(min(0, 1.04 * min(x$exp_var$gamma)), 1.04 * max(x$exp_var$gamma)),
                xlim = c(0, 1.04 * max(x$exp_var$dist)), xlab = "Distance", ylab = "Semi-variance",
                main = "Experimental variogram and fitted variogram model", mode = "direct",...)
   if (plotit) print(vario) else vario
}