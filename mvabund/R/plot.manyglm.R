################################################################################
# PLOT.MANYLM, PLOT.MANYGLM:                                                   #
# Plot for evaluation of goodness of fit for lm.mvabund objects                #
################################################################################

plot.manyglm <- function(x, res.type="pit.norm", which = 1, caption = c("Residuals vs Fitted","Normal Q-Q", "Scale-Location", "Cook\'s distance"), overlay=TRUE, n.vars=Inf, var.subset=NULL, sub.caption = NULL, ... ) 
{	
   default.plot.manyglm(x, res.type=res.type, which = which, caption = caption, overlay=overlay, n.vars=n.vars, var.subset=var.subset, sub.caption = sub.caption, ... )
   return(invisible())
}
 
