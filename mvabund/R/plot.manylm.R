################################################################################
# PLOT.MANYLM, PLOT.MANYGLM:                                                   #
# Plot for evaluation of goodness of fit for lm.mvabund objects                #
################################################################################

plot.manylm  <- function(x, res.type="pearson", which = 1:4, caption = c("Residuals vs Fitted","Normal Q-Q", "Scale-Location", "Cook\'s distance"), overlay=TRUE, n.vars=Inf, var.subset=NULL, sub.caption = NULL, studentized=TRUE,	... ) 
{	
   dots <- list(...)
   if ("res.type" %in% names(dots))
   {
     if ( res.type!="pearson" & res.type!="pit.norm" ) 
       warning("Only Pearson residuals implemented for manylm functions.")
   }
   default.plot.manylm(x, which = which, caption = caption, overlay=overlay, n.vars=n.vars, var.subset=var.subset, sub.caption = sub.caption, studentized=studentized, ... )
   
   return(invisible())
}
