summary.tscgm <- function(object, ...)
{
  precision.estimate <- object$theta
	autoReg.estimate <- object$gamma
	optimal.tuning.prec <-  object$lam1.opt
	optimal.tuning.autoReg <- object$lam2.opt
	min.value.model.selection <- object$min.ic
	model.selection <- object$tun.ic
	precision.sparsity <- object$s.theta
	autoReg.sparsity <-  object$s.gamma
	
	output = list(precision.estimate = precision.estimate,
    autoReg.estimate = autoReg.estimate,
    optimal.tuning.prec = optimal.tuning.prec,
    optimal.tuning.autoReg = optimal.tuning.autoReg,
	  min.value.model.selection  = min.value.model.selection, 
    model.selection = model.selection, 
    precision.sparsity = precision.sparsity,
    	autoReg.sparsity = 	autoReg.sparsity
  )
  class(output) <- "summary.tscgm"
  output
	
}

    