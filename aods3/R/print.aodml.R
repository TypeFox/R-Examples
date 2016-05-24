print.aodml <- function(x, ...)
  
  print(
		list(
			call = x$call,
  		b = x$b, phi = x$phi, phi.scale = x$phi.scale,
  		#df.model = df.model, df.residual = df.residual,
  		varparam = x$varparam, logL = x$logL,
  		iterations = x$iterations, code = x$code
			)
    )
