`plot.mlcm` <-
function(x, standard.scale = FALSE, transpose = FALSE, SD.scale = FALSE, ...){
   par(ask = FALSE)
   m <- if (transpose) t(x$pscale) else x$pscale
   if (standard.scale){
   	dm <- dim(m)
   	mx <- max(m[dm[1], ])
   	matplot(m/mx, ...)
   	} else
   	matplot((if (SD.scale) 2 else 1) * m, ...)
}

`lines.mlcm` <- function(x, standard.scale = FALSE,
	transpose = FALSE, SD.scale = FALSE, ...){
	m <- if (transpose) t(x$pscale) else x$pscale
   if (standard.scale){
   	dm <- dim(m)
   	mx <- max(m[dm[1], ])
   	matlines(m/mx, ...)
   	} else
   	matlines((if (SD.scale) 2 else 1) * m, ...)		
}

`points.mlcm` <- function(x, standard.scale = FALSE,
	transpose = FALSE, SD.scale = FALSE, ...){
	m <- if (transpose) t(x$pscale) else x$pscale
   if (standard.scale){
   	dm <- dim(m)
   	mx <- max(m[dm[1], ])
   	matpoints(m/mx, ...)
   	} else
   	matpoints((if (SD.scale) 2 else 1) * m, ...)		
}
