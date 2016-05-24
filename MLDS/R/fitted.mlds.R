`fitted.mlds` <-
function(object, ...) {
#object, object of class mlds
	if (object$method == "glm")
		ans <- fitted(object$obj) else
		{fam <- binomial(link = object$link)
		 d <- object$data
		 s <- object$sigma
		 psc <- object$pscale
		 del <-  matrix(psc[unlist(d[, -d$resp])],
#		 matrix(psc[unlist(subset(d, select = -resp))], 
		 	ncol = 4) %*% c(1, -1, -1, 1)
		z <- del/s
		ans <- fam$linkinv(z)
		}
	as.vector(ans)
	}

`fitted.mlbs` <- function(object, ...){
	 ans <- fitted(object$obj)
	 as.vector(ans)
	}