coef.ordinalgmifs <-
function(object, model.select="AIC", ...) {
    if (model.select=="AIC") {
    	model.select = object$model.select
    } else if (model.select=="BIC") {
    	model.select = which.min(object$BIC)
    }
	if (is.null(object$x)) {
		if (dim(object$w)[2]!=0) {
			if (object$probability.model=="Stereotype") {
				coef<-c(object$alpha, object$phi, object$zeta)
			} else {
				coef<-c(object$alpha, object$zeta)
			}
		} else {
			if (object$probability.model=="Stereotype") {
				coef<-c(object$alpha, object$phi)
			} else {
				coef<-c(object$alpha)
			}
		}
	} else if (dim(object$w)[2]!=0) {
			if (object$probability.model=="Stereotype") {
				if (is.null(dim(object$phi))) {
					phi <- object$phi[model.select]
				} else {
					phi <- object$phi[model.select,]
				}
				if (is.null(dim(object$zeta))) {
					zeta <- object$zeta[model.select]
				} else {
					zeta <- object$zeta[model.select,]
				}
				coef<-c(object$alpha[model.select,], phi, zeta, object$beta[model.select,])
			} else {
				if (is.null(dim(object$zeta))) {
					zeta <- object$zeta[model.select]
				} else {
					zeta <- object$zeta[model.select,]
				}
				coef<-c(object$alpha[model.select,], zeta, object$beta[model.select,])
			}
	} else if (object$probability.model=="Stereotype") {
			if (is.null(dim(object$phi))) {
				phi <- object$phi[model.select]
			} else {
				phi <- object$phi[model.select,]
			}
			coef<-c(object$alpha[model.select,], phi, object$beta[model.select,])
			} else {
				coef<-c(object$alpha[model.select,], object$beta[model.select,])
			}
	coef
}
