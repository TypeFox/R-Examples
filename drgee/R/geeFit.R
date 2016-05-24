geeFit <-
    function(y, x, link = c("identity", "log", "logit")) {

	link <- match.arg(link)

        eq.x <- x

	if (link == "identity") {
            fit <- try( glm.fit(x = x, y = y, family = gaussian()) )
	} else if (link == "log") {
            fit <- try( glm.fit(x = x, y = y, family = quasipoisson()) )
	} else if(link == "logit") {
            fit <- try( glm.fit(x = x, y = y, family = quasibinomial()) )
	}

	if (class(fit) == 'try-error') {

            d.res = matrix(rep(NA,nrow(x) * ncol(x)), nrow = nrow(x))

            return( list(coefficents = rep(NA, ncol(x)),
                         res = rep(NA,nrow(y)),
                         d.res = d.res,
                         eq.x = eq.x,
                         optim.object = NULL,
                         naive.var = NULL) )

        } else {

            if (link == "identity") {
		d.res <- -x
            } else if (link == "log") {
		d.res <- -x * fit$fitted.values
            } else if (link == "logit") {
		d.res <- -x * fit$fitted.values / (1 + exp(fit$linear.predictors))
            }

            colnames(d.res) <- colnames(x)

            return(list(coefficients = fit$coefficients,
                        res = as.vector(y - fit$fitted.values),
                        d.res = d.res,
                        eq.x = x,
                        optim.object = NULL,
                        naive.var = NULL))
        }

    }
