`predict.mlds` <-
function(object, newdata = NULL,
		 type = "link", ...) {
#object, obj of class mlds
	# don't need type "terms"
	miss <- missing(newdata)
	if (object$method == "glm") {
		if (miss) {
			ans <- predict(object$obj, type = type, ...) } else
			{
			if (length(newdata) == 5) newdata <- make.ix.mat(newdata)	
			 ans <- predict(object$obj, newdata = newdata,
			 	 type = type, ...)
			}
	} else
	{
	 fam <- binomial(link = object$link)
	 psc <- object$pscale
	 s <- object$sigma
	 if (miss) d <- object$data else
	 					  d <- newdata
	 del <-  matrix(psc[unlist(d[, -d$resp])],
#	 matrix(psc[unlist(subset(d, select = -resp))], 
		ncol = 4) %*% c(1, -1, -1, 1)
	 z <- del/s
	 if (type == "link") ans <- z else
	 	ans <- fam$linkinv(z)
			}
	 as.vector(ans)
	}

`predict.mlbs`  <- function (object, newdata = NULL, type = "link", ...) {
	miss <- missing(newdata)
	if (miss) {
            ans <- predict(object$obj, type = type, ...)
        }
        else {
        	mx <- max(newdata[, -1])
    d <- within(newdata, {S1 <- factor(S1, levels = seq_len(mx))
    		S2 <- factor(S2, levels = seq_len(mx))
    		S3 <- factor(S3, levels = seq_len(mx))
   		})
   	m.lst <- lapply(names(d[, -1]), function(nm) {
   			f <- as.formula(paste("~", nm))
   			m <- model.matrix(f, d)
   			if (nm == "S2") m <- -2 * m
   			m
    		})
 	m <- Reduce("+", m.lst)
 	dsInc.df <- data.frame(newdata[, 1], m[, -1])
 	names(dsInc.df) <- c("resp", paste("S", 2:mx, sep = ""))
            ans <- predict(object$obj, newdata = dsInc.df, type = type, 
                ...)
        }
	as.vector(ans)
}