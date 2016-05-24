FHtesticp.formula <-
function (formula, data, subset, na.action,...)
{
	call <- match.call()
	m <- match.call(expand.dots = FALSE)
	m[[1]] <- as.name("model.frame")
	m$... <- NULL
	m <- eval(m, parent.frame())
	Y <- model.extract(m, "response")
	nmf <- length(m)
	if (nmf != 2)
      stop("formula should be in the form y~x, where y may be a Surv object")
	if (!is.Surv(Y)){
	  if (is.numeric(Y) & is.vector(Y)){
        Y <- Surv(Y, rep(1, length(Y)))
      } else
          stop("Response must be a survival object or numeric vector")
	}
	LR <- SurvLR(Y)
	if(!is.numeric(m[[2]])){
      group <- as.character(m[[2]])
    } else
        group <- m[[2]]
	out <- do.call("FHtesticp", c(list(L = LR$L, R = LR$R, group = group), list(...)))
	group.name <- as.character(names(m)[2])
	formula.name <- as.character(names(m)[1])
	if((sum(nchar(group.name))+sum(nchar(formula.name))) > 60){
		group.name <- "group"
		formula.name <- "Surv(L, R, type = \"interval2\")"
	}
	out$data.name <- paste("Data:", paste(formula.name, group.name, sep = " by "))
	out$n = table(group)
	names(out$n) <- paste(group.name, "=", names(out$n), sep="")
	out$call <- call
	out
}
