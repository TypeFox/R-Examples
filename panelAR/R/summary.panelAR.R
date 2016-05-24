summary.panelAR <- function(object,...){
	rdf <- object$df.residual
	rank <- object$rank
	N <- length(object$residuals)
	k <- length(object$aliased)
	df <- c(rank,rdf,k)
    
	# SE
	se <- sqrt(diag(object$vcov))
	
	# test statistics
	coef <- object$coefficients
	t.stat <- (coef)/se
	p.val <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
	tab <- cbind(coef,se,t.stat,p.val)
	dimnames(tab) <- list(names(coef), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    
    # set up hypothesis for wald test
    hyp <- names(coef)
    hyp <- hyp[hyp!="(Intercept)"]
    
    # wald test
    lh <- linearHypothesis(object, hyp, test=c("Chisq", "F"), vcov.=object$vcov, singular.ok=FALSE)
    wald <- c(lh$Chisq[2],lh$Df[2],lh[["Pr(>Chisq)"]][2])
    names(wald) <- c("value","df","Pr(>Chisq)") # wald stat, model dof, p stat

    # check if balanced
    N.panel <- nrow(object$panelStructure$obs.mat)
    N.time <- ncol(object$panelStructure$obs.mat)
    balanced <- ifelse(N.panel*N.time==N,T,F)
    
	# calculate number of obs per panel, number of time observations, balanced vs. unbalanced 
	N.per.panel <- rowSums(object$panelStructure$obs.mat)
    N.min <- min(N.per.panel)
    N.max <- max(N.per.panel)
    N.avg <- N/N.panel # average units per panel # 
    
    # create list with variables that describe panel structure
    panelStruct <- list(N=N,N.panel=N.panel,N.time=N.time,balanced=balanced,N.min=N.min,N.max=N.max,N.avg=N.avg,N.per.panel=N.per.panel)
    
    out <- list(call=object$call,terms=object$terms,coefficients=tab,residuals=object$residuals, aliased=object$aliased, df=df, rho = object$panelStructure$rho, Sigma=object$panelStructure$Sigma, r2=object$r2, wald=wald, vcov=object$vcov, na.action=object$na.action,panelStructure=panelStruct)
    
	class(out) <- "summary.panelAR"
	out
	}