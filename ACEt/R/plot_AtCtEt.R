plot_AtCtEt <- function(AtCtEt, boot=FALSE)
{
	if(class(AtCtEt)!='AtCtEt_model')
	{
		stop('The first parameter must be an object obtained from the AtCtEt function.')
	}

	if(boot==TRUE)
	{
		if(is.null(AtCtEt$boot)==TRUE)
		{
		stop('Please first run the AtCtEt model with bootstrapping.')
		}
	}
	
	model_cur <- AtCtEt
	
	l_a <- model_cur$n_beta_a
	l_c <- model_cur$n_beta_c
	l_e <- model_cur$n_beta_e

	#pheno_m <- c(t(data_m[,1:2]))
	#pheno_d <- c(t(data_d[,1:2]))
	#T_m <- rep(data_m[,3], each=2)
	#T_d <- rep(data_d[,3], each=2)
  n_beta <- 1:(l_a+l_c+l_e)
	order <- 3
	x <- seq(from=model_cur$min_t, to=model_cur$max_t, length.out=500)
	if(model_cur$n_beta_a>1)
	{
		bb_a <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	}else{
		bb_a <- splineDesign(model_cur$knots_a, x = x, ord=1, outer.ok = TRUE)
		if(model_cur$beta_a[1]==-Inf)
		{n_beta[1]=0}
	}
	
	if(model_cur$n_beta_c>1)
	{
		bb_c <- splineDesign(model_cur$knots_c, x = x, ord=order, outer.ok = TRUE)
	}else{
		bb_c <- splineDesign(model_cur$knots_c, x = x, ord=1, outer.ok = TRUE)
		if(model_cur$beta_c[1]==-Inf)
		{n_beta[model_cur$n_beta_a+1]=0}
	}
	
	if(model_cur$n_beta_e>1)
	{
		bb_e <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
	}else{
		bb_e <- splineDesign(model_cur$knots_e, x = x, ord=1, outer.ok = TRUE)
	}

	points_a <- exp(bb_a%*%model_cur$beta_a)
	points_c <- exp(bb_c%*%model_cur$beta_c)
	points_e <- exp(bb_e%*%model_cur$beta_e)

	#fisher <- solve(model_cur$hessian[2:(1+l_a+l_c),2:(1+l_a+l_c)])
	fisher <- solve(model_cur$hessian[n_beta,n_beta])
	
	max_v <- max(points_c, points_a, points_e)*1.2
	plot(range(x), c(0,max_v), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance curves of the A, C, and E components")
	
	index <- 1
	# bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	lines(x, points_a, col = "red", lwd = 2)
	lines(x, points_c, col = "blue", lwd = 2)
	lines(x, points_e, col = "pink", lwd = 2)

	if(boot == TRUE)
	{
		lines(AtCtEt$boot$x, AtCtEt$boot$lower.ci_a, col = "orange" ,lty = 2 , lwd = 0.6)
		lines(AtCtEt$boot$x, AtCtEt$boot$upper.ci_a, col = "orange" ,lty = 2 , lwd = 0.6)
		polygon(c(AtCtEt$boot$x, rev(AtCtEt$boot$x)),c(AtCtEt$boot$upper.ci_a, rev(AtCtEt$boot$lower.ci_a)),col='grey',border = NA, lty=3, density=20)
		lines(AtCtEt$boot$x, AtCtEt$boot$lower.ci_c, col = "green" ,lty = 2 , lwd = 0.6)
		lines(AtCtEt$boot$x, AtCtEt$boot$upper.ci_c, col = "green" ,lty = 2 , lwd = 0.6)
		polygon(c(AtCtEt$boot$x, rev(AtCtEt$boot$x)),c(AtCtEt$boot$upper.ci_c, rev(AtCtEt$boot$lower.ci_c)),col='grey',border = NA, lty=3, density=20)
		lines(AtCtEt$boot$x, AtCtEt$boot$lower.ci_e, col = "yellow" ,lty = 2 , lwd = 0.6)
		lines(AtCtEt$boot$x, AtCtEt$boot$upper.ci_e, col = "yellow" ,lty = 2 , lwd = 0.6)
		polygon(c(AtCtEt$boot$x, rev(AtCtEt$boot$x)),c(AtCtEt$boot$upper.ci_e, rev(AtCtEt$boot$lower.ci_e)),col='grey',border = NA, lty=3, density=20)
	
	}else{
	
	if((l_a>1)|(model_cur$beta_a[1]!=-Inf))
	{
	
	fisher_a <- fisher[index:l_a,index:l_a]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))
	flag <- 0

	for(i in 1:length(x))
	{
		delta <- t(bb_a[i,])%*%fisher_a%*%bb_a[i,]
		if(delta>=0)
		{
		sd[i] <- sqrt(delta)
		lower[i] <- sum(bb_a[i,]*model_cur$beta_a) - 1.96*sd[i]
		upper[i] <- sum(bb_a[i,]*model_cur$beta_a) + 1.96*sd[i]
		}else{flag <- 1}
	}
	#if(l_a>1)
	#{
		lower <- exp(lower)
		upper <- exp(upper)
	#}
	lower <- ifelse(lower<0, 0, lower)
	upper <- ifelse(upper>max_v, max_v, upper)
	index <- index + l_a
	if(flag == 0)
	{
	lines(x, lower, col = "orange" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "orange" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)
	}else{warning('The variance of one of the estimates from the Delta method for the A component is negative. Please try the bootstrap method for the confidence interval or use a different model.')}
	}
	
	
	if((l_c>1)|(model_cur$beta_c[1]!=-Inf))
	{
	fisher_c <- fisher[index:(index+l_c-1),index:(index+l_c-1)]
	
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))
	flag <- 0
	for(i in 1:length(x))
	{
		delta <- t(bb_c[i,])%*%fisher_c%*%bb_c[i,]
		if(delta>=0)
		{
		sd[i] <- sqrt(delta)
		lower[i] <- sum(bb_c[i,]*model_cur$beta_c) - 1.96*sd[i]
		upper[i] <- sum(bb_c[i,]*model_cur$beta_c) + 1.96*sd[i]
		}else{flag <- 1}
	}
	#if(l_c>1)
	#{
		lower <- exp(lower)
		upper <- exp(upper)
	#}
	lower <- ifelse(lower<0, 0, lower)
	upper <- ifelse(upper>max_v, max_v, upper)
	index <- index + l_c
	if(flag == 0)
	{
	lines(x, lower, col = "green" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "green" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)
	}else{warning('The variance of one of the estimates from the Delta method for the C component is negative. Please try the bootstrap method for the confidence interval or use a different model.')}

	}
	
	fisher_e <- fisher[index:(index+l_e-1),index:(index+l_e-1)]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(t(bb_e[i,])%*%fisher_e%*%bb_e[i,])
		lower[i] <- sum(bb_e[i,]*model_cur$beta_e) - 1.96*sd[i]
		upper[i] <- sum(bb_e[i,]*model_cur$beta_e) + 1.96*sd[i]
	}
	#if(l_e>1)
	#{
		lower <- exp(lower)
		upper <- exp(upper)
	#}
	lower <- ifelse(lower<0, 0, lower)
	upper <- ifelse(upper>max_v, max_v, upper)
	lines(x, lower, col = "yellow" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "yellow" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)
	
	} # boot == FALSE
	
	legend(x[1], max_v, c('Additive genetic component','Common environmental component', 'Unique environmental component'), col = c('red','blue','pink'), lty=c(1,1,1), lwd=c(2,2,2))

}