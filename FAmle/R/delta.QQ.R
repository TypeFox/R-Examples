delta.QQ <-
function(model,alpha=.1,ln=FALSE)
	sapply(as.list(model$x.info[,'Emp']),function(h)
		distr(c(alpha/2,1-alpha/2),'norm',as.numeric(delta.Q(h,model,ln)),type='q'))

