Diff.3 <-
function(i,model,p,ln=FALSE)
	Diff.1(model$par.hat[i],function(g) Diff.2(g,i,model,p,ln))

