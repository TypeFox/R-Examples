plot.time.dep <-
function(x,y, lp, test, file.name, col,...)
{
	roc.fit = NULL
	surv.sorted = sort(x)

	for (t in surv.sorted){
		tmp = survivalROC (Stime = surv.sorted, status =y[order(x)], marker=lp[order(x)], predict.time = t, span = 0.25*length(test)^(-0.20))$AUC
		roc.fit = c(roc.fit, round(tmp,2))

	}

	if(col == "black")
		plot(surv.sorted,roc.fit,xlim = c(0, 120), ylim = c(0,1), xlab= "Time (months)", ylab="AUC", main = file.name, col = col,type = "l") 
	else			
		lines(surv.sorted,roc.fit,col = col)     
}
