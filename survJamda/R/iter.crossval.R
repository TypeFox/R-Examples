iter.crossval<-
function (data, surv, censor,ngroup=10, plot.roc = 0, method = "none", zscore=0,gn.nb = 100,gn.nb.display = 0){
	res = NULL

	file.name=deparse(substitute(data)) 
	if (plot.roc)
		init.plot(file.name)

   	data =data[!is.na(surv),]
   	censor= censor[!is.na(surv)]
   	surv= surv[!is.na(surv)]

	cat ("Iteration\tAUC\tHR(CI)\t\tP-val\n")
	for (i in 1:ngroup){
		new.lst = cross.val.surv(data, surv, censor,ngroup, i, method, zscore,gn.nb,gn.nb.display,plot.roc)
		res = rbind (res, new.lst)
	}

	if(ngroup != length(surv)){
		cat ("Avg AUC+/-SD\tHR(CI)\n")
		if (plot.roc)
			legend (0.55,0.1, legend = paste("AUC+/-SD =", sprintf("%.2f",as.numeric(mean(res[,1], na.rm = TRUE))), "+/-", sprintf("%.2f",sd (res[,1],na.rm = TRUE)), sep = " "), bty = "n")
		cat (sprintf("%.2f",as.numeric(mean(res[,1], na.rm = TRUE))),  "+/-", sprintf("%.2f",sd (res[,1],na.rm = TRUE)), "\t", gm(res[,2]), "(", sprintf("%.2f",ci.gm(res[,2])[1]), "-", sprintf("%.2f",ci.gm(res[,2])[2]), ")\n", sep = "")
	}
}

