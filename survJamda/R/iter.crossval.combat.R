iter.crossval.combat <-
function (data, surv, censor,batchID,ngroup=10,plot.roc = 0,method = "none",gn.nb = 100){
	if(!exists("batchID"))
		stop("\rSet batchID", call.=FALSE)
	
	niter = ifelse(ngroup == length(surv), 1,10)
	res = NULL

	file.name=deparse(substitute(data)) 
	if (plot.roc)
		init.plot(file.name)

 	data =data[!is.na(surv),]
   	censor= censor[!is.na(surv)]
   	surv= surv[!is.na(surv)]

	cat ("Iteration\tAUC\tHR(CI)\t\tP-val\n")
	for (i in 1:niter){
		new.lst = cross.val.combat(data, surv, censor,batchID,method = "none", gn.nb, plot.roc, ngroup, i)
		res = rbind (res, new.lst)
	}

if(ngroup != length(surv)){
cat ("Avg AUC+/-SD\tHR(CI)\n")
if (plot.roc)
legend (0.55,0.1, legend = paste("AUC+/-SD =", sprintf("%.2f",as.numeric(mean(res[,1], na.rm = TRUE))), "+/-", sprintf("%.2f",sd (res[,1],na.rm = TRUE)), sep = " "), bty = "n")
cat (sprintf("%.2f",as.numeric(mean(res[,1], na.rm = TRUE))),  "+/-", sprintf("%.2f",sd (res[,1],na.rm = TRUE)), "\t", gm(res[,2]), "(", sprintf("%.2f",ci.gm(res[,2])[1]), "-", sprintf("%.2f",ci.gm(res[,2])[2]), ")\n", sep = "")
}
}

