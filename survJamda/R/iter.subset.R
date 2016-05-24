iter.subset <-
function (data, surv, censor,method = "none", gn.nb = 50, train.nb = 100){
 	data =data[!is.na(surv),]
   	censor= censor[!is.na(surv)]
   	surv= surv[!is.na(surv)]
	
	res = NULL
	iteration.nb = 100
	
	cat ("Iteration\tAUC\tHR(CI)\t\tP-val\n")

	for (i in 1:iteration.nb){
		new.lst = eval.subset(data, surv, censor,i, method, gn.nb, train.nb=train.nb)
		res = rbind (res, new.lst)
	}

	cat ("Avg AUC+/-SD\tHR(CI)\n")

	cat (sprintf("%.2f",mean(res[,1], na.rm = TRUE)),  "+/-", sprintf("%.2f",sd (res[,1],na.rm = TRUE)), "\t", sprintf("%.2f",gm(res[,2])), "(", sprintf("%.2f",ci.gm(res[,2])[1]), "-", sprintf("%.2f",ci.gm(res[,2])[2]), ")\n", sep = "")
        
}

