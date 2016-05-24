plot.roc.curves <-
function(x,y, lp, test, file.name, col,normalization, ...)
{

	tmp = survivalROC (Stime = x, status =y, marker=lp, predict.time = 120, span = 0.25*length(test)^(-0.20))

	if(col == "black"){
		plot(tmp$FP,tmp$TP,xlim = c(0, 1), ylim = c(0,1), xlab= "1-specificity", ylab="Sensitivity", main = file.name, col = col,type = "l") 
		abline(0,1, col = "red")
	}
	else			    
		lines(tmp$FP,tmp$TP,col = col)  
	
	switch (normalization,
		"ComBat"=(adjust=0),
		"Zscore1"=(adjust = .05),
		"Zscore2"=(adjust = .1)
	)

		
	legend(0.3,0.2-adjust,legend = paste(normalization, ",AUC=", round(tmp$AUC,2), sep = ""), text.col = col, bty = "n", cex = 1)	
}
