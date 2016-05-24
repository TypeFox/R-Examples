init.plot <-
function (file.name){
	if  (file.name == "r0002")
		file.name = "Vijver"
	else{
		sub = substr(file.name,1,3)	
		if(sub == "gse")
			file.name = toupper(file.name)
	}	
	plot(c(0,10),c(-1,10),type = "n",xlim = c(0,1),ylim = c(0,1), xlab = "1-specificity", ylab = "Sensitivity", main = paste(file.name, "10-fold Cross-validation", sep = "\n"), lty = 3) 
	abline(0,1, col = "red")
}

