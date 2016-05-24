continueResampling <- function(cycle.time){
	iter.time <- cycle.time/60 #seconds for iters into minutes.	
	print(paste("It is estimated that your iterations will take",round(iter.time,digits=2),"minutes.",sep=" "))
	if(interactive()){
		valueCaptured <- ""
		if(iter.time > 1){ #greater than 1 minute.
			while(tolower(valueCaptured) != "n" && tolower(valueCaptured) != "y"){
				valueCaptured <- readline("Do you want to proceed: Y/N")	
			}
			# if(tolower(valueCaptured)=="y"){
				# print("Proceeding with resample-based tests. Please take note of the progress bar.")
				# return(TRUE)
			# }
			if(tolower(valueCaptured)=="n"){
				print(paste("Resample-based tests exiting. Skipping tests.",sep=" "))
				return(FALSE)
			}
		}
		return(TRUE)
	}else{
		print("R is not in interactive() mode. Resample-based tests will be conducted. Please take note of the progress bar.")
		return(TRUE)
	}
}