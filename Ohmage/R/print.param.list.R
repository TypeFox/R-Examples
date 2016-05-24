print.param.list <- function(HTTPPARAMS){
	for(i in 1:length(HTTPPARAMS)){
		if(class(HTTPPARAMS[[i]]) == "FileUploadInfo"){
			cat(names(HTTPPARAMS[i]), " (file):\n", HTTPPARAMS[[i]]$filename, "\n\n")
		} else {
			cat(names(HTTPPARAMS[i]), ":\n", HTTPPARAMS[[i]], "\n\n")
		}
	}
}
