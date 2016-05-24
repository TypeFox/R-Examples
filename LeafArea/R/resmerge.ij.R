resmerge.ij <- function(path,prefix="\\.|-"){
	temp0 <- readtext.ij(path)
	temp <- sapply(temp0,sum)
	temp.data <- data.frame(file.name=names(temp),size=temp)
	temp.data$file.name <-factor(sapply(strsplit(as.character(temp.data$file.name), prefix), "[",1))
	res <- aggregate(temp.data["size"],temp.data["file.name"],sum)
	names(res) <- c("sample","total.leaf.area")
	return(res)	
}
