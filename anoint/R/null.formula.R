make.formula <- function(string){
	formula(paste(string[2],string[1],string[3],sep="",collapse=""))
}

make.null.formula <- function(formula){
	string <- sub("\\*","+",formula)
make.formula(string)
}
