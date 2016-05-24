select.by.ind <-
function(data.obj, geno.or.pheno = "pheno", expr){
	
	mat.to.sub.test <- grep("p", geno.or.pheno)
	if(length(mat.to.sub.test) > 0){
		sub.mat <- data.obj$pheno
		}else{
			sub.mat <- data.obj$geno
			}
	
	expr.pieces <- strsplit(expr, "\\ ")
	if(length(expr.pieces[[1]]) != 3){
		stop("Expression must be in the format 'colname comparison value'")
		}
	
	if(geno.or.pheno == "geno"){
		col.locale <- which(data.obj$marker.names == as.character(expr.pieces[[1]][1]))
		}else{
		col.locale <- which(colnames(sub.mat) == as.character(expr.pieces[[1]][1]))
		}
	if(length(col.locale) == 0){
		stop("I can't find the column name: ", expr.pieces[[1]][1])
		}
		
	if(length(col.locale) > 1){
		stop("There is more than one column that matches the string: ", expr.pieces[[1]][1])
		}

	#if the final element in the expression is a number
	if(!is.na(as.numeric(expr.pieces[[1]][3]))){
		cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,col.locale]), as.numeric(expr.pieces[[1]][3]))
		}else{ #otherwise
			cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,col.locale]), expr.pieces[[1]][3])
			}
			
	vals.locale <- which(eval(cl))	
	
	if(length(vals.locale) == 0){
		stop("There are no individuals that match this expression.")
		}
	
	if(length(vals.locale) < dim(data.obj$pheno)[1]){
		cat(dim(data.obj$pheno)[1] - length(vals.locale), "individuals were removed.\n")
		}
	
	data.obj$pheno <- data.obj$pheno[vals.locale,]
	data.obj$geno <- data.obj$geno[vals.locale,]
	if(!is.null(data.obj$raw.pheno)){
		data.obj$raw.pheno <- data.obj$raw.pheno[vals.locale,]
		}
	if(!is.null(data.obj$ET)){
		data.obj$ET <- data.obj$ET[vals.locale,]
		}	
	
	return(data.obj)
	
}
