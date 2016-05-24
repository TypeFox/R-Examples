fda.matchid <-
function(mat, acov, type, grouplab){
if(missing(mat) || missing(acov) || missing(type)) 
stop("Missing Arguments")
if(class(mat) != "data.frame" && class(mat) != "matrix") 
stop("mat must be a data.frame or a matrix")
if(class(acov) != "data.frame" && class(acov) != "matrix") 
stop("acov must be a data.frame or a matrix")

coln <- colnames(mat)
for(i in 1:length(coln)){
if(substr(coln[i], 1, 1) == "X"){
coln <- sub("X", "", coln)
colnames(mat) <- coln
}
}
numcolnames <- coln

acov <- acov[!is.na(acov[[2]]),]
id <- acov[[1]]
bthnum <- is.numeric(id) + is.numeric(coln)

if(bthnum == 1){
commonid <- intersect(as.numeric(id), as.numeric(numcolnames))
sltmat <- as.matrix(mat[, charmatch(as.numeric(commonid), as.numeric(numcolnames))])
}else{
commonid <- intersect(as.character(id), as.character(numcolnames))
sltmat <- as.matrix(mat[, charmatch(commonid, numcolnames)])
}
acov <- acov[charmatch(commonid, id),]

if(length(commonid) < 1)
stop("IDs do not match")

a <- acov
if(tolower(type) == "factor"){
a[[2]] <- as.factor(a[[2]])
m <- data.frame(id=a[[1]], model.matrix(~a[[2]]))

if(missing(grouplab))
stop("grouplab is required for 'factor' type.")
colnames(m) <- c("id", grouplab)
}else if(tolower(type) == "contin"){
a[[2]] <- as.numeric(a[[2]])
m <- data.frame(id=a[[1]], model.matrix(~a[[2]]))
contcovname <- names(a[2])
colnames(m) <- c("id", "intercept", paste(contcovname))
}else{
stop("type must be 'factor' or 'contin'.")
}

return(list(mat=sltmat, cov=m))
}
