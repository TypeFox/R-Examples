#
# general databel classes util
#

convert_intlogcha_index_to_int <- function(i,object, margin) 
{
	length_1 <- dim(object)[margin]
	if (missing(i)) return(1:length_1)
	#print(c(length_1,i,class(i)))
	if (class(i) == "numeric" || class(i) == "integer") {
		nobs <- as.integer(i);
	} else if (class(i) == "logical") {
		if (length(i) > length_1) stop("length of 'i' (logical) too long")
		new_i <- i
		if (length(i) < length_1) {
			new_i <- c(rep(i,floor(length_1 / length(i))))
			if ((length_1 %% length(i))>0) new_i <- c(new_i,i[1:(length_1 %% length(i))])
		}
		nobs <- as.integer(which(new_i));
		#print(nobs)
		rm(new_i)
	} else if (class(i) == "character") {
		tmp <- i
		new_i <- match(i,dimnames(object)[[margin]])
		if (any(is.na(new_i))) stop(paste("following IDs were not found:",tmp[which(is.na(new_i))],"\n"))
		nobs <- new_i;
		rm(new_i)
	} else {
		stop("class of 'i' must be numeric or logical or character");
	}
	if (any(is.na(nobs))) stop(paste("some names not found:",i[is.na(nobs)],class(i)));
	if (length(nobs) < 1) stop("no matching rows (observations) found");
	if (min(nobs)<=0) stop("all 'i's must be positive integer");
	if (max(nobs)>length_1) stop("'i' out of range")
	return(as.integer(nobs))
}

filevector_type <- function(type)
{
	#internal format data types
	#define UNSIGNED_SHORT_INT 1
	#define SHORT_INT          2
	#define UNSIGNED_INT       3
	#define INT                4
	#define FLOAT              5
	#define DOUBLE             6
	#define CHAR               7
	#define UNSIGNED_CHAR      8
	if (type == "UNSIGNED_SHORT_INT") intype <- 1
	else if (type == "SHORT_INT") intype <- 2
	else if (type == "UNSIGNED_INT") intype <- 3
	else if (type == "INT") intype <- 4
	else if (type == "FLOAT") intype <- 5
	else if (type == "DOUBLE") intype <- 6
	else if (type == "CHAR") intype <- 7
	else if (type == "UNSIGNED_CHAR") intype <- 8
	else stop(paste("Unknown type:",type))
	return(intype)
}