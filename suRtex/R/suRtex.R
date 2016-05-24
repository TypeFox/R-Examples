## Title: LaTeX descriptive statistic reporting for survey data
## Author: Dustin Landers
## Version: 1.0
## Date: July 20, 2013

suRtex <- function(data,mean=FALSE,median=FALSE,sd=FALSE,n=TRUE,sub='',digits=2,startdoc=FALSE,enddoc=FALSE) {

## Returns N or the number of non-NA values in the matrix or data frame for the ith column
return.n <- function(data,i) {
	n <- sum(!is.na(data[,i]))
	return(n)
}

## Returns the ith column proportion in the form of a percent based on the data frame or matrix and jth level specified
return.prop <- function(data,i,j) {
	tab <- prop.table(table(data[,i]))*100
	prop <- tab[c(paste(j))]
	if (is.na(prop)) {
		prop <- 0.0
	}
	return(prop)
}

## Tests to see if all levels of column x are equivalent to column y
test.true.ind <- function(x,y) {
	j <- levels(as.factor(x))
	k <- levels(as.factor(y)) 
	d <- j==k
	final <- TRUE
	for (i in 1:length(d)) {
		if (d[i]==FALSE) {
			final <- FALSE
		}
	}
	return(final)
}

## Tests whether all columns in data frame or matrix data have equivalent levels
test.true.whole <- function(data) {
	final <- TRUE
	j <- ncol(data)
	for (i in 1:ncol(data)) {
		if (test.true.ind(data[,i],data[,j])==FALSE) {
			final <- FALSE
		}
	}
	return(final)
}

## Returns column i in data frame or matrix data with the most number of levels
find.most <- function(data) {
	most <- 0
	col <- NA
	for (i in 1:ncol(data)) {
		len <- length(levels(as.factor(data[,i])))
		if (len > most) {
			most <- len
			col <- i
		}
	}
	return(col)
}
	
## Test if data is a matrix or data frame and stop if not
if (class(data)!="data.frame") {
	if (class(data)!="matrix") {
	stop('Data must be matrix or data frame')
	}
}
	
## Test if all columns in data frame have the same number of levels
test <- test.true.whole(data)
if (test==FALSE) {
	warning(paste('Variables in',class(data),'have different level vector lengths'))
}

## Create statlist for warnings
statlist <- c('')
if (mean==TRUE) {
	statlist <- paste(statlist,'Mean ')
}
if (median==TRUE) {
	statlist <- paste(statlist,'Median ')
}
if (sd==TRUE) {
	statlist <- paste(statlist,'Standard deviation ')
}

## If mean, median, or sd is chosen, warn if variables are not numeric
for (i in 1:ncol(data)) {
	if (is.numeric(data[,i])==FALSE) {
		warning(paste('Variable with name:',names(data)[i],'is not numeric and must be coerced to numeric for the following statistics: ',statlist))
	}
}

## Set digit length
options(digits=digits)

## Convert all columns in data frame or matrix to numeric for calculation
for (i in 1:ncol(data)) {
	data[,i] <- as.factor(data[,i])
}

## Find column in data frame or matrix with the most number of levels and record vector of level names from that column as well as a count for the length of that vector
most <- find.most(data)
lev <- levels(as.factor(data[,most]))
level.count <- length(lev)

## Find number of stats indicated as TRUE
stat.count <- 0
if (mean==TRUE) {
	stat.count <- stat.count + 1
}
if (median==TRUE) {
	stat.count <- stat.count + 1
}
if (sd==TRUE) {
	stat.count <- stat.count + 1
}
if (n==TRUE) {
	stat.count <- stat.count + 1
}

## Begin cat functions for LaTeX code
## Cat functions for if startdoc==TRUE
if (startdoc==TRUE) {
	cat('\\documentclass[10pt]{article} \n')
	cat('\\usepackage[english]{babel} \n')
	cat('\\usepackage{amsmath} \n')
	cat('\\usepackage{graphicx} \n')
	cat('\\usepackage{rotating} \n')
	cat('\\usepackage{natbib} \n')
	cat('\\begin{document} \n')
}

cat('\\begin{sidewaystable}\n')
cat('\\centering\n')
cat('\\begin{tabular}{l',rep('c',level.count),'|',rep('c',stat.count),'}\n')

## Cat functions for various statistics on header
cat('Survey Item')
	for (i in 1:level.count) {
		cat(' &',lev[i])
	}
	if (mean==TRUE) {
		cat(' & Mean')
	}
	if (median==TRUE) {
		cat(' & Median')
	}
	if (sd==TRUE) {
		cat(' & StdDev')
	}
	if (n==TRUE) {
		cat(' & N')
	}
cat(' \\\\ \n')
cat('\\hline \n')
cat('\\hline \n')

## Begin loop for cat functions for each column in the data frame or matrix
for (i in 1:ncol(data)) {

## Cat name of column
cat(names(data)[i],' &',sep='')	
	
	## Begin loop for each j level
	for (j in 1:level.count) {
	## Cat proportion
		cat(return.prop(data,i,lev[j]),'\\%',sep='')
		if (j!=level.count) {
			cat(' & ')
		}
	## End j loop	
	}	

## Put extra & if additional statistics are present
if (mean==TRUE | median==TRUE | sd==TRUE | n==TRUE) {
	cat(' & ')
}	
	
## Cat mean
if (mean==TRUE) {	
	cat(mean(as.numeric(data[,i]),na.rm=TRUE),sep='')
	if (median==TRUE | sd==TRUE | n==TRUE) {
		cat(' &',sep='')
	}
}

## Cat median
if (median==TRUE) {
	cat(median(as.numeric(data[,i]),na.rm=TRUE),sep='')
	if (sd==TRUE | n==TRUE) {
		cat(' &',sep='')
	}
}	

## Cat standard deviation
if (sd==TRUE) {
	cat(sd(as.numeric(data[,i]),na.rm=TRUE),sep='')
	if (n==TRUE) {
		cat(' &',sep='')
	}
}

## Cat N
if (n==TRUE) {
	cat(return.n(data,i),sep='')
}

## To finish each ith column, or row of table
cat(' \\\\ \n',sep='')
cat('\\hline \n')
	
## End main i loop	
}

## Finish LaTeX code
cat('\\end{tabular} \n')
cat('\\caption{',sub,'} \n',sep='')
cat('\\end{sidewaystable} \n')

## Cat function if ending document
if (enddoc==TRUE) {
	cat('\\end{document} \n')
}
## End of suRtex function
}