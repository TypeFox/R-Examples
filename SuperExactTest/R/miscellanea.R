#Compute overlap test and visualize intersections between multiple sets
#Author: Minghui Wang
#minghui.wang@mssm.edu
#
#use barcode (character strings of 0/1) to denote the overlap/non-overlap sets
mkBarcode=function(n){
#n, number of sets
	if(n<1) stop('n must be positive integer')
	barcode=paste(rep('0',n),collapse='')
	for(i in 1:n){
		barcode0=barcode
		substr(barcode0,i,i)='1'
		barcode=c(barcode,barcode0)
	}
	sort(barcode[-1])
}
#
mkBarcode.degree=function(n,degree=NULL){
#n, number of sets
	if(n<1) stop('n must be positive integer\n')
	if(is.null(degree)) degree=1:n
	if(any(degree < 1) || any(degree > n)) stop('Invalid degree value\n')
	maxdeg=max(degree)
	barcode1=paste(rep('0',n),collapse='')
	barcode=c()
	for(i in 1:n){
		barcode0=barcode1
		substr(barcode0,i,i)='1'
		n1=sapply(barcode0,function(a) countCharOccurrences('1',a))
		barcode=c(barcode,barcode0[n1 %in% degree])
		barcode1=c(barcode1,barcode0[n1 < maxdeg])
	}
	sort(barcode)
}
barcodeRelation=function(barcode1,barcode2){
#barcode, character strings of 0/1
#return a matrix, row specifies whether the row entry is a subset of the column entries.
	Mat=matrix(FALSE,length(barcode1),length(barcode2))
	Mat=sapply(barcode1,function(a)
		sapply(barcode2,function(b) is.subSet(a,b))
	)
	t(Mat)
}
is.subSet=function(a,b){
#a,b barcode
#return TRUE if a is b's subset; FALSE otherwise. Eg, a='00011' is a subset of b='00001', while a='00111' is not a subset of b='10001'.
	a1=strsplit(a,'')[[1]] == '1'
	b1=strsplit(b,'')[[1]] == '1'
	if(length(a1) != length(b1)) stop('Input a and b have different numbers of characters\n')
	all(a1[which(b1)])
}
#count number of occurrences a char appears in a string
countCharOccurrences=function(char, s) { #by Uli Kohler
	s2 <- gsub(char,"",s)
	nchar(s) - nchar(s2)
}
