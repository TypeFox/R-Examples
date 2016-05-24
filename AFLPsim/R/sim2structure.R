sim2structure <-
  function (x,filename) {
  x$S<-NULL
  x$SelMarkers<-NULL
  x<-x[!is.na(x)]
if(class(x[[1]])=="list"){
  x<-unlist(x,recursive=F)
  x <- lapply(x, function(x)x[,-1 ])
}
N<-length(x)
N2<-N
raw<- do.call("rbind", x)
ind<-rownames(raw)
popnames<-names(x)
n<-numeric(N2)
for (i in 1:N2){ 
  n[i]<-nrow(x[[i]])
}

pop<-rep(popnames[1:N2],n[1:N2])
col<-c("ind","pop",colnames(raw))
mat<-as.data.frame(cbind(ind,pop,raw))
colnames(mat)<-col
  manb <- dim(mat)[2]-2
  innb <- dim(mat)[1]
  mat <- mat[1:(innb), ]
  mat <- mat[order(mat[ ,2]), ]
  popnames <- as.vector(mat[ , 2])
	indnames <- as.vector(mat[, 1])
	matm <- as.matrix(mat[ , (3:(manb+2))])
	poptab <- table(popnames)
	nbpop <- length(poptab)
	numbers <- 0
	for (n in 1:nbpop){
	 numbers <- c(numbers, c(rep(n, poptab[n])))} 
   numbers <- numbers[2:(innb + 1)]
		cat(c(rep(0, manb)), "\n", file= filename,  sep="\t")
 	for (i in 1: innb) {
	cat(indnames[i], numbers[i], matm[i , ], "\n", file= filename, append=TRUE, sep="\t")	
	cat(indnames[i], numbers[i], matm[i , ], "\n", file= filename, append=TRUE, sep="\t")	
	}
}
