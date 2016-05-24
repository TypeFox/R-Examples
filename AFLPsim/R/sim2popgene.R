sim2popgene <-
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
	pops <- as.vector(mat[ , 2])
	inds <- as.vector(mat[ , 1])
	matm <- as.matrix(mat[ , (3:(manb+2))])
	popsizes <- table (pops)
	npop <- length(popsizes)
	popnames <- vector(mode="character", npop)
		n <- 0
		
	cat("/* AFLP Data Set */", "\n", "Number of populations = ", npop, "\n", file=filename, sep="")
	cat("Number of loci = ", manb, "\n", "Locus name :", "\n", file=filename, sep="", append=TRUE) 
  cat(col[3:(manb+2)], "\n", "\n", file=filename, sep=" ", append=TRUE)   
	
	for (i in 1:npop) {
		cat("name = ", pops[n + 1], "\n", "fis = 0", "\n", file=filename, sep="", append=TRUE) 	
		for (j in 1: popsizes[i]) {
			cat(matm[(n + j), ], "\n", file=filename, sep="", append=TRUE)}
		cat("\n", "\n",  file=filename, append=TRUE)
		popnames[i] <- pops[n+1]
		n <- n + popsizes [i]
		}
}
