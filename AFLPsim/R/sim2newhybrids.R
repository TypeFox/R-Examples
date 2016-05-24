sim2newhybrids <-
function (x,filename) {
  x<-x[!is.na(x)]
  if(class(x[[1]])=="list"){
    x<-unlist(x,recursive=F)
    x <- lapply(x, function(x)x[,-1 ])
  }
  N<-length(x)
  N2<-N-1
  if (names(x)[N]=="S")
    x$S<-NULL
  admix.gen <-do.call("rbind",x[c(3:N2)])
	raw<-rbind(x$PA,x$PB,admix.gen)
  raw[raw==1]<-"+"
  raw[raw==0]<-"-"
  ind<-rownames(raw)
	ninds<-length(ind)
	ind2<-seq(1,ninds,1)
	category<-c(rep("z0",nrow(x$PA)),rep("z1",nrow(x$PB)),rep("",nrow(admix.gen)))
	pop<-c(rep("1",nrow(x$PA)),rep("2",nrow(x$PB)),rep("3",nrow(admix.gen)))
	col<-c("ind","pop",colnames(raw))
	mat<-as.data.frame(cbind(ind2,pop,raw))
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
	
	cat("NumIndivs ", ninds, "\n", file=filename, sep="")
	cat("NumLoci ", manb, "\n", "Digits 1", "\n", "Format Lumped", "\n", file=filename, sep="", append=TRUE)
	cat("LocusNames", col[3:(manb+2)], "\n", "\n", file=filename, sep="\t", append=TRUE)   
	
	for (i in 1:npop) {
	  for (j in 1: popsizes[i]) {
	    cat(inds[n + j], "\t", category[n + j], "\t", matm[(n + j), ], "\n", file=filename, append=TRUE)}
	  popnames[i] <- pops[n+1]
	  n <- n + popsizes [i]
	}
}
