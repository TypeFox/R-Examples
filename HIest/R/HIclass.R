HIclass <-
function(G,P,type){
	N <- nrow(G)
	Nindex <- 1:N
	if(type=="codominant") { N <- N/2; Nindex <- rep(1:N,each=2) }
	g.out <- matrix(nrow=N,ncol=6)
		for(i in 1:N){
		if(i/10 == round(i/10)) {cat("\n",i)}else{cat(" *")}
		g <- G[Nindex==i,]
		class100 <- HILL(par=c(0,0),g,P,type)
		class010 <- HILL(par=c(0.5,1),g,P,type)
		class121 <- HILL(par=c(0.5,0.5),g,P,type)
		class110 <- HILL(par=c(0.25,0.5),g,P,type)
		class011 <- HILL(par=c(0.75,0.5),g,P,type)
		class001 <- HILL(par=c(1,0),g,P,type)
		g.out[i,] <- c(class100,class010,class121,class110,class011,class001)
		}
	Class <- apply(g.out,1,which.max)
	fun.fn <- function(x){x[rank(x,ties.method="first")==6]-x[rank(x,ties.method="first")==5]}
	LLD <- apply(g.out,1,fun.fn)
	
	colnames(g.out) <- c("class100","class010","class121","class110","class011","class001")
	data.frame(g.out,Best=c("class100","class010","class121","class110","class011","class001")[Class],LLD=LLD)
	}

