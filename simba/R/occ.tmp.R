"occ.tmp" <- 
function(x, y, adjust=TRUE, gen.occ=FALSE, perc=TRUE, nc.acc=FALSE, ...) {	
	if (ncol(x)==3) {
		x <- mama(x)
        x <- x[order(names(x))]
	}
	if (is.vector(y)){
	   stop("y has to be a table with three columns or a species matrix")
	}
    if (ncol(y)==3) {
		y <- mama(y)
		y <- y[order(names(y))]
    }
	df1 <- ifelse(x>0, 1, 0)
	df2 <- ifelse(y>0, 1, 0)
#	df1 <- as.matrix(x)
#	df2 <- as.matrix(y)
    nms1 <- names(data.frame(df1))
    nms2 <- names(data.frame(df2))
	n.spc1 <- ncol(df1)
	n.spc2 <- ncol(df2)
	plots <- row.names(df1)
	if (adjust) {
    df.ges <- merge(df1, df2, by=0, suffixes=c(".xxx", ".yyy"))
    row.names(df.ges) <- df.ges[,1]
    df.ges <- df.ges[,-1]
    spc.nms <- names(df.ges)
    names(df.ges) <- c(1:ncol(df.ges))
    df1 <- df.ges[,-grep(".yyy", spc.nms)]
    df2 <- df.ges[,-grep(".xxx", spc.nms)]
#    vor der rückbenennung müssen noch die arten die aus der jeweils anderen matrix übernommen wurden mit nullen versehen werden:
    only1 <- data.frame((sapply(nms1, grep, spc.nms)))
    only1.slct <- only1[1, apply(only1, 2, diff) == 0]
    only2 <- data.frame((sapply(nms2, grep, spc.nms)))
    only2.slct <- only2[1, apply(only2, 2, diff) == 0]
    df1[,as.character(only2.slct)] <- 0
    df2[,as.character(only1.slct)] <- 0
    names(df1) <- gsub(".xxx", "", spc.nms[as.numeric(names(df1))])
    names(df2) <- gsub(".yyy", "", spc.nms[as.numeric(names(df2))])
    df1 <- df1[,sort(names(df1))]
    df2 <- df2[,sort(names(df2))]
	plots <- row.names(df.ges)
	}
    n.spc <- ncol(df1)
    n.plots <- nrow(df1)
    spec.occ <- ifelse(colSums(df1) > 0, 1, 0) - ifelse(colSums(df2) > 0, 1, 0)
    specvek0 <- names(df1)[spec.occ==0]
    specvek1 <- names(df1)[spec.occ==1]
    specvek2 <- names(df1)[spec.occ==-1]
    if (gen.occ) {
        bac <- summary(as.factor(colSums(df2)-colSums(df1)))
    }
    else {
	a <- colSums(df1*df2)
	b <- colSums((df1==1) & (df2==0))
	c <- colSums((df1==0) & (df2==1))
	if (nc.acc){
	   a <- sum(a)
	}
	else {
       a <- sum(ifelse(a>0, 1, 0))
	}
    attr(a, "names") <- 0
    b <- sapply(split(b, b), length)[-1]
    attr(b, "names") <- as.numeric(attr(b, "names"))*(-1)
    c <- sapply(split(c, c), length)[-1]
	bac <- c(b, a, c)
	bac <- bac[order(as.numeric(names(bac)))]
    }
    if(perc){
        if(nc.acc){
            bac <- bac/n.spc/n.plots
        }
        else{
            bac <- bac/n.spc
        }
    }
    res <- list(bac=bac, n.plots=n.plots, n.spc=n.spc, n.spc1=n.spc1, n.spc2=n.spc2, n.spc1o=length(specvek1), n.spc2o=length(specvek2), spec.names1o=specvek1, spec.names2o= specvek2)
#    noch die artenanzahlen in den ursprungsmatrizen einbeziehen in das ergebnis und auch den call mit zurückgeben?
    return(res)
}