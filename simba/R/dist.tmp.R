"dist.tmp" <- 
function(x, y, method = "euclidean", margin = 1, adjust=TRUE) {	
	if (!is.na(pmatch(method, "euclidean"))) 
        method <- "euclidean"
    METHODS <- c("manhattan", "euclidean", "bray", "canberra", 
        "kulczynski", "gower", "jaccard")
    method <- pmatch(method, METHODS)
    if(margin==2){
        x <- t(x)
        y <- t(y)
    }
    if(is.vector(x)){
        x <- data.frame(x)
    }
    if(is.vector(y)){
        y <- data.frame(y)
    }
    df1 <- x
	df2 <- y
	nms1 <- names(data.frame(df1))
    nms2 <- names(data.frame(df2))
	n.spc1 <- ncol(df1)
	n.spc2 <- ncol(df2)
#	this rownames are valid if adjust=FALSE
    plots <- rownames(data.frame(df1))
    if (adjust) {
        df.ges <- merge(df1, df2, by=0, suffixes=c(".xxx", ".yyy"))
        row.names(df.ges) <- df.ges[,1]
        df.ges <- df.ges[,-1]
        spc.nms <- names(df.ges)
        names(df.ges) <- c(1:ncol(df.ges))
        df1 <- df.ges[,-grep(".yyy", spc.nms)]
        df2 <- df.ges[,-grep(".xxx", spc.nms)]
#   before re-renaming, the variables (species) from the each other matrix have to be filled in an assigned with zeros
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
	    plots <- rownames(df.ges)
	}
    pa1 <- ifelse(df1>0, 1, 0)
    pa2 <- ifelse(df2>0, 1, 0)
    a <- rowSums(pa1*pa2)
	b <- rowSums((pa1==1) & (pa2==0))
	c <- rowSums((pa1==0) & (pa2==1))
	d <- rowSums(!((pa1==0) & (pa2==0)))
	o1 <- (rowSums(pa1)!=0) & (rowSums(pa2)==0)
	o2 <- (rowSums(pa1)==0) & (rowSums(pa2)!=0)
    
    if(method==1){
        dis <- rowSums(abs(df1-df2))
    }
    if(method==2){
        dis <- sqrt(rowSums((df1-df2)^2))
    }
    if(method==3){
        dis <- rowSums(abs(df1-df2))/rowSums(df1+df2)
        dis[d==0] <- 0
        ##dis[o1==TRUE] <- -dis[o1==TRUE]
    }
    if(method==4){
        dis <- (rowSums(abs(df1-df2))/rowSums(df1+df2))/apply((df1-df2), 1, function(x) sum(x!=0))
    }
	if(method==5){
	   dis <- 1 - (0.5*((rowSums(pmin(df1,df2))/rowSums(df1)) + (rowSums(pmin(df1,df2))/rowSums(df2))))
	}
	if(method==6){
	   maxmin <- apply(pmin(df1,df2), 1, max) - apply(pmin(df1,df2), 1, min)
	   maxmin <- matrix(maxmin, nrow(df1), ncol(df2))
	   dis <- rowSums(abs(df1-df2)/maxmin)
	}
	if(method==7){
        dis <- rowSums(abs(df1-df2))/rowSums(df1+df2)
        dis <- (2*dis)/(1+dis)
        dis[d==0] <- 0
        dis[o1==TRUE] <- -dis[o1==TRUE]
	}
	   
    return(dis)
}