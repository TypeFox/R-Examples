"sim.tmp" <- 
function(x, y, method="soer", normalize = FALSE, adjust=TRUE, ...) {
	if (!is.na(pmatch(method, "jaccard"))) 
        method <- "jaccard"
	METHODS <- c("soerensen", "jaccard", "ochiai", "mountford", "whittaker", "lande", "wilsonshmida", "cocogaston", "magurran", "harrison", "cody", "williams", "williams2", "harte", "simpson", "lennon", "weiher", "ruggiero", "lennon2", "rout1ledge", "rout2ledge", "rout3ledge", "sokal1", "dice", "kulcz1insky", "kulcz2insky", "mcconnagh", "manhattan", "simplematching", "margaleff", "pearson", "roger", "baroni", "dennis", "fossum", "gower", "legendre", "sokal2", "sokal3", "sokal4", "stiles", "yule")
	method <- pmatch(method, METHODS)
	if (is.na(method)){
		stop("invalid similarity method")
		}
	if (method == -1){
		stop("ambiguous similarity method")
		}
	if (ncol(x)==3) {
		x <- mama(x)
        x <- x[order(names(x))]
	}
    if (ncol(y)==3) {
		y <- mama(y)
		y <- y[order(names(y))]
    }
	df1 <- ifelse(x>0, 1, 0)
	df2 <- ifelse(y>0, 1, 0)
	nms1 <- names(data.frame(df1))
    nms2 <- names(data.frame(df2))
	n.spc1 <- ncol(df1)
	n.spc1 <- ncol(df1)
#	this rownames are valid if adjust=FALSE
	plots <- row.names(df1)
	if (adjust) {
        df.ges <- merge(df1, df2, by=0, suffixes=c(".xxx", ".yyy"))
        row.names(df.ges) <- df.ges[,1]
        df.ges <- df.ges[,-1]
        spc.nms <- names(df.ges)
        names(df.ges) <- c(1:ncol(df.ges))
        df1 <- df.ges[,-grep(".yyy", spc.nms)]
        df2 <- df.ges[,-grep(".xxx", spc.nms)]
#   before re-renaming, the species from the each other matrix have to be filled in and assigned with zeros
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
	a <- rowSums(df1*df2)
	b <- rowSums((df1==1) & (df2==0))
	c <- rowSums((df1==0) & (df2==1))
#	still need to think about d
	d <- n.spc - a - b - c
	if (normalize) {
		an <- a/(a+b+c)
		bn <- b/(a+b+c)
		cn <- c/(a+b+c)
		a <- an
		b <- bn
		c <- cn
		}
	if (method == 1) {
        dis <- (2*a)/((2*a) + b + c)
    }
    else if (method == 2) {
        dis <- a / (a + b + c)
    }
    else if (method == 3) {
        dis <- a / sqrt((a + b) * (a + c))
    }
    else if (method == 4) {
        dis <- (2 * a) / (a * (b + c) + (2 * b * c))
    }
    else if (method == 5) {
        dis <- ((a + b + c) / ((2 * a + b + c)/2))-1
    }
    else if (method == 6) {
        dis <- (b + c)/2
    }
    else if (method == 7) {
        dis <- (b+c)/((2*a)+b+c) 
    }
    else if (method == 8) {
        dis <- (b+c)/(a+b+c)
    }
    else if (method == 9) {
        dis <- ((2*a)+b+c)*(1-(a/(a+b+c)))
    }
    else if (method == 10) {
        dis <- pmin(b,c) / (pmax(b,c) + a)
    }
    else if (method == 11) {
        dis <- 1-((a*((2*a)+b+c))/(2*(a+b)*(a+c)))
    }
    else if (method == 12) {
        dis <- pmin(b,c) / (a+b+c) 
    }
    else if (method == 13) {
        dis <- ((b*c)+1) / ((((a+b+c) * exp(2)) - (a+b+c)) / 2)
    }
    else if (method == 14) {
        dis <- 1-((2*a) / ((2*a) + b + c)) 
    }
    else if (method == 15) {
        dis <- pmin(b,c) / (pmin(b,c) + a) 
    }
    else if (method == 16) {
        dis <- (2 * abs(b-c)) / ((2*a) + b +c)
    }
    else if (method == 17) {
        dis <- b + c
    }
    else if (method == 18) {
        dis <- a / (a+c)
    }
    else if (method == 19) {
        dis <- 1 - (log((2*a + b + c)/(a + b + c)) / log(2))
    }
    else if (method == 20) {
        dis <- (((a+b+c) * exp(2))/(((a+b+c)*exp(2))-(2*b*c)))-1
    }
    else if (method == 21) {
        dis <- log(2*a+b+c)-((1/(2*a+b+c))*2*a*log(2))-((1/(2*a+b+c))*((a+b)*log(a+b)+(a+c)*log(a+c))) 
    }
    else if (method == 22) {
        dis <- log(2*a+b+c)-((1/(2*a+b+c))*2*a*log(2))-((1/(2*a+b+c))*((a+b)*log(a+b)+(a+c)*log(a+c)))
    	dis <- exp(dis)-1
    }
    else if (method == 23) {
        dis <- a / (a + 2*(b + c))
    }
    else if (method == 24) {
        dis <- a / (pmin((b+a),(c+a))) 
    }
    else if (method == 25) {
        dis <- a / (b+c)
    }
    else if (method == 26) {
        dis <- ((a/2) * ((2*a) + b +c)) / ((a+b)*(a+c)) 
    }
    else if (method == 27) {
        dis <- ((a^2)-(b*c)) / ((a+b)*(a+c))
    }
    else if (method == 28) {
        dis <- (b+c) / (a+b+c+d) 
    }
    else if (method == 29) {
        dis <- (a+d) / (a+b+c+d) 
    }
    else if (method == 30) {
        dis <- (a * (a+b+c+d)) / ((a+b) * (a+c))
    }
    else if (method == 31) {
        dis <- ((a*d) - (b*c)) / sqrt((a + b)*(a + c)*(d + b)*(d + c))
    }
    else if (method == 32) {
        dis <- (a + d) / (a + 2*(b + c) +d)
    }
    else if (method == 33) {
        dis <- ((sqrt(a*d))+c) / ((sqrt(a*d))+b+c+a)
    }
    else if (method == 34) {
        dis <- ((a+d) - (b*c)) / (sqrt((a+b+c+d)*(a+b)*(a+c)))
    }
    else if (method == 35) {
        dis <- ((a+b+c+d) * (-1 * ((a/2)^2))) / ((a+b)*(a+c))
    }
    else if (method == 36) {
        dis <- (a - (b+c)+d) / (a+b+c+d)
    }
    else if (method == 37) {
        dis <- a / (a+b+c+d)
    }
    else if (method == 38) {
        dis <- a / (a*d) / sqrt((a+b)*(a + c)*(d + b)*(d + c))
    }
    else if (method == 39) {
        dis <- ((2*a)+(2*d)) / (a+d+(a+b+c+d))
    }
    else if (method == 40) {
        dis <- (a+d) / (b+c)
    }
    else if (method == 41) {
        dis <- log(((a+b+c+d) * (( abs((a*d)-(b*c)) - ( (a+b+c+d) / 2))^2) / ((a+b)*(a+c) *(b+d)*(c+d))))
    }
    else if (method == 42) {
        dis <- ((a*d)-(b*c)) / ((a+d)+(b*c))
    }
#	d <- as.dist(dis)
#	attr(d, "Size") <- anz
#    attr(d, "Labels") <- zeina
#    attr(d, "method") <- METHODS[method]
#    attr(d, "call") <- match.call()
#    class(d) <- "dist"
#    if (listout) {
#        d <- liste(d, entry=METHODS[method])
#    	}
    return(dis)
}