
plot.venn <- function(x, y, ...,
                      small=0.7,
                      showSetLogicLabel=FALSE,
                      simplify=FALSE
                      )
  {
    drawVennDiagram(
                    data=x,
                    small=small,
                    showSetLogicLabel=showSetLogicLabel,
                    simplify=simplify
                    )
  }

## data should be a matrix.
##   - The first column of the matrix is the
##     count of the number of objects with the specified pattern.
##   - The second and subsequent columns contain 0-1 indicators
##     giving the pattern of group membership


drawVennDiagram <-function(data,
                           small=0.7,
                           showSetLogicLabel=FALSE,
                           simplify=FALSE)
  {
	numCircles<-NA
	data.colnames<-NULL
	data.rownames<-NULL

	if(is.matrix(data)) {
		numCircles<-ncol(data)-1
		data.colnames<-colnames(data)[2:(ncol(data))]
		# Order is reverted since later indexing starts with
		# the "lowest bit" and that is expected at the left
		data.rownames<-rownames(data)
	}
	else {
		if (is.list(data)) {
			stop("gplots.drawVennDiagram: This internal function is used wrongly. ",
                             "Please call the function 'venn' with the same arguments, instead.\n")
		}

		warning("drawVennDiagram: Testing only, presuming first argument to specify",
		     "the number of circles to draw.\n")
		numCircles<-data
	}


	m<-(0:(-1+2^numCircles))

	if (! is.matrix(data)) {
		##cat("prepare randomised data\n")
		data<-t(sapply(X=m,FUN=function(v){
			l<-baseOf(v,2,numCircles)
			#print(l)
			return(l)
		}))

		#print(data)

		#data.names<-apply(data,1,function(X){
		#	return(paste(X),collapse="")
		#})
		for(i in m) {
			n<-paste(data[i+1,],collapse="")
			if (is.null(data.rownames)) {
				data.rownames<-n
			}
			else {
				data.rownames<-c(data.rownames,n)
			}
		}
		#print(data.rownames)
		data<-cbind(sample(1:100,size=2^numCircles,replace=TRUE),data)
		#print(data)
		rownames(data)<-data.rownames
		data.colnames<-LETTERS[1:numCircles]
		colnames(data)<-c("num",data.colnames)
	}


	plot.new()
	h<-400
	plot.window(c(0,h), c(0,h), ylab="", xlab="")

	if ((2 <= numCircles && numCircles <= 3) || (4 == numCircles && simplify)) {


		circle <- function(x,y=NULL,r=1) {
			elps=cbind(r*cos(seq(0,2*pi,len=1000)), r*sin(seq(0,2*pi,len=1000)));
			if (!is.null(y)) {
				if (length(x) != length(y))
				  stop("circle: both x and y need to be of same length")
				if (is.matrix(x) && ncol(x)>1)
				  stop("circle: if y is not NULL, then x must not be a matrix")
				x<-cbind(x,y)
			}
			for(i in 1:nrow(x)) {
				ax<-elps[,1]+rep(x[i,1],1000)
				ay<-elps[,2]+rep(x[i,2],1000)
				polygon(ax,ay)
			}
		}

		#nolongerreguired#require(grid)

		##cat("drawing circles\n")
		# draw circles with radius 1.7 equally distributed
		# with centers on a circle of radius 1

		degrees<-2*pi/numCircles*(1:numCircles)

		# scaling factor
		s<-1/8*h
		# radius for circles
		r<-3/12*h

		x<-sapply(degrees,FUN=sin)*s + 0.5*h
		y<-sapply(degrees,FUN=cos)*s + 0.5*h

		##cat("filling data\n")
		circle(x,y,r)

		distFromZero<-rep(NA,2^numCircles)
		degrees<-rep(NA,2^numCircles)

		degrees[(2^numCircles)]<-0
		distFromZero[(2^numCircles)]<-0

		for (i in 0:(numCircles-1)) {
			distFromZero[2^i+1] <- r
			degrees[2^i+1] <- 2*pi/numCircles*i
			d<-degrees[2^i+1]

			#print(data.colnames)

			text(
				# starting from the lowest bit, hence reading
				# lables from the right
				label=data.colnames[numCircles - i],
				x=sin(d)*5/12*h+0.5*h,
				y=cos(d)*5/12*h+0.5*h
			)

		}

		if (4==numCircles) {
			for (i in 0:(numCircles-1)) {
				# Current set bit plus the bit left of it and the bit right of it
				distFromZero[2^i
						+2^((i+numCircles-1)%%numCircles)
						+2^((i+1)%%numCircles)+1] <- 2/12*h
				distFromZero[2^i+1] <- 3.5/12*h
				degrees[2^i
						+2^((i+numCircles-1)%%numCircles)
						+2^((i+1)%%numCircles)+1] <- degrees[2^i+1]
			}
		}

				#degrees[2^i+1] + degrees[2^((i+1)%%numCircles)+1])/2

		if (3 <=numCircles) {
			for (i in 0:(numCircles-1)) {
				distFromZero[(2^i+2^((i+1)%%numCircles))+1]<- 2.2/12*h
				distFromZero[2^i+1] <- 3/12*h
				if (i == (numCircles-1)) {
					degrees[(2^i+2^((i+1)%%numCircles))+1] <- (
						degrees[2^i+1] + 2*pi+ degrees[1+1])/2
				}
				else {
					degrees[(2^i+2^((i+1)%%numCircles))+1] <- (
						degrees[2^i+1] + degrees[2^((i+1)%%numCircles)+1])/2
				}
			}
		}

		for(i in 1:2^numCircles) {
			n<-paste(baseOf((i-1),2,numCircles),collapse="")
			v<-data[n,1]
			d<-degrees[i]
			if (1 == length(d) && is.na(d)) {
				if (v>0) warning("Not shown: ",n," contains ",v,"\n")
			}
			else {
				l<-distFromZero[i]
				x<-sin(d)*l+0.5*h
				y<-cos(d)*l+0.5*h
				#cat("i=",i," x=",x," y=",y," label=",n,"\n")
				l<-v
				if (showSetLogicLabel) l<-paste(n,"\n",v,sep="")
				text(label=l,x=x,y=y)
			}
		}

	}
	else if ( (4 == numCircles && !simplify) || numCircles <= 5 ) {

		# Function to turn and move ellipses/circles
		relocate_elp <- function(e, alpha, x, y){
			phi=(alpha/180)*pi;
			xr=e[,1]*cos(phi)+e[,2]*sin(phi)
			yr=-e[,1]*sin(phi)+e[,2]*cos(phi)
			xr=x+xr;
			yr=y+yr;
			return(cbind(xr, yr))
		}

		lab<-function (identifier, data, showLabel=showSetLogicLabel) {
			r<-data[identifier,1]
			if (showLabel) {
				return(paste(identifier,r,sep="\n"))
			}
			else {
				return(r)
			}
		}

	    if (4 == numCircles) {
	        elps=cbind(162*cos(seq(0,2*pi,len=1000)), 108*sin(seq(0,2*pi,len=1000)));

		#plot(c(0, 400), c(0, 400), type="n", axes=F, ylab="", xlab="");
		polygon(relocate_elp(elps, 45,130,170));
		polygon(relocate_elp(elps, 45,200,200));
		polygon(relocate_elp(elps,135,200,200));
		polygon(relocate_elp(elps,135,270,170));

		text( 35, 315, data.colnames[1],cex=1.5)
		text(138, 347, data.colnames[2],cex=1.5)
		text(262, 347, data.colnames[3],cex=1.5)
		text(365, 315, data.colnames[4],cex=1.5)

	        elps <- cbind(130*cos(seq(0,2*pi,len=1000)),
			80*sin(seq(0,2*pi,len=1000)))

		text( 35, 250, lab("1000",data));
		text(140, 315, lab("0100",data));
		text(260, 315, lab("0010",data));
		text(365, 250, lab("0001",data));

		text( 90, 280, lab("1100",data), cex=small)
		text( 95, 110, lab("1010",data) )
		text(200,  50, lab("1001",data), cex=small)
		text(200, 290, lab("0110",data))
		text(300, 110, lab("0101",data))
		text(310, 280, lab("0011",data), cex=small)

		text(130, 230, lab("1110",data))
		text(245,  75, lab("1101",data),cex=small)
		text(155,  75, lab("1011",data),cex=small)
		text(270, 230, lab("0111",data))

		text(200,150,lab("1111",data))
	    }
	    else if (5 == numCircles) {

	        elps <- cbind(150*cos(seq(0,2*pi,len=1000)),
			60*sin(seq(0,2*pi,len=1000)))

		polygon(relocate_elp(elps, 90,200, 250))
		polygon(relocate_elp(elps, 162,250, 220))
		polygon(relocate_elp(elps, 234,250, 150))
		polygon(relocate_elp(elps, 306,180, 125))
		polygon(relocate_elp(elps, 378,145, 200))

		text( 20, 295, data.colnames[1],cex=1.5)
		text(140, 380, data.colnames[2],cex=1.5)
		text(350, 318, data.colnames[3],cex=1.5)
		text(350,   2, data.colnames[4],cex=1.5)
		text( 50,  10, data.colnames[5],cex=1.5)

		text( 61, 228, lab("10000",data));
		text(194, 329, lab("01000",data));
		text(321, 245, lab("00100",data));
		text(290,  81, lab("00010",data));
		text(132,  69, lab("00001",data));

		text(146, 250, lab("11000",data), cex=small)
		text(123, 188, lab("10100",data), cex=small)
		text(275, 152, lab("10010",data), cex=small)
		text(137, 146, lab("10001",data), cex=small)
		text(243, 268, lab("01100",data), cex=small)
		text(175, 267, lab("01010",data), cex=small)
		text(187, 117, lab("01001",data), cex=small)
		text(286, 188, lab("00110",data), cex=small)
		text(267, 235, lab("00101",data), cex=small)
		text(228, 105, lab("00011",data), cex=small)

		text(148, 210, lab("11100",data),cex=small)
		text(159, 253, lab("11010",data),cex=small)
		text(171, 141, lab("11001",data),cex=small)
		text(281, 175, lab("10110",data),cex=small)
		text(143, 163, lab("10101",data),cex=small)
		text(252, 145, lab("10011",data),cex=small)
		text(205, 255, lab("01110",data),cex=small)
		text(254, 243, lab("01101",data),cex=small)
		text(211, 118, lab("01011",data),cex=small)
		text(267, 211, lab("00111",data),cex=small)

		text(170, 231,lab("11110",data),cex=small)
		text(158, 169,lab("11101",data),cex=small)
		text(212, 139,lab("11011",data),cex=small)
		text(263, 180,lab("10111",data),cex=small)
		text(239, 232,lab("01111",data),cex=small)

		text(204,190,lab("11111",data))
	    }
	}
	else {
		stop(paste("Venn diagrams for ",numCircles," dimensions are not yet supported.\n"))
	}

}
