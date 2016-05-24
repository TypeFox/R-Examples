
# Function: add.cutline
# Reads in output from W-NOMINATE and adds a cutting line to existing plot
#   INPUTS: a numeric vector of length 4, cutData
#       midpoint1d<-cutData[1]
#       spread1d<-cutData[2]
#       midpoint2d<-cutData[3]
#       spread2d<-cutData[4]

add.cutline <- function(cutData,weight,lwd=2) {

    slope <- -cutData[2]/(cutData[4]*weight)
    if (is.na(slope)) {
        x <- c(cutData[1],cutData[1])
            y <- c(sqrt(1-cutData[1]^2),-sqrt(1-cutData[1]^2))
                slope <- NA
                intercept <- NA
        }
        else {
                intercept <- -slope*cutData[1]+cutData[3]
                x <- c( (-slope*intercept + sqrt( (slope*intercept)^2 -
            (1+slope*slope)*(intercept*intercept-1)))/(1+slope*slope),
                    (-slope*intercept - sqrt( (slope*intercept)^2 - 
            (1+slope*slope)*(intercept*intercept-1)))/(1+slope*slope) )
            if (is.na(x[1])) {
                warning("Couldn't solve for points on the unit circle!\n")
                x<-NA
                y<-NA
                slope<-NA
                intercept<-NA  
            }             
            else {
                y <- intercept + slope*x
                y[y < -1] <- -sqrt(1-x[y<1]^2)
                y[y >  1] <-  sqrt(1-x[y>1]^2)
            }
        }
    lines(x,y,lwd=lwd)
}

plot.angles <- function(x, main.title="Cutting Line Angles",
        x.title="Angle in Degrees", y.title="Count",dims=c(1,2),...) {

    if(!class(x)=="nomObject")
        stop("Input is not of class 'nomObject'.")
    if(x$dimensions==1)
        stop("All angles in 1D NOMINATE are 90 degrees.")
    if(length(dims)!=2)
        stop("'dims' must be an integer vector of length 2.")

    weight<-x$weight[dims[2]]/x$weight[dims[1]]

    contrained <- ((abs(x$rollcalls[,paste("spread",dims[1],"D",sep="")]) > 0.0 |
                 abs(x$rollcalls[,paste("spread",dims[2],"D",sep="")]) > 0.0)
                 & (x$rollcalls[,paste("midpoint",dims[1],"D",sep="")]**2 +
                 x$rollcalls[,paste("midpoint",dims[2],"D",sep="")]**2) < .95)

    cutvector1 <- na.omit(x$rollcalls[contrained,paste("spread",dims[2],"D",sep="")]*weight/
                    sqrt(x$rollcalls[contrained,paste("spread",dims[1],"D",sep="")]^2
                    + weight^2*x$rollcalls[contrained,paste("spread",dims[2],"D",sep="")]^2))
    cutvector2 <- -1*na.omit(x$rollcalls[contrained,paste("spread",dims[1],"D",sep="")]/
                    sqrt(x$rollcalls[contrained,paste("spread",dims[1],"D",sep="")]^2
                    + weight^2*x$rollcalls[contrained,paste("spread",dims[2],"D",sep="")]^2))
    cutvector1[cutvector2<0] <- -cutvector1[cutvector2<0]
    cutvector2[cutvector2<0] <- -cutvector2[cutvector2<0]
    angles <- atan2(cutvector2,cutvector1)*180/pi
    
    suppressWarnings(hist(angles, breaks=seq(0,180,10),
        main=main.title,
        xlab=x.title,
        ylab=y.title,
        cex.main=1.2,
        cex.lab=1.2,
        font.main=2,
        at=seq(0,180,10)
        ,...))
}

plot.cutlines <- function(x, main.title="Cutting Lines",
        d1.title="First Dimension", d2.title="Second Dimension",
        lines=50,dims=c(1,2),lwd=2,...) {

    if(!class(x)=="nomObject")
        stop("Input is not of class 'nomObject'.")
    if(x$dimensions==1)
        stop("All angles in 1D NOMINATE are 90 degrees.")
    if(length(dims)!=2)
        stop("'dims' must be an integer vector of length 2.")
    if(lines<1)  stop("'Lines' must be less than 1.")

    constrained <- ((abs(x$rollcalls[,"spread1D"]) > 0.0 | abs(x$rollcalls[,"spread2D"]) > 0.0)
        & (x$rollcalls[,"midpoint1D"]**2 + x$rollcalls[,"midpoint2D"]**2) < .95)
    
    cutlineData <- cbind(x$rollcalls[constrained,paste("midpoint",dims[1],"D",sep="")],
                     x$rollcalls[constrained,paste("spread",dims[1],"D",sep="")],
                     x$rollcalls[constrained,paste("midpoint",dims[2],"D",sep="")],
                     x$rollcalls[constrained,paste("spread",dims[2],"D",sep="")])
    cutlineData <- na.omit(cutlineData)

    suppressWarnings(symbols(x=0, y=0, circles=1, inches=FALSE, asp=1,
        main=main.title,
        xlab=d1.title,
        ylab=d2.title,
        xlim=c(-1.0,1.0),
        ylim=c(-1.0,1.0),
        cex.main=1.2,
        cex.lab=1.2,
        font.main=2,
        lwd=2,
        fg="grey",
        frame.plot=FALSE,...))

    if(lines<dim(cutlineData)[1])
	cutlineData <- cutlineData[sample(1:dim(cutlineData)[1],lines),]
    
    suppressWarnings(apply(cutlineData, 1, add.cutline,
        weight=x$weights[dims[2]]/x$weights[dims[1]],lwd=lwd))

}

plot.coords <- function (x, main.title="W-NOMINATE Coordinates",
    d1.title="First Dimension", d2.title="Second Dimension", dims=c(1,2),
    plotBy="party", color=TRUE, shape=TRUE, cutline=NULL, Legend=TRUE,
    legend.x=0.8,legend.y=1,...) {
   
    if(!class(x)=="nomObject")
        stop("Input is not of class 'nomObject'.")
    if(!any(colnames(x$legislators)==plotBy)){
        warning("Variable '", plotBy ,"' does not exist in your W-NOMINATE object.")
	types <- rep("Leg",dim(x$legislators)[1])
    } else {
        types <- x$legislators[,plotBy]
    }
    if(length(dims)!=2 & x$dimensions!=1)
        stop("'dims' must be an integer vector of length 2.")
   
    # determine number of parties
    nparties <- length(unique(types))
    
    # set default colors and shapes
    colorlist <- c("darkblue", "firebrick", "darkcyan", "darkgreen", "darkmagenta", "darkolivegreen", 
    "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen", "darkslateblue", 
    "darkslategray", "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dodgerblue")
    shapes <- rep(c(16,15,17,18,19,3,4,8),3)
    
    # color and shape options
    if (color==FALSE) colorlist <- sample(colors()[160:220],50)
    if (shape==FALSE) shapes <- rep(16,50)

    if(x$dimensions==1){   
        coord1D <- x$legislators[,"coord1D"]
        ranking <- rank(x$legislators[,"coord1D"])
        plot(seq(-1,1,length=length(coord1D)),
                1:length(coord1D),
                type="n",
                cex.main=1.2,
                cex.lab=1.2,
                font.main=2,
                xlab="First Dimension Nominate",
                ylab="Rank",
                main="1D W-NOMINATE Plot")

        if(Legend)  legend(0.67, 0.7*length(coord1D), unique(types), pch=shapes[1:nparties],
                            col=colorlist[1:nparties], cex=0.7)
        for(i in 1:nparties) suppressWarnings(points(coord1D[types==unique(types)[i]],
            ranking[types==unique(types)[i]],pch=shapes[i],col=colorlist[i],cex=1.1,lwd=2))  
    } else {

    #2 Dimensional Case begins here
    coord1D <-  x$legislators[,paste("coord",dims[1],"D",sep="")]
    coord2D <-  x$legislators[,paste("coord",dims[2],"D",sep="")]
    
    # Plotting
    suppressWarnings(symbols(x = 0, y = 0, circles = 1, inches = FALSE,
            asp = 1,
            main=main.title,
            xlab=d1.title,
            ylab=d2.title,
            xlim=c(-1.0,1.0),
            ylim=c(-1.0,1.0),
            cex.main=1.2,
            cex.lab=1.2,
            font.main=2,
            lwd=2,
            fg="grey",
            frame.plot=FALSE,...))

    if(!is.null(cutline)) {
        for(i in 1:length(cutline)){
        if(all(is.na(x$rollcalls[cutline[i],])))
            stop("Roll call for cutline did not meet minimum lopsidedness requirements.")
        add.cutline(c(x$rollcalls[cutline[i],paste("midpoint",dims[1],"D",sep="")],
                    x$rollcalls[cutline[i],paste("spread",dims[1],"D",sep="")],
                    x$rollcalls[cutline[i],paste("midpoint",dims[2],"D",sep="")],
                    x$rollcalls[cutline[i],paste("spread",dims[2],"D",sep="")]),
            weight=x$weights[dims[2]]/x$weights[dims[1]],
            lwd=2)
        }
    }
   
    if(Legend)
        legend(legend.x, legend.y, unique(types), pch=shapes[1:nparties],
	bty="n",col=colorlist[1:nparties], cex=0.7)

    for(i in 1:nparties) suppressWarnings(points(coord1D[types==unique(types)[i]],
        coord2D[types==unique(types)[i]],pch=shapes[i],col=colorlist[i],cex=1.1,lwd=2))
    }
}

plot.scree <- function(x, main.title="Scree Plot", x.title="Dimension",
                        y.title="Eigenvalue",...) {

    if(!class(x)=="nomObject")
        stop("Input is not of class 'nomObject'.")
    if(is.null(x$eigenvalues))
    stop("No eigenvalues exist in this W-NOMINATE object.")
    suppressWarnings(plot(1:20,
        x$eigenvalues[1:20],
    type='o',
        main=main.title,
        xlab=x.title,
        ylab=y.title,
        cex.main=1.2,
        cex.lab=1.2,
        font.main=2,
        lwd=1,
        pch=16,
        at=1:20,...))

}               

plot.nomObject <- function(x,dims=c(1,2),...) {
    if(!class(x)=="nomObject")
        stop("Input is not of class 'nomObject'.")
    if(length(dims)!=2 & x$dimensions!=1)
        stop("'dims' must be an integer vector of length 2.")

    if(x$dimensions==1) {
        par(mfrow=c(1,2))       
        suppressWarnings(plot.coords(x,dims=dims))
        suppressWarnings(plot.scree(x,dims=dims))
    } else {
        par(mfrow=c(2,2))
        suppressWarnings(plot.coords(x,dims=dims))
        suppressWarnings(plot.angles(x,dims=dims))
        suppressWarnings(plot.scree(x,dims=dims))
        suppressWarnings(plot.cutlines(x,dims=dims,lwd=1))
    }
}               

summary.nomObject<-function(object,verbose=FALSE,...){

    if(!class(object)=="nomObject")
        stop("Input is not of class 'nomObject'.")

    cat("\n\nSUMMARY OF W-NOMINATE OBJECT")
    cat("\n----------------------------\n")
    cat("\nNumber of Legislators:\t  ", dim(na.omit(object$legislators))[1],
    " (", dim(object$legislators)[1]-dim(na.omit(object$legislators))[1],
    " legislators deleted)", sep="")
    cat("\nNumber of Votes:\t  ", dim(na.omit(object$rollcalls))[1],
    " (", dim(object$rollcalls)[1]-dim(na.omit(object$rollcalls))[1],
    " votes deleted)", sep="")
    cat("\nNumber of Dimensions:\t ", object$dimensions)

    correctYea<-sum(as.numeric(object$legislators[,"correctYea"]),na.rm=TRUE)
    allYea<-correctYea+sum(as.numeric(object$legislators[,"wrongNay"]),na.rm=TRUE)
    correctNay<-sum(as.numeric(object$legislators[,"correctNay"]),na.rm=TRUE)
    allNay<-correctNay+sum(as.numeric(object$legislators[,"wrongYea"]),na.rm=TRUE)
    cat("\nPredicted Yeas:\t\t  ", correctYea, " of ", allYea, " (", round(100*correctYea/allYea,1), "%) predictions correct", sep="")
    cat("\nPredicted Nays:\t\t  ", correctNay, " of ", allNay, " (", round(100*correctNay/allNay,1), "%) predictions correct", sep="")
    cat("\nCorrect Classifiction:\t ", paste(round(object$fits[1:object$dimensions],2),"%",sep=""), sep=" ")
    cat("\nAPRE:\t\t\t ", round(object$fits[(object$dimensions+1):(2*object$dimensions)],3), sep=" ")
    cat("\nGMP:\t\t\t ", round(object$fits[(2*object$dimensions+1):(3*object$dimensions)],3), "\n\n\n", sep=" ")


                
    if(!verbose) {
        cat("The first 10 legislator estimates are:\n")
	if(object$dimensions!=1) {
	round(object$legislators[1:10,paste("coord",1:object$dimensions,"D",sep="")],3)
	} else{
	round(object$legislators[1:10,c("coord1D","se1D")],3)
	}
    } else {
	if(object$dimensions!=1) {
        round(object$legislators[,paste("coord",1:object$dimensions,"D",sep="")],3)
	} else{
	round(object$legislators[,c("coord1D","se1D")],3)
	}
    }


}
