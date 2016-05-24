
# Function: add.cutline
# Reads in output from W-NOMINATE and adds a cutting line to existing plot
#   INPUTS: a numeric vector of length 4, cutData
#       normVector1D<-cutData[1]
#       normVector2D<-cutData[2]
#       midpoints<-cutData[3]

add.OCcutline <- function(cutData,lwd=2) {

    slope <- -cutData[1]/cutData[2]
    if (is.na(slope)) {
        x <- c(cutData[1]*cutData[3],cutData[1]*cutData[3])
            y <- c(sqrt(1-(cutData[1]*cutData[3])^2),-sqrt(1-(cutData[1]*cutData[3])^2))
                slope <- NA
                intercept <- NA
        }
        else {
                intercept <- -slope*cutData[1]*cutData[3]+cutData[2]*cutData[3]
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

plot.OCangles <- function(x, main.title="Cutting Line Angles",
        x.title="Angle in Degrees", y.title="Count",dims=c(1,2),...) {

    if(!class(x)=="OCobject")
        stop("Input is not of class 'OCobject'.")
    if(x$dimensions==1)
        stop("All angles in 1D Optimal Classification are 90 degrees.")
    if(length(dims)!=2)
        stop("'dims' must be an integer vector of length 2.")

    cutvector1 <- na.omit(x$rollcalls[,paste("normVector",dims[2],"D",sep="")])
    cutvector2 <- na.omit(x$rollcalls[,paste("normVector",dims[1],"D",sep="")])
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

plot.OCcutlines <- function(x, main.title="Cutting Lines",
        d1.title="First Dimension", d2.title="Second Dimension",
        lines=50,dims=c(1,2),lwd=2,...) {

    if(!class(x)=="OCobject")
        stop("Input is not of class 'OCobject'.")
    if(x$dimensions==1)
        stop("All angles in 1D Optimal Classification are 90 degrees.")
    if(length(dims)!=2)
        stop("'dims' must be an integer vector of length 2.")
    if(lines<1)  stop("'Lines' must be less than 1.")

#    constrained <- ((abs(x$rollcalls[,"spread1D"]) > 0.0 | abs(x$rollcalls[,"spread2D"]) > 0.0)
#        & (x$rollcalls[,"midpoint1D"]**2 + x$rollcalls[,"midpoint2D"]**2) < .95)
    
    cutlineData <- cbind(x$rollcalls[,paste("normVector",dims[1],"D",sep="")],
                     x$rollcalls[,paste("normVector",dims[2],"D",sep="")],
                     x$rollcalls[,"midpoints"])
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
    
    suppressWarnings(apply(cutlineData, 1, add.OCcutline,lwd=lwd))

}

plot.OCcoords <- function (x, main.title="OC Coordinates",
    d1.title="First Dimension", d2.title="Second Dimension", dims=c(1,2),
    plotBy="party", color=TRUE, shape=TRUE, cutline=NULL, Legend=TRUE,
    legend.x=0.8,legend.y=1,...) {
   
    if(!class(x)=="OCobject")
        stop("Input is not of class 'OCobject'.")
    if(!any(colnames(x$legislators)==plotBy)){
        warning("Variable '", plotBy ,"' does not exist in your OC object.")
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
        ranking <- x$legislators[,"coord1D"]
        plot(1:length(coord1D),
                1:length(coord1D),
                type="n",
                cex.main=1.2,
                cex.lab=1.2,
                font.main=2,
                xlab="Rank",
                ylab="Rank",
                main="1D Optimal Classification Plot")

        if(Legend)  legend(0.9*length(coord1D), 0.7*length(coord1D), unique(types),
                            pch=shapes[1:nparties], col=colorlist[1:nparties], cex=0.7)
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
        add.OCcutline(c(x$rollcalls[cutline[i],paste("normVector",dims[1],"D",sep="")],
                    x$rollcalls[cutline[i],paste("normVector",dims[2],"D",sep="")],
                    x$rollcalls[cutline[i],"midpoints"]),lwd=2)
        }
    }
   
    if(Legend)
        legend(legend.x, legend.y, unique(types), pch=shapes[1:nparties],
    bty="n",col=colorlist[1:nparties], cex=0.7)

    for(i in 1:nparties) suppressWarnings(points(coord1D[types==unique(types)[i]],
        coord2D[types==unique(types)[i]],pch=shapes[i],col=colorlist[i],cex=1.1,lwd=2))
    }
}

plot.OCskree <- function(x, main.title="Skree Plot", x.title="Dimension",
                        y.title="Eigenvalue",...) {

    if(!class(x)=="OCobject")
        stop("Input is not of class 'OCobject'.")
    if(is.null(x$eigenvalues))
    stop("No eigenvalues exist in this OC object.")
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

plot.OCobject <- function(x,dims=c(1,2),...) {
    if(!class(x)=="OCobject")
        stop("Input is not of class 'OCobject'.")
    if(length(dims)!=2 & x$dimensions!=1)
        stop("'dims' must be an integer vector of length 2.")

    if(x$dimensions==1) {
        par(mfrow=c(1,2))       
        suppressWarnings(plot.OCcoords(x,dims=dims))
        suppressWarnings(plot.OCskree(x,dims=dims))
    } else {
        par(mfrow=c(2,2))
        suppressWarnings(plot.OCcoords(x,dims=dims))
        suppressWarnings(plot.OCangles(x,dims=dims))
        suppressWarnings(plot.OCskree(x,dims=dims))
        suppressWarnings(plot.OCcutlines(x,dims=dims,lwd=1))
    }
}               

summary.OCobject<-function(object,verbose=FALSE,...){

    if(!class(object)=="OCobject")
        stop("Input is not of class 'OCobject'.")

    cat("\n\nSUMMARY OF OPTIMAL CLASSIFICATION OBJECT")
    cat("\n----------------------------\n")
    cat("\nNumber of Legislators:\t ", dim(na.omit(object$legislators))[1],
    " (", dim(object$legislators)[1]-dim(na.omit(object$legislators))[1],
    " legislators deleted)", sep="")
    cat("\nNumber of Votes:\t ", dim(na.omit(object$rollcalls))[1],
    " (", dim(object$rollcalls)[1]-dim(na.omit(object$rollcalls))[1],
    " votes deleted)", sep="")
    cat("\nNumber of Dimensions:\t", object$dimensions)

    correctYea<-sum(as.numeric(object$legislators[,"correctYea"]),na.rm=TRUE)
    allYea<-correctYea+sum(as.numeric(object$legislators[,"wrongNay"]),na.rm=TRUE)
    correctNay<-sum(as.numeric(object$legislators[,"correctNay"]),na.rm=TRUE)
    allNay<-correctNay+sum(as.numeric(object$legislators[,"wrongYea"]),na.rm=TRUE)

    cat("\nPredicted Yeas:\t\t ", correctYea, " of ", allYea, " (", round(100*correctYea/allYea,1), "%) predictions correct", sep="")
    cat("\nPredicted Nays:\t\t ", correctNay, " of ", allNay, " (", round(100*correctNay/allNay,1), "%) predictions correct\n\n", sep="")
                
    if(!verbose) {
        cat("The first 10 legislator estimates are:\n")
    if(object$dimensions!=1) {
    round(object$legislators[1:10,paste("coord",1:object$dimensions,"D",sep="")],3)
    } else{
    round(object$legislators[1:10,c("coord1D"),drop=FALSE],3)
    }
    } else {
    if(object$dimensions!=1) {
        round(object$legislators[,paste("coord",1:object$dimensions,"D",sep="")],3)
    } else{
    round(object$legislators[,c("coord1D"),drop=FALSE],3)
    }
    }


}
