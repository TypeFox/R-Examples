idclass <- function(data,type=c("PKa","SBC","KHa","KH","PK"),a.in=NULL,
                    outplot=c("summary","detail","none"),plot.focus=NULL){
# Time series categorisation for intermittent demand
#
# Inputs
#   data          Time series dataset. Each column is a series.
#                 Alternatively this can be a single series.
#   type          Type of categorisation:
#                   "SBC" - Syntetos Boylan Croston;
#                   "KH"  - Kostenko Hyndman (exact*);
#                   "KHa" - Kostenko Hyndman (approximate);
#                   "PK"  - Petropoulos Kourentzes (exact*);
#                   "PKa" - Petropoulos Kourentzes (approximate).
#                 * These are computationally expensive, as SBA is optimised
#                   for each time series. 
#   a.in          Vector of SBA demand interval smoothing parameters. This
#                 must be same length as number of series. This is used for
#                 categorisations "KH" and "PK". If a.in == NULL then the
#                 parameters are calculated internally using MAR as a cost
#                 function.
#   outplot       Plot results of categorisation:
#                   "summary" - simlified plot that reports number of series
#                               in each class and cut-off points;
#                   "detail"  - scatterplot between average inter-demand 
#                               interval (p) and squared coefficient of
#                               variation of non-zero demand (CV^2). Series
#                               that are categorised for SBA or SES are plotted
#                               in shaded areas;
#                   "none"    - do not produce plot.
#   plot.focus    Only relevant to outplot == "detail". Can be used to specify
#                 the maximum p and CV^2 to plot, so that the scatterplot can 
#                 be focused on the separation area between the categories. Use
#                 vector of two elements. First one is max p and second one is 
#                 max CV^2. Example: plot.focus=c(1.5,1.5). If NULL then maximums
#                 are defined from the dataset. 
#
# Outputs
#   idx.croston   Index of series that are categorised under Croston.
#   idx.sba       Index of series that are categorised under SBA.
#   idx.ses       Index of series that are categorised under SES. Provided only 
#                 for "PK" and "PKa" types.
#   cv2           Coefficient of variation squared of non-zero demand.
#   p             Inter-demand interval. 
#   summary       Summary of number of series under each category. 
#
# Example:
#   # Create/load some data. Each column is a time series
#   dataset <- simID(100,60,idi=1.15,cv2=0.3)
#   idclass(dataset)
#
# Notes:
# Classification schemes described in:
# F. Petropoulos and N. Kourentzes, 2014, Journal of Operational Research Society. 
# http://dx.doi.org/10.1057/jors.2014.62
# http://kourentzes.com/forecasting/2014/05/13/forecast-combinations-for-intermittent-demand/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

type <- type[1]
type <- toupper(type)
outplot <- outplot[1]

# Initialise
if ((sum(dim(data)==1)>0 | is.vector(data)) && class(data)!="data.frame"){
  data <- matrix(data,nrow=length(data))
  N <- 1
  p <- NA
  v <- NA
  if (is.null(a.in)){
    a <- NA
  } else {
    a <- a.in[1]
  }
} else {
  N <- dim(data)[2]
  p <- vector("numeric",N)    # Average interval
  v <- vector("numeric",N)    # Squared CV
  if (is.null(a.in)){
    a <- vector("numeric",N)  # Interval alpha with SBA
  } else {
    a <- a.in
    if (length(a)!=N){
      stop("Length of a.in must be equal to number of series.")
    }
  }
}

for (i in 1:N){
  # Get time series
  temp <- data[,i]
  temp <- temp[!is.na(temp)]
  # Croston decomposition
  nzd <- which(temp!=0)
  k <- length(nzd)
  z <- temp[nzd]                        # Demand
  x <- c(nzd[1],nzd[2:k]-nzd[1:(k-1)])  # Intervals
  p[i] <- mean(x)
  v[i] <- (sd(z)/mean(z))^2
  if (is.null(a.in)){
    if (type=="KH" | type=="PK"){
      fit <- crost(temp,h=0,type="sba")
      a[i] <- fit$weights[2]
    }
  }
}

# Classify between SBA & Croston
if (type=="KH" | type=="PK"){
  # Exact Konstenko & Hyndman
  use.sba <- v > (4*p*(2-p)-a*(4-a)-p*(p-1)*(4-a)*(2-a))/(p*(4-a)*(2*p-a))
} else if(type=="SBC"){
  # Syntetos, Boylan & Croston
  use.sba <- (p > 1.32 | v > 0.49)
} else {
  # Approximate Konstenko & Hyndman
  use.sba <- v > (2-(3/2)*p)
}
use.croston <- !use.sba

# Classify between SES and Croston/SBA
if (type=="PK" | type=="PKA"){
  use.ses <- p <= 1
  use.sba[use.ses] <- FALSE
  use.croston[use.ses] <- FALSE
  idx.ses <- which(use.ses)
} else {
  use.ses <- NULL
  idx.ses <- NULL
}
idx.croston <- which(use.croston)
idx.sba <- which(use.sba)

# Summarise results
if (is.null(use.ses)){
  summary <- rbind(sum(use.croston),sum(use.sba))
  rownames(summary) <- c("Croston","SBA")
} else {
  summary <- rbind(sum(use.croston),sum(use.sba),sum(use.ses))
  rownames(summary) <- c("Croston","SBA","SES")
}
colnames(summary) <- "Series"

if (type=="PKA"){
  type.plot = "PKa"
} else if(type=="KHA"){
  type.plot = "KHa"
} else {
  type.plot = type
}

# Produce detailed plots
if (outplot == "detail"){
  # Set space in plot for PK and PKa types
  if (type=="PK" | type=="PKA"){
    if (!is.null(plot.focus)){
      xmin <- 1-(plot.focus[1]-1)*0.1
    } else {
      xmin <- 1-(max(c(max(p)+diff(range(p))*0.1),1.5)-1)*0.1
    }
  } else {
    xmin <- 1
  }
  # Set limits depending if plot.focus is given or not
  if (!is.null(plot.focus)){
    xx <- c(xmin,plot.focus[1])
    yy <- c(0,plot.focus[2])
    plot.out <- (p > plot.focus[1] | v > plot.focus[2])
  } else {
    xx <- c(xmin,max(c(max(p)+diff(range(p))*0.1),1.5))
    yy <- c(0,max(c(max(v)+diff(range(v))*0.1),1.5))
    plot.out <- vector("logical",N)
  }
  # Plot SBA series in focus
  plot(p[use.sba & !plot.out],v[use.sba & !plot.out],
       type="p",pch=20,xaxs="i",yaxs="i",
       xlim=xx,ylim=yy,cex=0.8,xlab="p",ylab=parse(text="CV^2"),
       main=paste(type.plot,"classification"))
  # Plot SBA series out-of-focus
  if (sum(plot.out)>1){
    p.temp <- p[plot.out]
    p.temp[p.temp>xx[2]] <- xx[2]
    v.temp <- v[plot.out]
    v.temp[v.temp>yy[2]] <- yy[2]
    points(p.temp,v.temp,pch=21,bg="darkgray",cex=1)
  }
  # Create shaded area for Croston series
  if (type == "SBC"){
    polygon(c(1,1.32,1.32,1),c(0,0,0.49,0.49),col="lightgrey",border=NA)
  } else if(type == "KHA" | type == "PKA"){
    polygon(c(1,4/3,1),c(0,0,0.5),col="lightgrey",border=NA)
  } else if(type == "KH" | type == "PK"){
    p.x <- c(seq(1,4/3,0.02),4/3)
    a.p <- 1
    v.y <- (4*p.x*(2-p.x)-a.p*(4-a.p)-p.x*(p.x-1)*(4-a.p)*(2-a.p))/(p.x*(4-a.p)*(2*p.x-a.p))
    polygon(c(p.x,1),c(v.y,0),col="lightgrey",border=NA)
    a.p <- 0
    v.y2 <- (4*p.x*(2-p.x)-a.p*(4-a.p)-p.x*(p.x-1)*(4-a.p)*(2-a.p))/(p.x*(4-a.p)*(2*p.x-a.p))
    v.y2[v.y2<0] <- 0
    polygon(c(p.x,p.x[18:1]),c(v.y,v.y2[18:1]),col=gray(0.2,0.3),border=NA)
  }
  # Plot Croston series
  points(p[use.croston],v[use.croston],pch=1,col="red",cex=0.6)
  # Create area for ses
  if (type == "PKA" | type == "PK"){
    polygon(c(0,1,1,0),c(0,0,yy[2],yy[2]),col=gray(0.4),border=NA)
    points(rep(xmin+(1-xmin)/2,sum(use.ses & !plot.out)),v[use.ses & !plot.out],
           col="white",pch=20,cex=0.8)
    if (sum(use.ses & plot.out)>1){
      points(xmin+(1-xmin)/2,yy[2],col="grey",pch=21,bg="white")
    }
  }
} 
# Produce simplified plots
if (outplot=="summary"){
  # Set plot area
  if (type=="PK" | type=="PKA"){
    xx <- c(0.94,1.5)
  } else if(type=="SBC"){
    xx <- c(1,1.64)
  } else {
    xx <- c(1,1.5)
  }
  yy <- c(0,1)
  # Start a plot
  plot(0,0,xlim=xx,ylim=yy,xaxs="i",yaxs="i",xlab="p",ylab=parse(text="CV^2"),
       main=paste(type.plot,"classification"),xaxt="n",yaxt="n")
  # Create shaded area for Croston series
  if (type == "SBC"){
    polygon(c(1,1.32,1.32,1),c(0,0,0.49,0.49),col="lightgrey",border=NA)
  } else if(type == "KHA" | type == "PKA"){
    polygon(c(1,4/3,1),c(0,0,0.5),col="lightgrey",border=NA)
  } else if(type == "KH" | type == "PK"){
    p.x <- c(seq(1,4/3,0.02),4/3)
    a.p <- 1
    v.y <- (4*p.x*(2-p.x)-a.p*(4-a.p)-p.x*(p.x-1)*(4-a.p)*(2-a.p))/(p.x*(4-a.p)*(2*p.x-a.p))
    polygon(c(p.x,1),c(v.y,0),col="lightgrey",border=NA)
    a.p <- 0
    v.y2 <- (4*p.x*(2-p.x)-a.p*(4-a.p)-p.x*(p.x-1)*(4-a.p)*(2-a.p))/(p.x*(4-a.p)*(2*p.x-a.p))
    v.y2[v.y2<0] <- 0
    polygon(c(p.x,p.x[18:1]),c(v.y,v.y2[18:1]),col=gray(0.2,0.3),border=NA)
  }
  # Label ticks
  if (type != "SBC"){
    axis(1,at=c(1,4/3),labels=c("1","4/3"))
    axis(2,at=c(0,0.5),labels=c("0","1/2"))
  } else {
    # Special plot for SBC
    axis(1,at=c(1,1.32),labels=c("1","1.32"))
    axis(2,at=c(0,0.49),labels=c("0","0.49"))
    lines(c(1,xx[2]),c(0.49,0.49),lty=3,col="grey")
    lines(c(1.32,1.32),c(0,yy[2]),lty=3,col="grey")
    text(1.16,0.75,"Erratic",col=gray(0.4))
    text(1.16,0.25,"Smooth",col=gray(0.4))
    text(1.48,0.75,"Lumpy",col=gray(0.4))
    text(1.48,0.25,"Intermittent",col=gray(0.4))
  }
  # Report number of series in each area
  text(xx[2],1,paste("SBA: ",sum(use.sba),sep=""),adj=c(1.2,1.2))
  text(1,0,paste("Croston: ",sum(use.croston),sep=""),adj=c(-0.2,-1.2))
  # Shaded area and text for PK and PKa classifications
  if (type == "PKA" | type == "PK"){
    polygon(c(0,1,1,0),c(0,0,yy[2],yy[2]),col=gray(0.6),border=NA)
    text(0.97,0,paste("SES: ",sum(use.ses),sep=""),adj=c(-0.2,0.2),srt=90)
    lines(c(0.94,1),c(0.5,0.5),lty=3)
  }
}

return(list(idx.croston=idx.croston,idx.sba=idx.sba,idx.ses=idx.ses,
            cv2=v,p=p,summary=summary))

}