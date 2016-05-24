

### intersphere ###


#' @export
#' @name intersphere
#' @aliases LSD.intersphere
#' @title Intersphere: a fancy Venn diagram
#' @description Create circles for visualizing overlaps between up to 4 datasets.
#' @param data a list with n entries having elements that can be represented as sets (have union and intersect methods).
#' @param colors a character vector of  R build-in colors for circles representing different sets.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param cex a numeric value giving the character expansion factor for intersect size text inside each circle.
#' @param expand.circles a numeric value giving the expansion factor of circles (multiplicative).
#' @param expand.lims a numeric value giving the expansion of x and y limits (additive).
#' @param main title(s) of the plot, standard graphics parameter.
#' @param onlySets vectors, which n-overlaps should be shown, default to all 1 < n < length(data).
#' @author Sebastian Duemcke, Bjoern Schwalb
#' @seealso \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples data = list(
#' 	"A" = sample(1:200,100),
#' 	"B" = sample(1:200,150),
#' 	"C" = sample(1:200,50))
#' 
#' intersphere(data,colors = c("orange","skyblue","green"))
#' 
#' data = list(
#' 	"A" = sample(1:200,100),
#' 	"B" = sample(1:200,150),
#' 	"C" = sample(1:200,50),
#' 	"D" = sample(1:200,75))
#' 
#' colors = c("orange","skyblue","green","purple")
#' intersphere(data,colors,expand.circles = 0.5,expand.lims = 0.5)
#' @keywords intersphere, Venn


intersphere = function(data,colors = NULL,alpha = 25,cex = 1,expand.circles = 1,expand.lims = 1.5,main = "intersphere: overlap diagram",onlySets = seq(length(data),2,by=-1))
{
    
    # number of cicles and overlap sizes #
    
    nb.circles = 2^length(data) - 1
    get.n.overlap = function(n,data){combn(1:length(data),n,FUN=function(x,d){length(Reduce(intersect,d[x]))},simplify=TRUE,data)}
    get.overlap.size = function(d,...){length(Reduce(intersect,d[...]))}
    circle.sizes = unlist(lapply(1:length(data),get.n.overlap,data))
    stopifnot(length(circle.sizes) == nb.circles)
    
    # scale of the radii #
    
    radius.scale = max(circle.sizes)/expand.circles
    radius.scale[radius.scale == 0] = 1
    
    # default colors #
    
    if(is.null(colors)){colors = distinctcolors(nb.circles)}
    
    # get coordinates of points on unit circle (original sets) #
    
    xcoor = cos(seq(0, 2 * pi , length = length(data)+1)+pi/4)/0.5
    ycoor = sin(seq(0, 2 * pi , length = length(data)+1)+pi/4)/0.5
    coords = cbind(xcoor,ycoor)
    
    # define plotting region #
    
    xmax = max(xcoor)
    ymax = max(ycoor)
    plot(0,type="n",axes=FALSE,xlim=c(-xmax-expand.lims,xmax+expand.lims),ylim=c(-ymax-expand.lims,ymax+expand.lims),ylab="",xlab="",main="")
    mtext(paste(main),3,2,cex=1.25)
    
	# draw.circles function #
	
	draw.circles = function(x,y,radius,border = "black",col = "white"){
		roots = seq(0,2*pi - 2*pi/100,by = 2*pi/100)
		for (circle in 1:length(radius)){
			xv = cos(roots)*radius[circle] + x[circle]
			yv = sin(roots)*radius[circle] + y[circle]
			polygon(xv,yv,border = border,col = col[circle])
		}
	}
	
    # draw function #
    
    draw = function(points,coords,data){
        mass = apply(coords[points,],2,sum) / length(points)
        # shift mass point. new mass = third of way between center and pairwise intersect of neighbouring #
        if(length(points) == 2 && abs(points[1] - points[2]) == 2 && length(data) != 3){
            mass = ((apply(coords[1:length(data),],2,sum) / length(data)) + (apply(coords[c(points[1],points[1]+1),],2,sum) / 2)) * 0.8
        }
        # draw line segments to mass point in same color as origin #
        sapply(points,function(p){points(rbind(coords[p,],mass),col=colors[p],type="l")})
        # draw mass points #
        draw.circles(mass[1],mass[2],r=sqrt(get.overlap.size(data,points)/(pi*radius.scale)),col="white",border="black")
        draw.circles(mass[1],mass[2],r=sqrt(get.overlap.size(data,points)/(pi*radius.scale)),col=convertcolor("grey50",alpha),border="grey50")
        # add text #
        text(mass[1],mass[2],labels=get.overlap.size(data,points),cex=cex)
    }
    
    # draw mass points, line segments and text #
    
    sapply(onlySets,function(n,draw,coords,data){combn(1:length(data),n,FUN=draw,simplify=FALSE,coords,data)},draw,coords,data)
    sapply(1:length(data),function(p){draw.circles(coords[p,1],coords[p,2],r=sqrt(length(data[[p]])/(pi*radius.scale)),col="white",border="black")})
    sapply(1:length(data),function(p){draw.circles(coords[p,1],coords[p,2],r=sqrt(length(data[[p]])/(pi*radius.scale)),col=convertcolor(colors[p],alpha),border=colors[p])})
    sapply(1:length(data),function(p){text(coords[p,1],coords[p,2],labels=length(data[[p]]),cex=cex)})
    
    # return invisible #
    
    invisible()
}


### aliases ###


LSD.intersphere = intersphere



