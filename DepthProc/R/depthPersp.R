#' @name depthPersp
#' @title Perspective plot for depth functions
#' @export
#' @importFrom lattice wireframe
#' @importFrom colorspace heat_hcl
#' @description Draws a perspective plot of depth function over x-y plane.
#' 
#' @param x bivariate data
#' @param plot_method there are two options "lattice", and "rgl" - see details
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#' @param n number of points that will be used to create plot (n^2)
#' @param xlab description of x-axis
#' @param ylab description of y-axis
#' @param plot_title plot title (default NULL means paste(method, "depth"))
#' @param colors function for colors pallete (e.g. gray.colors).
#' @param ... arguments passed to depth function
#' 
#' @details
#' 
#' plot_method - rgl package is not in depends list beacuse it may cause problems when OpenGL is not supported.  To use plot_method = "rgl" you must load this package on your own. 
#' 
#' @author Daniel Kosiorowski, Mateusz Bocian, Anna Wegrzynkiewicz and Zygmunt Zawadzki from Cracow University of Economics.
#' 
#' @examples
#'  x = mvrnorm(100,c(0,0),diag(2))
#'  depthPersp(x, method = "Euclidean")
#' 
#' # EXAMPLE 2
#' data(inf.mort,maesles.imm)
#' data1990=na.omit(cbind(inf.mort[,1],maesles.imm[,1]))
#' 
#' \dontrun{
#' require(rgl)
#' depthPersp(data1990, method = "Projection",plot_method= "rgl")
#' }
#' 
depthPersp<-function(x, plot_method = "lattice", xlim = extendrange(x[,1],f=0.1), ylim = extendrange(x[,2],f=0.1),n=50,
	xlab = "x", ylab = "y", plot_title=NULL, colors = heat_hcl, ...)
{
	if(dim(x)[2]==2)
	{
	
			axis_x = seq(xlim[1],xlim[2],length.out = n)
 			axis_y = seq(ylim[1],ylim[2],length.out = n)
				
			xy_surface = expand.grid(axis_x,axis_y)

		
			xy_surface=matrix(unlist(xy_surface),ncol=2)
  
      
			depth_params = .extractDepthParams(xy_surface,x,...)
			z_surface = do.call(depth,depth_params)
			method = depth_params$method
      if(is.null(plot_title)) plot_title = paste(method, "depth")
			graph_params = .removeDepthParams(...)
			#z_surface = depth(xy_surface,x,method = method,...)

			ztmp = z_surface * 100/max(z_surface)  
    
			#colors = rev(rainbow(100,start=0,end=1/4))
			colors = colors(100)
      col = colors[ztmp]
			
			if(plot_method == "rgl") 
        {
          rgl::persp3d(axis_x,axis_y,z_surface,color=col,back = "lines",
				xlab = xlab, ylab = ylab,zlab=plot_title)
			  }
			
			if(plot_method == "lattice") wireframe(z_surface~xy_surface[,1]+xy_surface[,2],colorkey = TRUE,drape = TRUE,
					xlab = xlab, ylab = ylab,zlab="",col.regions = colors,lwd = 0.4,main  = plot_title,scales = list(arrows = FALSE),...)
			else {print=c("Wrong plot.method")}
		}
		else print("Wrong matrix!")
}
