"intr.plot" <-
function(b.0, b.x, b.z, b.xz, x.min=NULL, x.max=NULL, z.min=NULL, z.max=NULL, n.x=50, n.z=50, 
x=NULL, z=NULL, col="lightblue", hor.angle=-60, vert.angle=15,
xlab="Value of X", zlab="Value of Z", ylab="Dependent Variable", expand=0.5, lines.plot=TRUE, 
col.line="red", line.wd=2, gray.scale=FALSE, ticktype="detailed",  ...)
{


#########  Situation 1
if (!is.null(x.min) && !is.null(x.max))
{
if(x.min>x.max)
stop ("'x.max' must be bigger than 'x.min'!")

if(!is.null(x) ) stop ("Since you entered the range of x, please do not enter the specific vector x.")
} 

if(!is.null(z.min) && !is.null(z.max))
{
if(z.min>z.max)
stop ("'z.max' must be bigger than 'z.min'!")

if(!is.null(z) ) stop ("Since you entered the range of z, please do not enter the specific vector z.")
}

###################################
    
###########  Situation 2
if (!is.null(x)) 
{if(!is.null(x.min)|| !is.null(x.max) ) 
stop("Since you entered the sepcific vector x, please do not enter the range of x.")  }


if( !is.null(z) )
        {
        if( !is.null(z.min) || !is.null(z.max) )
        stop("Since you entered the sepcific vector z, please do not enter the range of z.")      
        }

#######################################

if(!is.null(x.min) && is.null (x.max) ) stop ("Since you entered 'x.min', 'x.max' is also needed.")
if(!is.null(x.max) && is.null (x.min) ) stop ("Since you entered 'x.max', 'x.min' is also needed.")
if(!is.null(z.min) && is.null (z.max) ) stop ("Since you entered 'z.min', 'z.max' is also needed.")
if(!is.null(z.max) && is.null (z.min) ) stop ("Since you entered 'z.max', 'z.min' is also needed.")


#######  Generate vectors
if(is.null(x) )
{ x <- seq(x.min, x.max, length.out=n.x) }

if(is.null(z) )
{ z <- seq(z.min, z.max, length.out=n.z) }

############################################

if(gray.scale)
{col="white"
col.line="black"
}


    Pred.Eq <- function(x, z, b0=b.0, b1=b.x, b2=b.z, b3=b.xz) 
        {
        return(b0 + b1*x + b2*z + b3*x*z)
        }
               
    y.plot <- outer(x, z, Pred.Eq)

    persp(x, z, y.plot, theta=hor.angle, phi=vert.angle, expand=expand, col=col, 
    xlab = xlab, ylab = zlab, zlab = ylab, ticktype=ticktype, ...)-> surface


if(lines.plot){
        z0<- mean(z)
        sd.z<- sd(z)
        z1<- z0+sd.z
        z2<- z0+2*sd.z
        z_1<- z0-sd.z
        z_2<- z0-2*sd.z
        lines.legend<- paste("Plotted regression lines are \n")
        
        if( z_2>min(z) && z_2<max(z) ){
        lines(trans3d(x, y=z_2, z=b.0+b.x*x+b.z*z_2+b.xz*x*z_2, surface), col=col.line,lwd=line.wd)
        lines.legend<- paste(lines.legend, "-2, " )
        
        }
        
        if(z_1>min(z) && z_1<max(z) ){
        lines(trans3d(x, y=z_1, z=b.0+b.x*x+b.z*z_1+b.xz*x*z_1, surface), col=col.line,lwd=line.wd)
         lines.legend<- paste(lines.legend, "-1, " )
        }
           
   
        if(z1>min(z) && z1<max(z) ){
        lines(trans3d(x, y=z1, z=b.0+b.x*x+b.z*z1+b.xz*x*z1, surface), col=col.line,lwd=line.wd)
         lines.legend<- paste(lines.legend, "1, " )
        }
        
        if(z2>min(z) && z2<max(z) ){
        lines(trans3d(x, y=z2, z=b.0+b.x*x+b.z*z2+b.xz*x*z2,surface), col=col.line,lwd=line.wd)
        lines.legend<- paste(lines.legend, "2, " )
        }
        
        lines(trans3d(x, y=z0, z=b.0+b.x*x+b.z*z0+b.xz*x*z0, surface), col=col.line,lwd=line.wd)
        lines.legend<- paste(lines.legend, "and 0 \n" )
        
        lines.legend<- paste(lines.legend, "standard deviations above z's mean." )

}

angle<- paste("horizonal angle=", hor.angle,";", "vertical angle=", vert.angle)

title(main=angle, sub=lines.legend)
}
