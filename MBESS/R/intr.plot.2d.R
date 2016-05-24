"intr.plot.2d" <-
function(b.0, b.x, b.z, b.xz,x.min=NULL, x.max=NULL, x=NULL, 
n.x=50, mean.z=NULL, sd.z=NULL, z=NULL,xlab="Value of X",  
ylab="Dependent Variable", sd.plot=TRUE, sd2.plot=TRUE, sd_1.plot=TRUE, 
sd_2.plot=TRUE, type.sd=2, type.sd2=3, type.sd_1=4, type.sd_2=5, 
legend.pos="bottomright", legend.on=TRUE, ... )
{

## situation X1, to input a specific X vector
if(!is.null(x))
{
        if(!is.null(x.min) || !is.null(x.max) )
        stop("Only the vector x is needed. Please do not enter x's range.")
}   



## situation X2, to input X's range
if(!is.null(x.min) ) {
        if (is.null(x.max)) stop("Since you entered x.min, please enter x.max too.") }
        
if(!is.null(x.max) ) {
        if (is.null(x.min)) stop("Since you entered x.max, please enter x.min too.") }  

if (!is.null(x.min) && !is.null(x.max) )
        {
        if ( !is.null(x) )
        stop("Since you entered the range of x, Please do not enter the specific vector x.")  
       
        x <- seq(x.min, x.max, length.out=n.x)
        }


## situation Z1, to input mean.z and sd.z
if(!is.null(mean.z) && is.null(sd.z)) stop ("Since you entered mean.z, please enter sd.z too.")
if(!is.null(sd.z) && is.null(mean.z)) stop ("Since you entered sd.z, please enter mean.z too.")

if(!is.null(mean.z) && !is.null(sd.z) ) {
    if(!is.null(z)) stop ("Since you entered mean.z and sd.z, please do not enter the specific z vector.")
    
   
    }


## situation Z2, to input a specific vector z
if(!is.null(z)) {
    if(!is.null(mean.z) || !is.null(sd.z) ) stop ("Only vector z is needed. Please do not enter mean.z 
    and/or sd.z.")
    
    mean.z<- mean(z)
    sd.z<- sd(z)
    }


z0<-mean.z
y0 <- b.0 + b.x*x + b.z*z0 + b.xz*x*z0
plot(x,y0, type="l", xlab = xlab, ylab=ylab, ...)
legend.p<-c("at z's mean")
type.p<-c(1)
            
if(sd.plot)
{
z1 <- z0+sd.z
y1 <- b.0 + b.x*x + b.z*z1 + b.xz*x*z1
lines (x, y1, lty=type.sd)
legend.p<-c(legend.p, "1 sd above z's mean")
type.p<-c(type.p, type.sd)
}
        
if (sd2.plot)
{
z2 <- z0+2*sd.z
y1 <- b.0 + b.x*x + b.z*z2 + b.xz*x*z2
lines (x, y1, lty=type.sd2)
legend.p<-c(legend.p, "2 sd above z's mean")
type.p<-c(type.p, type.sd2)
}
        
if (sd_1.plot)
{
z_1 <- z0-sd.z
y_1 <- b.0 + b.x*x + b.z*z_1 + b.xz*x*z_1
lines (x, y_1, lty=type.sd_1)
legend.p<-c(legend.p, "1 sd below z's mean")
type.p<-c(type.p, type.sd_1)
}
        
if (sd_2.plot)
{
z_2 <- z0-2*sd.z
y_2 <- b.0 + b.x*x + b.z*z_2 + b.xz*x*z_2
lines (x, y_2, lty=type.sd_2)
legend.p<-c(legend.p, "2 sd below z's mean")
type.p<-c(type.p, type.sd_2)
}

if(legend.on)    
{legend(legend.pos, legend=legend.p, lty=type.p)}



}

