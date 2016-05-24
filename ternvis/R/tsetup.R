tsetup <-
function(L=diag(c(1,1,1))/sqrt(2), # use Brier score by default
                   q=cbind(1,1,1)/3          # use tercile climatology by default
)
{
oB    <- cbind(1,0,0)
oN    <- cbind(0,1,0)
oA    <- cbind(0,0,1)
#
# linear algebra rules to
# map column vector p <-> vector x living on simplex in R^3
#
# define lengths of sides of triange in R2
    a <- as.double(sqrt((oB-oN)%*%t(L)%*%L%*%t(oB-oN)))
    b <- as.double(sqrt((oA-oN)%*%t(L)%*%L%*%t(oA-oN)))
    n <- as.double(sqrt((oB-oA)%*%t(L)%*%L%*%t(oB-oA)))
  phi <- acos((a^2+n^2-b^2)/(2*a*n))


M32 <- rbind(cbind(0 ,a*cos(phi),n),
             cbind(0 ,a*sin(phi),0)   )

M23 <- rbind(cbind(-a*sin(phi) ,  a*cos(phi)-n ),
             cbind(          0 ,  n            ),
             cbind( a*sin(phi) , -a*cos(phi))  )   / (a*n*sin(phi))

colnames(M23) <- NULL	   
colnames(M32) <- NULL	

    xmin <- min(0,a*cos(phi))
    xmax <- max(n,a*cos(phi))
    ymin <- 0
    ymax <- a*sin(phi)	
        
    height <- ymax - ymin
    width  <- xmax - xmin
    maxhw  <- max(height,width) 
    extra  <- 0.03*maxhw
    
    xlim <- c(xmin-extra, xmin + maxhw+extra)
    ylim <- c(ymin-extra, ymin + maxhw+extra)
    
    xoff <- 0
    yoff <- 0.5*(maxhw-height)

Q  <- M32 %*% t(q)
OB <- M32 %*% t(oB)
ON <- M32 %*% t(oN)
OA <- M32 %*% t(oA)

out <- list(
       a=a,
	   b=b,
	   n=n,
	   phi=phi,
	   M32=M32,
	   M23=M23,
	   xmin=xmin,
	   ymin=ymin,
	   xmax=xmax,
	   ymax=ymax,
	   height=height,
	   width=width,
	   maxhw=maxhw,
	   extra=extra,
	   xlim=xlim,
	   ylim=ylim,
	   xoff=xoff,
	   yoff=yoff,
       oB=oB,
       oN=oN,
       oA=oA,      
       q=q,
       OB=OB,
       ON=ON,
       OA=OA,
       Q=Q      
           )
out
}
