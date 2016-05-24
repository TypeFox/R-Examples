tpalette <-
function(q=cbind(1,1,1)/3,
                     bars = TRUE,
                     m  = 0.7, 
                     flip  = FALSE,
                     lmain = FALSE,
                     theta0 = 0,
                     dich = "none",
                     cex =1,
		             nhist = 12)
{
xoff <- tsetup()$xoff
yoff <- tsetup()$yoff
width <- tsetup()$width
height <- tsetup()$height

if (bars) {dat <- tgrid(nhist)}
    else  {nn <- 80; dat <- tgrid(nn)}

B <- dat[,1]
N <- dat[,2]
A <- dat[,3]

x <- xf(dat)
y <- yf(dat)

if (lmain)
{
main <- substitute("Palette when "*a*" = "*b*""*e*" , "*c*" = "*d*"",
                  list(a=parse(text="theta[0]")[[1]],
		       b=parse(text=sprintf("%.2f",theta0))[[1]],
		       c=parse(text="m")[[1]],
		       d=parse(text=sprintf("%.2f",m))[[1]],
                       e=parse(text="degree")[[1]]
		      )
		 )
	      
qstr <- paste("(",
              sprintf("%.2f",q[1]),
	      ", ",
	      sprintf("%.2f",q[2]),
	      ", ",
	      sprintf("%.2f",q[3]),
	      ")",
	      sep="")

sub <- bquote(bold(q) * minute == .(qstr) )
}
else
{
main <- ""
sub <-  ""
}
tplot(q,L = diag(c(1,1,1))/sqrt(2),col="transparent",cex=0,main=main,sub=sub,newpage=FALSE,grid=FALSE)
if (bars)
{
    pushViewport(viewport(width  = 0.8, 
                          height = 0.8, 
			              xscale = tsetup()$xlim,
                          yscale = tsetup()$ylim))


for (i in 1:length(B))
{
pushViewport(viewport(x=xoff + x[i],
                      y=yoff + y[i],
		      just=c("centre","centre"),
		      default.units="native",
		      width=width/nhist,
		      height=height/nhist       
		      )
	     )                       # create viewport for legend

   grid.polygon(x=c(0,1,1,0),
                y=c(0,0,1,1),
	        gp=gpar(col="black",
	                fill=tcolour(p=rbind(dat[i,]),
			             q=q,
				     m=m,
				     flip=flip,
				     theta0=theta0,
				     dich=dich)
		       )
                )
                
   grid.polygon(x=c(0,1/3,1/3,0),
                y=c(0,0,B[i],B[i]),
	        gp=gpar(col="black",fill="grey"))
   grid.polygon(x=c(1/3,2/3,2/3,1/3),
                y=c(0,0,N[i],N[i]),
	        gp=gpar(col="black",fill="grey"))
   grid.polygon(x=c(2/3,1,1,2/3),
                y=c(0,0,A[i],A[i]),
	        gp=gpar(col="black",fill="grey"))
   grid.polygon(x=c(0,1,1,0),
                y=c(0,0,1,1))	     
	     
popViewport()
} # end of i loop
popViewport()
} # if (bars)
else
{
   col <- matrix(0,nrow=nn,ncol=1)  # initialise vector of colours
   for (i in 1:length(B))
   {	 
      col[i] = tcolour(p=cbind(B[i],N[i],A[i]),
                       q=q,
		       m=m, 
		       flip = flip,
		       theta0 = theta0,
                       dich = dich)
   }
tplot( cbind(B,N,A),
         L = diag(c(1,1,1))/sqrt(2),
         dimnames=NULL,
         col=col, 
         pch=19,  
         cex=cex,
         main=main,
         sub=sub,
         newpage=FALSE,
         bg=FALSE,
         grid=FALSE,
         border=FALSE
       )   
} # end of if (bars) 
# now overlay climatology
tplot(q,
    L = diag(c(1,1,1))/sqrt(2),
	col="blue",
	pch=4,
	border="transparent",
	main="",
	bg=FALSE,
	grid=FALSE,
	newpage=FALSE)
}
