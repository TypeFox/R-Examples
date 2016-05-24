"plot.diffmich" <- 
function(x, y, which=3, two=2, ...) {
	perm <- x
   par(mar=c(5,5,4,2)) 
   if(which==3) {
       switch(two, nada=0 , disp=par(mfrow=c(1,2)), disp=par(mfrow=c(2,1)))
   }
   plota <- function(perm){
   	with(perm, {
   	xlaa <- mean(range(c(diffa, permsa)))-(diff(range(c(diffa, permsa)))+diff(range(c(diffa, permsa)))*0.2)/2
   	xlba <- mean(range(c(diffa, permsa)))+(diff(range(c(diffa, permsa)))+diff(range(c(diffa, permsa)))*0.2)/2
   	fig.dist <- hist(permsa, xlim=(c(xlaa, xlba)), 
                  main="Is difference in a significant?")
   	abline(v=diffa) 
   	text(diffa, diff(range(fig.dist$counts))/2, adj = -0.5, expression(bold(a0)), cex=1.5 )})
   }
   plotb <- function(perm){
   	with(perm, {    
   	xlab <- mean(range(c(diffb, permsb)))-(diff(range(c(diffb, permsb)))+diff(range(c(diffb, permsb)))*0.2)/2
   	xlbb <- mean(range(c(diffb, permsb)))+(diff(range(c(diffb, permsb)))+diff(range(c(diffb, permsb)))*0.2)/2
   	fig.distF <- hist(permsb, xlim=(c(xlab, xlbb)), 
                  main="Is difference in b significant?")
   	abline(v=diffb) 
   	text(diffb, diff(range(fig.distF$counts))/2, adj = -0.5, expression(bold(b0)), cex=1.5 )})
   }
   plotB <- function(perm){
       plota(perm)
       plotb(perm)
   }
   switch(which, plota(perm), plotb(perm), plotB(perm))
}
