"plot.dmn" <- 
function(x, y, which=3, two=2, ...) {
	perm <- x
   par(mar=c(5,5,4,2)) 
   if(which==3) {
       switch(two, nada=0 , disp=par(mfrow=c(1,2)), disp=par(mfrow=c(2,1)))
   }
   plotM <- function(perm) {
   with(perm, {
   xla <- mean(range(c(diff, bootsM)))-(diff(range(c(diff, bootsM)))+diff(range(c(diff, bootsM)))*0.2)/2
   xlb <- mean(range(c(diff, bootsM)))+(diff(range(c(diff, bootsM)))+diff(range(c(diff, bootsM)))*0.2)/2
   fig.dist <- hist(bootsM, xlim=(c(xla, xlb)), 
                  main="Are means different?")
   abline(v=diff) 
   text(diff, diff(range(fig.dist$counts))/2, adj = -0.5, expression(bold(r)), cex=1.5 )})
   }
   plotF <- function(perm) {
   with(perm, {    
   xlaF <- mean(range(c(Fval, bootsF)))-(diff(range(c(Fval, bootsF)))+diff(range(c(Fval, bootsF)))*0.2)/2
   xlbF <- mean(range(c(Fval, bootsF)))+(diff(range(c(Fval, bootsF)))+diff(range(c(Fval, bootsF)))*0.2)/2
   fig.distF <- hist(bootsF, xlim=(c(xlaF, xlbF)), 
                  main="Is Fval significant?")
   abline(v=Fval) 
   text(Fval, diff(range(fig.dist$counts))/2, adj = -0.5, expression(bold(r)), cex=1.5 )})
   }
   plotB <- function(perm){
       plotM(perm)
       plotF(perm)
   }
   switch(which, plotM(perm), plotF(perm), plotB(perm))
}