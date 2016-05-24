plot.itemanal <- function(x, ...) {

# Plotting star plots
dev.new()
stg<-cbind(x$Scale.Stats[,2],x$Scale.Stats[,4],x$Alpha.Stats[,2:4])
rownames(stg) <- x$Alpha.Stats[,1]
colnames(stg) <- c("Mean", "Variance", "Total Corr", "Square Mult Corr", "Alpha")
stars(stg, key.labels = colnames(stg),
main = "If item deleted scale statistics", col.segments = 3:(ncol(stg)+3), draw.segments = TRUE, cex = .8)
legend(locator(1), pch=17, cex=.9, col=3:(ncol(stg)+3), bty="n",legend=colnames(stg))


# Plotting "If Item Deleted" alphas
dev.new()
idal <- matrix(stg[,5],,1)
rownames(idal) <- rownames(stg)
ylma <- round(max(idal),digits=2)
ylmi <- round(min(idal),digits=2)
plot(round(idal, digits=2), type="b", pch=19, lty=1, lwd=2, col=3,axes = FALSE, frame.plot = TRUE,
ylim=c(ylmi, ylma), ylab="Cronbach's Alpha", xlab="Item", main="Alpha if item deleted" )
axis(1, at=c(1:nrow(idal)), labels=c(rownames(idal)) )
axis(2, at=round(idal,digits=2) )


# Plotting alpha bootsrap

plot.multi.dens <- function(s) { 
junk.x = NULL
junk.y = NULL
for(i in 1:length(s))
{ junk.x = c(junk.x, density(s[[i]])$x)
junk.y = c(junk.y, density(s[[i]])$y) }
xr <- range(junk.x)
yr <- range(junk.y)
plot(density(s[[1]]), xlim = xr, ylim = yr, main = "Alpha and standardized alpha bootstraps", xlab="Alphas")
for(i in 1:length(s))
{ lines(density(s[[i]]), xlim = xr, ylim = yr, col=2+i, lwd=2  ) }}

dev.new()
paramal <- c(x$Alpha[,2],x$Alpha.Bootstrap[,3],x$Alpha.Bootstrap[,4])
paramstdal <- c(x$Std.Alpha[,2],x$Std.Alpha.Bootstrap[,3],x$Std.Alpha.Bootstrap[,4])
comp <- list(x$Bootstrap.Simmulations, x$Bootstrap.Std.Simmulations)
plot.multi.dens(comp)
abline(v=paramal, col=3, 
lty=c(1,3,3), lwd=c(1,1,1))
abline(v=paramstdal, col=4, 
lty=c(1,3,3), lwd=c(1,1,1))
legend(locator(1), c("Alpha", "Lower", "Upper", 
"Std. Alpha", "Std. Lower", "Std. Upper"), col=c(3,3,3,4,4,4), 
lty=c(1,3,3,1,3,3), lwd=c(2,2,2,2,2,2), box.lty=0)


}