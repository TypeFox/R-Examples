plot.smooth <-
function(x,...)
{
tab_smooth <- x[-1,]                           
plot(rownames(tab_smooth),rep(0,length(rownames(tab_smooth))),type="l",ylab="Smooth",xlab="TL",col="blue",ylim=range(0,max(tab_smooth)))

for (compteur in 3:ncol(tab_smooth))
lines(rownames(tab_smooth),tab_smooth[,compteur],type="l",col="blue")    
}

