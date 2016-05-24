panel.cor.res <- function(x, y, digits=2, meth="pearson", cex.cor=1)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- round(cor(x, y,method=meth),digits)
         cr<-round(cor.test(x,y,method=meth,alternative="t")$p.value,digits)
         
         text(0.5, 0.65,paste("r =",r), cex = cex.cor)
         
         if(!cr==0){
         text(0.5, 0.35,paste("p =",cr), cex = cex.cor)}
         else text(0.5, 0.35,paste("p < 0.01"), cex = cex.cor) 
         
     }

#############################################################################

panel.lm<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1,col.line = 2,lty = par("lty")) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
    beta.hat <- as.numeric(lm(y~x)$"coefficients")
    abline(beta.hat[1],beta.hat[2],col=col.line,lty=lty)	
    }

