# 2-3-2008 MRC-Epid JHZ

KCC <- function(model,GRR,p1,K)
{
   model.idx <- charmatch(model,c("multiplicative","additive","recessive","dominant","overdominant"))
   if(is.na(model.idx)) stop("Invalid model type")
   if(model.idx == 0) stop("Ambiguous model type")
   multiplicative <- c(1,GRR,GRR*GRR)
   additive <- c(1,GRR,2*GRR-1)
   recessive <- c(1,1,GRR)
   dominant <- c(1,GRR,GRR)
   overdominant <- c(GRR,1,GRR)
   f <- switch(model.idx,multiplicative,additive,recessive,dominant,overdominant)
   scale <- K/(f[1]*(1-p1)^2+f[2]*2*p1*(1-p1)+f[3]*p1^2)
   f <- f*scale
   if(f[3]>1) stop("misspecified model")
   pprime <- (f[3]*p1^2+f[2]*p1*(1-p1))/K
   p <- ((1-f[3])*p1^2+(1-f[2])*p1*(1-p1))/(1-K)
   invisible(list(pprime=pprime,p=p))
}
