## rhombus formula
RhombusFormula <- function(cor.matrix, obs)
{
    if (obs<=1e-4) return(1)

    cor.matrix[cor.matrix>1-1e-10]  <- 1-1e-10
    cor.matrix[cor.matrix<-1+1e-10] <- -1+1e-10

    k <- ncol(cor.matrix) # number of statistics
    len.matrix <- acos(cor.matrix)

    loc.g <- function(x)
    {
        x^2 #+ 2*x^4/3 + 17*x^6/45 + 41*x^8/320 + 61*x^10/1728 + 3721*x^12/518400
    }

    loc.f <- function(L, x)
    {
        if (L<=pi/2) qet <- 2*(pnorm(x*L/2)-0.5) + exp(-(x^2)*loc.g(L/2)/2)*(pnorm(x*(pi-L)/2)-pnorm(x*L/2))
	      else  qet <- exp(-(x^2)*loc.g((pi-L)/2)/2)*(pnorm(x*L/2)-pnorm(x*(pi-L)/2)) + 2*(pnorm(x*(pi-L)/2)-0.5)
        qet/x
    }

    loc.phi <- function(dex)
    {
        L <- len.matrix[dex[1],dex[2]]
        loc.f(L, obs)
    }

    loc.app <- function(x)
    {
        b <- 0
        for (i in k:2)
      	{
	          y <- cbind(x[1:(i-1)], x[i])
	          a <- apply(y,1,loc.phi)
	          b <- b + min(a)
	      }
	      b
    }
    
    comb <- matrix(unlist(combinat::permn(k)),ncol=k, byrow=TRUE)   #????????
    res  <- apply(comb, 1, loc.app)

    pv <- (k-2)*(pnorm(obs)-pnorm(-obs)-1) + 4*dnorm(obs)*min(res)
    
    min(max(pv, 0),1)
}
