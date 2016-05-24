VUS.CutPoint <- function(x,y,z,cut.seq=NULL,optimize=F)
  {
    ######################################################################################################################################################
    ########This function calculates the optimal cut-point in VUS analysis by calculating the distance of all possible coordinates in the data to the perfect classification coordinate (1,1,1) 
    ########
    ####1. Input:
    ###x,y,z:vectors are the test marker measurements for the 3 ordinal group (D-,D0,D+)
    ###cut.seq: the sequence provided by users from which to choose the optimal cut point (usually rough and short if optimize is set to be TRUE),default=unique collection of x, y and z, should give close estimation to optimize result
    ###optimize
    ###details: the function will find among the pair of values taken from cut.seq to find the optimal cut point corresponding to min dist to (1,1,1)
    ###then a optimimiztion will be implemented to get the optimal cut point
    ####2. Output:  a named numeric of length 2 "t.minus", "t.plus". optimal cut-point from VUS analyses which has the smallest dist to the perfect classification point (1,1,1)    
    ######################################################################################################################################################


    Sp <- function(x,t.minus){sum(x<=t.minus)/length(x)}##specificity Pr(x<=t-|D-)
    Se <- function(z,t.plus){sum(z>=t.plus)/length(z)}##sensitivity Pr(z>=t+|D+)
    Sm <- function(y,t.minus,t.plus){sum(y>t.minus & y<t.plus)/length(y)}## Pr(t-<=y<=t+|D0)

    Sqr.dist <- function(beta,x0,y0,z0)
      {##beta:a vector of length 2, 1st is t.minus, 2nd is t.plus
        t.minus <- beta[1]
        t.plus <- beta[2]
        
        sp <- Sp(x0,t.minus)
        se <- Se(z0,t.plus)
        sm <- Sm(y0,t.minus,t.plus)
        (1-sp)^2+(1-se)^2+(1-sm)^2##squared distance to (1,1,1)
      }
    
    x <- na.exclude(x)
    y <- na.exclude(y)
    z <- na.exclude(z)

    n.minus <- length(x)
    n0 <- length(y)
    n.plus <- length(z)

    if(missing(cut.seq) | is.null(cut.seq)) xyz.val <- sort(unique(c(x,y,z))) else xyz.val <- cut.seq
    
    nobs <- length(xyz.val)
    
    temp.cut <- sapply(1:(nobs-1),function(i) c(rep(xyz.val[i],nobs-i),xyz.val[(i+1):nobs]))##cut1<cut2
    temp.cut <- do.call("rbind",lapply(temp.cut,function(aa) matrix(as.numeric(aa),ncol=2,byrow=F)))
    temp.cut <- as.data.frame(temp.cut)
    names(temp.cut) <- c("t.minus","t.plus")
    rownames(temp.cut) <- NULL
    dist0 <- sapply(1:nrow(temp.cut),function(i) {t.minus <- temp.cut[i,1]; t.plus <- temp.cut[i,2];Sqr.dist(beta=c(t.minus,t.plus),x,y,z)})
    min.idx <- which.min(dist0)
    #min.idx <- min.idx[length(min.idx)]##
    optimal.cut <- temp.cut[min.idx,]##rarely but can be multiple, then need to choose 
    
    if(optimize)
      {
        ##the starting point better not at the boundry 
        optimal.cut <- optim(par=c(optimal.cut[1,1],optimal.cut[1,2]),fn=Sqr.dist,x0=x,y0=y,z0=z,control=list(trace=F,reltol=1e-20,maxit=2000),method="Nelder-Mead")$par       
        optimal.cut <- data.frame(t.minus=optimal.cut[1],t.plus=optimal.cut[2])
      }

    optimal.cut <- c(optimal.cut[1,1],optimal.cut[1,2])
    names(optimal.cut) <- c("t.minus","t.plus")
    
    return(optimal.cut)    
    
  }


