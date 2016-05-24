estimate.shiftvec <- function(predator, prey, num.samples=1000, num.preysamples=1, pred.distinput=FALSE, prey.distinput=TRUE) {
    ## given two data frames, one containing predator samples and one containing prey samples,
    ## randomly create shift vectors through the monte-carlo process
    ## d is the dimension (number of isotopes)
    ## num.samples is how many shift samples we want to return
    ## num.preysamples is the number of samples to draw from each of the prey types during each iteration
    ## if pred.distinput is true, then predator is a 1 x (2*d + 1) data frame
    ##    the first column contains the name of the predator source
    ##    the next 2*d columns contain the mean and standard distribution of each isotope, where column
    ##    2i contains the mean and column 2i+1 contains the sd of the ith isotope
    ## if prey.distinput is true, then let s be the number of prey types. prey is a s x (2*d+1) data frame
    ##    the first column contains the names of each of the prey sources
    ##    the next 2*d columns contain the mean and standard distribution of each isotope, organized in the same
    ##    way as the predator data frame
    ## randomly sample isotopes from predator and prey
    
    cur.samples = 0
    allshifts = c()
    ## the dimension of the data (= the number of isotopes tested)
    d = dim(predator)[2]
    ## assume predator is a m x d data frame
    ## assume prey is a n x d*s data frame, where
    ## s is the number of prey types
    while(cur.samples < num.samples) {
        ## we first partition the prey frame into it's
        ## prey subframes
        predator = na.omit(predator)
        ## sample a point from the predator gaussian distribution
        if (pred.distinput) {
            ## parse the input matrix and get the predator distribution
            parsed = parse.isotopeframe(predator)
            predmu = parsed$mu
            predsigma = array(0, dim=c(d,d))
            diag(predsigma) = parsed$sigma
        } else {
            ## compute the predator distribution from the data
            predmu = colMeans(predator)
            predsigma = var(predator)
            m = dim(predator)[1]
        }
        rpred = mvrnorm(1, mu=predmu, Sigma=predsigma**2)
        ## R sucks and is stupid. Look at what I have to do.
        rpred = data.frame(t(as.matrix(rpred)))
        if (prey.distinput) {
            s = dim(prey)[1]
        } else {
            s = dim(prey)[2]/d
        }
        rprey = array(0, dim=c(s*num.preysamples,d))
        ## sample a number of points from each prey type by sampling from a gaussian on the type
        for (j in 1:s) {
            ## calculate the mean and sd
            if (prey.distinput) {
                ## parse the input matrix and get the prey distribution
                parsed = parse.isotopeframe(prey)
                preymu = parsed$mu[j,]
                preysigma = array(0, dim=c(d,d))
                diag(preysigma) = parsed$sigma[j,]
            } else {
                ## compute the prey distribution from the data
                tmp = na.omit(prey[,c(2*j-1, 2*j)])
                preymu = colMeans(tmp)
                preysigma = var(tmp)
            }
            rprey[((j-1)*num.preysamples+1):(j*num.preysamples),] = mvrnorm(num.preysamples, mu=preymu, Sigma=preysigma**2)
        }
        out = shift.mean(rpred, data.frame(rprey))
        if (!all(is.na(out))) {
            allshifts = rbind(allshifts, out)
            cur.samples = cur.samples + 1
        }
    }
    allshifts = allshifts[-1,]
    rownames(allshifts) = seq(dim(allshifts)[1])
}
