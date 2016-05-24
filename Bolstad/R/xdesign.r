#' Monte Carlo study of randomized and blocked designs
#' 
#' Simulates completely randomized design and randomized block designs from a
#' population of experimental units with underlying response values \eqn{y} and
#' underlying other variable values \eqn{x} (possibly lurking)
#' 
#' 
#' @param x a set of lurking values which are correlated with the response
#' @param y a set of response values
#' @param corr the correlation between the response and lurking variable
#' @param size the size of the treatment groups
#' @param n.treatments the number of treatments
#' @param n.rep the number of Monte Carlo replicates
#' @return If the ouput of xdesign is assigned to a variable, then a list is
#' returned with the following components: \item{block.means}{a vector of the
#' means of the lurking variable from each replicate of the simulation stored
#' by treatment number within replicate number} \item{treat.means}{a vector of
#' the means of the response variable from each replicate of the simulation
#' stored by treatment number within replicate number} \item{ind}{a vector
#' containing the treatment group numbers. Note that there will be twice as
#' many group numbers as there are treatments corresponding to the simulations
#' done using a completely randomized design and the simulations done using a
#' randomized block design}
#' @keywords misc
#' @examples
#' 
#' # Carry out simulations using the default parameters 
#' 
#' xdesign()
#' 
#' # Carry out simulations using a simulated response with 5 treaments, 
#' # groups of size 25, and a correlation of -0.6 between the response 
#' # and lurking variable
#' 
#' xdesign(corr = -0.6, size = 25, n.treatments = 5)
#' 
#' @export xdesign
xdesign = function(x = NULL, y = NULL, corr = 0.8, size = 20,
                   n.treatments = 4, n.rep = 500){

    if(is.null(x)){ ## simulate some data
        nx = size*n.treatments
        x = rnorm(nx)
        y = rnorm(nx)

        y = sqrt(1 - corr^2) * y + corr * x
    }

    nx = size * n.treatments
    if(length(x) != length(y))
        stop("x and y must be of equal length")
    
    if(length(x) != size * n.treatments)
        stop("x and y must be equal to the same size times the number of treatments")
    
    if(corr < (-1) | corr > 1)
        stop("Correlation coeficient must be between -1 and 1")

    if(n.rep < 10)
        stop("Must have at least 10 Monte Carlo replicates")

    cat("Variable\tN\tMean\tMedian\tTrMean\tStDev\tSE Mean\n")
    cat(paste("X\t",length(x),
              round(mean(x),3),
              round(median(x),3),
              round(mean(x,trim=0.1),3),
              round(sd(x),3),
              round(sd(x)/sqrt(length(x)),3),sep="\t"))
    cat("\n")
    cat(paste("Y\t",length(y),
              round(mean(y),3),
              round(median(y),3),
              round(mean(y,trim=0.1),3),
              round(sd(y),3),
              round(sd(y)/sqrt(length(y)),3),sep="\t"))
    cat("\n\n")

    qx = quantile(x,c(0.25,0.75))
    qy = quantile(y,c(0.25,0.75))

    cat("Variable\tMinimum\tMaximum\tQ1\tQ3\n")
    cat(paste("X\t",round(min(x),3)
              ,round(max(x),3)
              ,round(qx[1],3)
              ,round(qx[2],3),sep="\t"))
    cat("\n")
    cat(paste("Y\t",round(min(y),3)
              ,round(max(y),3)
              ,round(qy[1],3)
              ,round(qy[2],3),sep="\t"))

    cat("\n\n")

    cat("The Pearson correlation between X and Y is: ");
    cat(paste(round(cor(x,y),3),"\n\n"))

    plot(x, y)

    ssx = rep(0,n.rep)
    ssy = rep(0,n.rep)

    treat.groupmean = matrix(0,ncol=n.treatments,nrow=n.rep)
    block.groupmean = matrix(0,ncol=n.treatments,nrow=n.rep)

    for(block in c(FALSE,TRUE)){
        ## block is indicator for blocking
        ## FALSE =  completely randomized design,
        ## TRUE = randomized block design

        for(i in 1:n.rep){
            if(!block){
                group = rep(1:n.treatments,size)

                z = rnorm(nx)
                o = order(z)
                z = z[o]
                group = group[o]

                x2 = x
                y2 = y
            }else{
                o = order(x)
                x2 = x[o]
                y2 = y[o]

                group = NULL

                for(j in 1:size){
                    gp = 1:n.treatments
                    z = rnorm(n.treatments)
                    gp = gp[order(z)]
                    group = c(group,gp)
                }
            }

            split.x = split(x2,group)
            split.y = split(y2,group)

            x.bar = sapply(split.x,mean)
            y.bar = sapply(split.y,mean)

            x.mean = mean(x.bar)
            y.mean = mean(y.bar)

            ssx[i] = sum((x.bar-x.mean)^2)
            ssy[i] = sum((y.bar-y.mean)^2)

            treat.groupmean[i,] = y.bar
            block.groupmean[i,] = x.bar
        }

        if(!block){
            treat.var0 = as.vector(treat.groupmean)
            block.var0 = as.vector(block.groupmean)
            index0 = rep(1:n.treatments,rep(n.rep,n.treatments))
        }else{
            treat.var1 = as.vector(treat.groupmean)
            block.var1 = as.vector(block.groupmean)
            index1 = rep(1:n.treatments,rep(n.rep,n.treatments))
        }
    }

    treat.var = c(treat.var0,treat.var1)
    block.var = c(block.var0,block.var1)
    index = c(index0,index1)

    ind = rep(1:2,c(length(treat.var0),length(treat.var1)))
    ind = n.treatments*(ind-1)+index

    par(ask=interactive())

    rng = range(block.var)
    y.lims = max(abs(c(rng[1]-0.1*diff(rng),rng[2]+0.1*diff(rng))))
    y.lims = c(-y.lims,y.lims)

    boxplot(block.var~ind,
            main="Boxplots of Lurking/Blocking variable group means",
            sub="Lurking variable in completely randomized design\nBlocking variable in randomized block design",
            col=rep(c("blue","red"), rep(n.treatments,2)),
            ylim = y.lims)
    legend("topright", bty = "n", cex = 0.7, 
           legend = c("Completely randomized design",
                      "Randomized block design"),
           fill = c("blue","red"))

    rng = range(treat.var)
    y.lims = max(abs(c(rng[1]-0.1*diff(rng),rng[2]+0.1*diff(rng))))
    y.lims = c(-y.lims,y.lims)
    boxplot(treat.var~ind
            ,main="Boxplots of treatment group means"
            ,col=rep(c("blue","red"),rep(n.treatments,2))
            ,ylim=y.lims)
    legend("topright", cex = 0.7, bty = "n",
           legend=c("Completely randomized design",
                    "Randomized block design"),
           fill = c("blue","red"))

    x = treat.var[ind<=n.treatments]
    y = treat.var[ind>n.treatments]
    cat("Variable\tN\tMean\tMedian\tTrMean\tStDev\tSE Mean\n")
    cat(paste("Randomized",length(x),
              round(mean(x),3),
              round(median(x),3),
              round(mean(x,trim=0.1),3),
              round(sd(x),3),
              round(sd(x)/sqrt(length(x)),3),sep="\t"))
    cat("\n")
    cat(paste("Blocked\t",length(y),
              round(mean(y),3),
              round(median(y),3),
              round(mean(y,trim=0.1),3),
              round(sd(y),3),
              round(sd(y)/sqrt(length(y)),3),sep="\t"))
    cat("\n\n")

    qx = quantile(x,c(0.25,0.75))
    qy = quantile(y,c(0.25,0.75))

    cat("Variable\tMinimum\tMaximum\tQ1\tQ3\n")
    cat(paste("Randomized",round(min(x),3)
              ,round(max(x),3)
              ,round(qx[1],3)
              ,round(qx[2],3),sep="\t"))
    cat("\n")
    cat(paste("Blocked\t",round(min(y),3)
              ,round(max(y),3)
              ,round(qy[1],3)
              ,round(qy[2],3),sep="\t"))

    cat("\n\n")

    invisible(list(block.means=block.var,treat.means=treat.var,ind=ind))
}
