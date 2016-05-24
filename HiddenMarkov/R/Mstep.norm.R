Mstep.norm <- function(x, cond, pm, pn){
    nms <- sort(names(pm))
    n <- length(x)
    m <- ncol(cond$u)
    if (all(nms==c("mean", "sd"))){
        #   estimate both
        mean <- as.numeric(matrix(x, nrow = 1) %*% 
                  cond$u)/apply(cond$u, MARGIN = 2, FUN = sum)
        sd <- sqrt(apply((matrix(x, nrow = n, ncol=m) - 
                   matrix(mean,
                   nrow = n, ncol=m, byrow=TRUE))^2 * cond$u, MARGIN=2,
                   FUN=sum)/apply(cond$u, MARGIN = 2, FUN = sum))
        return(list(mean=mean, sd=sd))
    }
    if (all(nms=="mean")){
        #   estimate mean only
        mean <- as.numeric(matrix(x, nrow = 1) %*% 
                  cond$u)/apply(cond$u, MARGIN = 2, FUN = sum)
        return(list(mean=mean))
    }
    if (all(nms=="sd")){
        #   estimate sd only
        sd <- sqrt(apply((matrix(x, nrow = n, ncol=m) - 
                   matrix(pn$mean, nrow = n, ncol=m,
                   byrow=FALSE))^2 * cond$u, MARGIN=2,
                   FUN=sum)/apply(cond$u, MARGIN = 2, FUN = sum))
        return(list(sd=sd))
    }
    stop("Invalid specification of parameters")
}

