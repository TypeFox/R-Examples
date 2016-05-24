require(lpSolve)

lp2 <- function(direction = "min", objective.in, const.mat, const.dir, 
                const.rhs, free.var = NULL){

  if(is.null(free.var))
    lp2.out <- lp(direction = direction, objective.in = objective.in,
                const.mat = const.mat, const.dir = const.dir,
                const.rhs = const.rhs)
  
  if(!is.null(free.var)){
    if(length(free.var) != length(objective.in))
      stop("The length of free.var should be equivalent to that of objective.in")

    free.var <- as.logical(free.var)
    nn <- length(objective.in)
    
    objective.in2 <- vector()
    const.mat2 <- matrix(, nrow = nrow(const.mat))
    my.idx <- vector()

    for(i in 1:nn){
      if(free.var[i] != 1){
        objective.in2 <- c(objective.in2, objective.in[i])
        const.mat2 <- cbind(const.mat2, const.mat[, i])
        if(i == 1)
          const.mat2 <- as.matrix(const.mat2[, -1], nrow = nrow(const.mat))
        my.idx <- c(my.idx, 0)
      }
      else{
        objective.in2 <- c(objective.in2, objective.in[i],
                           -objective.in[i])
        const.mat2 <-
          cbind(const.mat2, const.mat[, i], -const.mat[, i])
        if(i == 1)
          const.mat2 <- as.matrix(const.mat2[, -1], nrow = nrow(const.mat))
        my.idx <- c(my.idx, 1, -1)
      }
    }

    lp2.out <- lp(direction = direction, objective.in = objective.in2,
                  const.mat = const.mat2, const.dir = const.dir,
                  const.rhs = const.rhs)

    re <- list()
    objval <- sum(lp2.out$objval)
    tmp.solution <- lp2.out$solution
    
    solution <- rep(NA, nn)
    sol.0 <- tmp.solution[my.idx == 0]
    sol.p1 <- tmp.solution[my.idx == 1]
    sol.n1 <- tmp.solution[my.idx == -1]
    sol.1 <- sol.p1 - sol.n1
    solution[free.var == 0] <- sol.0
    solution[free.var == 1] <- sol.1
    lp2.out$solution.free <- solution
  }
  return(lp2.out)
}


