integration <- function(values, intervals){
    n <- length(values)
    sum(0.5 * (values[1:n-1] + values[2:n]) * diff(intervals))
}

getNC <- function(beta, subbetas, nvertex, ncolor,
                  edges, neighbors=NULL, blocks=NULL, 
                  algorithm=c("SwendsenWang", "Gibbs", "Wolff"), n, burn){

      if(length(subbetas) < 2)
          stop("There has to be at least 2 intermidiate betas.")
      if(max(subbetas) > beta)
          stop("The maxmimum of the intermidiate betas has to be less than or equal to beta")
      if(max(subbetas) < beta)
          subbetas <- c(subbetas, beta)
      if(is.null(edges))
          stop("'edges' are needed to get the normalizing constant.")
      algorithm <- match.arg(algorithm)
      algorithm <- switch(algorithm, Gibbs = 1, SwendsenWang = 2, Wolff = 3)
      if(algorithm == 1 && (is.null(neighbors) || is.null(blocks)))
          stop("'neighbors' and 'blcoks' are needed to run Gibbs sampling.")
      if(algorithm == 3 && is.null(neighbors))
          stop("'neighbors' are needed to run the Wolff algorithm.")

      if(algorithm ==1)
          p.body <- quote(BlocksGibbs(n, nvertex, ncolor, neighbors, blocks, beta=subbetas[i]))
      else if(algorithm == 2)
          p.body <- quote(SW(n, nvertex, ncolor, edges, beta=subbetas[i]))
      else
          p.body <- quote(Wolff(n, nvertex, ncolor, neighbors, beta=subbetas[i]))

      EUs <- rep(0, length(subbetas))
      for (i in 1 : length(subbetas)){
           colors <- eval(p.body)
           EUs[i] <- mean(apply(colors[,-(1:burn)], 2, function(x) sum(x[edges[,1]]==x[edges[,2]])))
       }

     exp(integration(EUs, subbetas) + nvertex * log(ncolor))
}

