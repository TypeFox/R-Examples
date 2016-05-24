#' @title Components of beta-diversity for a list of networks
#' @description
#' Given a list of networks, returns the pairwise beta-diversity components
#'
#' @param N a list of networks
#' @param complete (boolean) whether all combinations of networks should be tested
#' @param ... additional arguments to be passed to \link{betalink}
#'
#' @return A dataframe with the pairwise distances
#'
#' @export
network_betadiversity <- function(N, complete=FALSE, ...){
   N <- name_networks(N)
   beta <- NULL
   stop_at <- ifelse(complete, length(N), length(N)-1)
   for(i in c(1:stop_at))
   {
      start_at <- ifelse(complete, 1, i+1)
      inner_stop <- length(N)
      inner_steps <- c(start_at:inner_stop)
      inner_steps <- inner_steps[which(inner_steps!=i)]
      for(j in inner_steps)
      {
         b <- betalink(N[[i]], N[[j]], ...)
         b$i <- names(N)[i]
         b$j <- names(N)[j]
         beta <- rbind(beta, rbind(unlist(b)))
      }
   }
   if(NROW(beta) == 1) {
      beta <- data.frame(t(beta[,c('i', 'j', 'S', 'OS', 'WN', 'ST')]))
   } else {
      beta <- data.frame(beta[,c('i', 'j', 'S', 'OS', 'WN', 'ST')])
   }
   beta$OS <- as.numeric(as.vector(beta$OS))
   beta$S <- as.numeric(as.vector(beta$S))
   beta$WN <- as.numeric(as.vector(beta$WN))
   beta$ST <- as.numeric(as.vector(beta$ST))
   return(beta)
}
