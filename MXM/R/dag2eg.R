dag2eg <- function(dag, type = NULL) {

  ## dag is a square matrix corresponding to a DAG
  if ( nrow(dag) != ncol(dag) ) {
    essential <- paste("This matrix is not square")
    
  } else {
    if ( sum( dag >= 2 ) > 0 ) {
      typos = 1  
      g3 <- which( dag == 3 )
      dag[g3] <- 0
      g2 <- which( dag == 2 )
      dag[g2] <- 1
    } else  typos = 2
    
      essential <- ggm::essentialGraph(dag)

      if ( is.null(type) ) {
        
        if ( typos == 1 ) { 
          eg <- essential + t(essential)
          a <- which(eg == 2) 
          b <- which(eg == 1, arr.ind = TRUE)   
          b <- t( apply(b, 1, sort ) )         
          b <- unique(b )
          eg[cbind(b[, 2], b[, 1]) ] <- 3
          eg[cbind(b[, 1], b[, 2]) ] <- 2
          eg[ a ] <- 1
          essen <- eg
           
        } else if (typos == 2) {
          essen <- essential 
        }

      } else if ( typos == 1 ) {
        essen <- eg 
      } else {
        essen <- essential 
      }
  }    
 
  essen

} 
       
     
  
  