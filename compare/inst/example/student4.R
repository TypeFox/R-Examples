 id <- c(1, 2, 3, 4, 5, 6)
 
 age <- c(30, 32, 28, 39, 20, 25)
 
 edu <- c(0, 0, 0, 0, 0, 0)
 
 classNames <- c( "poor", "middle")
 class <- rep( classNames, c(3, 3)) 
 factor( class, levels= classNames)
 
 indianMothers <- data.frame(id=id,
                              age,
                              edu,
                              class) 
 
 
