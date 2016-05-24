TestPSOParameters <-
  function(Population,Iteration){
    # Check the PSO parameters 'Population' and 'Iteration'.
    # 
    # input
    #   Population: numeric. The number of particles.
    #   Iteration: numeric. The number of iterations.
    #
    #    For example, Population <- 1000,  Iteration <- 50
    #
    # Junliang Shang
    # 10.30/2014
    
    CharNum <- nchar(Population)
    DotNum <- 0
    for (j in 1:CharNum){
      status <- (substr(Population,j,j) %in% 
                   c("0","1","2","3","4","5","6","7","8","9","."))
      if (status==FALSE){
        stop("   Error! ",Population," should be a number.\n")
      }
      if (substr(Population,j,j)=="."){
        DotNum <- DotNum+1
      }
    }
    if (DotNum>1){
      stop("   Error! ",Population," should be a number.\n")
    }
    if (floor(as.numeric(Population))
        !=as.numeric(Population)){
      stop("   Error! ",Population," should be an integer.\n")
    } 
    
    
    CharNum <- nchar(Iteration)
    DotNum <- 0
    for (j in 1:CharNum){
      status <- (substr(Iteration,j,j) %in% 
                   c("0","1","2","3","4","5","6","7","8","9","."))
      if (status==FALSE){
        stop("   Error! ",Iteration," should be a number.\n")
      }
      if (substr(Iteration,j,j)=="."){
        DotNum <- DotNum+1
      }
    }
    if (DotNum>1){
      stop("   Error! ",Iteration," should be a number.\n")
    }
    if (floor(as.numeric(Iteration))
        !=as.numeric(Iteration)){
      stop("   Error! ",Iteration," should be an integer.\n")
    } 
    
  }