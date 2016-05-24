projectn<-
function (v0, mat, matsd = NULL, estamb = FALSE, estdem = FALSE, 
    equalsign = TRUE, stmat = NULL, fecundity1 = TRUE, nrep = 1, 
    time = 10, management=NULL, round=TRUE) 
{
          #management=matrix with adition or deletion of oindividuals to be applied in each generation
	  #round: should the projections be rounded to the nex integer (i.e. consider finite indiviuals)?

    if (sum(class(mat) == "tmatrix") == 0) 
        mat <- as.tmatrix(mat)
	
    # initialite list 	to store replication results
    vn <- NULL
    vm <- NULL # para almacenar los resultados del manejo (harvest)
    for (i in 1:nrep) {
        vn[[i]] <- cbind(v0, v0)
	vm[[i]] <- cbind(v0*NA, v0*NA)
    }
    for (i in 1:time) {
       
	for (ii in 1:nrep){
	
		 v <- project1(v0 = vn[[ii]][, i + 1], mat = mat, matsd = matsd,
		                    estamb = estamb, estdem = estdem, 
		                    equalsign = equalsign, stmat = stmat, 
		                    fecundity1 = fecundity1)
		if(round==TRUE) v<- round(v)
		v.0 <- v # v antes del manejo
		if(!is.null  (management)) {
		  #if mangement is a vector, transform it to a column matrix
                   management <- matrix(management, nrow=dim(mat)[1])
		   # apply management actions sequentialy
		   for ( j in 1:dim(management)[2]){
		         # which stages are managed on a proportion or on a individual basis
			  proportions <- abs(management[,j]) < 1 
			  # transfom proportions to individuals
		          management[proportions,j] <- (round(v*management))[proportions]
			  v <- v+management[,j]
			  v[v<0] <-0
				
			}
                }
		
		v.m <- v.0-v # el reultado del manejo es la diferencia entre el v resultante y el v antes del manejo
		harvest <- v.m
		harvest[v.m < 0] <- 0
	
		
	       vn[[ii]] <- cbind(vn[[ii]], v)
	       vm[[ii]] <- cbind(vm[[ii]], harvest)
	
          }	
	
    
        }
        vn <- lapply(vn, function(x) x[, -1])
        vm <- lapply(vm, function(x) x[, -1])

       vnm<-list(vn=vn, harvest=vm, mat=mat, management=management)
       class(vnm) <- c("rmas", class(vnm))
      return(vnm)


}
