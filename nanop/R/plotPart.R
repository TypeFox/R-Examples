plotPart <- function(nanop, radius=0.4, legend=TRUE, col=NA, box=FALSE, play=FALSE, atoms=NA, miller=NA, lattice=c(4.08)){
  sym <- attributes(nanop)$sym
  Ncore <- length(which(attributes(nanop)$atomType>0)) #number of atoms in the core
  N <- attributes(nanop)$nAtomTypes #number of atoms types
  if(length(radius) < N)
    radius <- rep(radius, N)
  if(is.na(col[1]) || length(col)!=N) 
    col=1:N 
  types <- 1:max(attributes(nanop)$atomType)
  if(min(attributes(nanop)$atomType) < 0)
    types <- append(types, -1:min(attributes(nanop)$atomType))
  
  colInd <- 0
  for(i in 1:N)   
    colInd[which(attributes(nanop)$atomType==types[i])] <- col[i]   
  
  rr <- radius
  for(i in 1:N)   
    radius[which(attributes(nanop)$atomType==types[i])] <- rr[i] 
   
  incl <- types
  extractPlane <- 1:nrow(nanop)
  extract <- 1:nrow(nanop)
  type <- 0
  namesCore <-attributes(nanop)$atomsCore
  if(!is.null(attributes(nanop)$atomsShell)){ 
    namesShell <- attributes(nanop)$atomsShell
  }else 
    namesShell <- NA
	
  if(!is.na(atoms[1])){
    for(i in 1:length(namesCore)){
	  if(atoms == paste("core", namesCore[[i]]) )
        type <- i	  
	}
	if(!is.na(namesShell[1])){
      for(i in 1:length(namesShell)){
	    if(atoms == paste("shell", namesShell[[i]]) )
          type <- -i	  
	  }
	}
    if(type==0)
	  stop("no atoms of given type in core/shell \n", immediate. = TRUE)
	
	extract <- which(attributes(nanop)$atomType==type) 
	incl <- type
  }
  if(!is.na(miller[1])){
    extractPlane <- plotPlane(nanop[1:Ncore, ], miller=miller, sym=sym, lattice=lattice, type=type)
	if(is.na(atoms[1])) incl <- types[which(types>0)]
  }
  extract <- intersect(extract, extractPlane)
  
  if(length(extract)<1) 
    stop("no atoms under specified conditions \n", immediate. = TRUE)
  
  r3dDefaults <- list()
  r3dDefaults <<- rgl::r3dDefaults
  rgl::open3d()

  if(!box)
    res <- rgl::spheres3d(nanop[extract,], radius=radius[extract], color=colInd[extract])
  else
    res <- rgl::plot3d(nanop[extract,], type="s", radius=radius[extract], xlab="x", ylab="y", zlab="z", col=colInd[extract])
	 
#   select=FALSE
#   if(select){
#     selectpoints3d(res["data"],
#	        multiple = function(x) {
#	                spheres3d(x, color = "red", alpha = 0.3,  radius = 0.6)
#	                TRUE}
#					)
#					
# }
  
  if(legend){
    x <- y <- z <- max(attributes(nanop)$r)
	if(attributes(nanop)$shape=="box" || attributes(nanop)$shape=="ellipse"){
	  x <- attributes(nanop)$r[1]
	  y <- attributes(nanop)$r[2]
	  z <- attributes(nanop)$r[3]	
	}
	
    names <- namesCore
    if(!is.na(namesShell[1]))
	  names <- c(namesCore, namesShell)
	  
    x = 1.1*x
    y =-1.1*y
    z = 1.1*z
    for(i in incl){
      rgl::plot3d(x=x, y=y, z=z, type="s", radius=min(radius), add=TRUE,  col=col[which(types==i)])
      rgl::text3d(x=x+0.4*x, y=y, z=z, font=4, texts=names[which(types==i)])
	  z=z-0.2*z
    }
  }
  
  if(play){
     M <- rgl::par3d("userMatrix")
     rgl::play3d(rgl::spin3d(axis = c(1, 1, 1), rpm = 2))   
   }
   
   
 }
 

##############################################################################################
plotPlane <- function(nanop, miller=c(1,1,1), sym, lattice, type=type){
  
  if(sym == "fcc" || sym == "bcc" || sym == "sc"){
    struct <- "cubic"
	#matrix describing new basis in terms of cartesian basis
	A <- matrix(c(c(1,0,0), 
	              c(0,1,0),
				  c(0,0,1)), byrow=TRUE, ncol=3)
	Ap <- solve(t(A)) # inverse and transpose
	                  # Ap -- convertion to new basis matrix
  }	
  if(sym == "hcp"){
    struct <- "hexagonal"
	A <- matrix(c(c(sqrt(3)/2,-0.5,0), 
	              c(0,1,0),
				  c(-sqrt(3)/2,-0.5,0),
				  c(0,0,1)), byrow=TRUE, ncol=3)
	Ap <-  matrix(c(c(1/sqrt(3),0,0), 
				  c(0,1,0),
	              c(-1/sqrt(3),1,0),
				  c(0,0,1)), byrow=TRUE, ncol=3)
				  
    Ap <-  matrix(c(c(sqrt(3)/2,-0.5,0), 
				  c(0,1,0),
	              c(-0.5,-sqrt(3)/2,0),
				  c(0,0,1)), byrow=TRUE, ncol=3)
	## coords in new basis = c(x,y,u,z) = Ap * coords in cartesian basis
  }
  
  nanopNB <- t(Ap %*% t(nanop))
  	
  extractPlane <- 0
  for(i in 1:nrow(nanopNB)){
    if(struct != "hexagonal"){
      dist <- nanopNB[i,1]*miller[1]/lattice[1] + nanopNB[i,2]*miller[2]/lattice[1] + nanopNB[i,3]*miller[3]/lattice[1]
	  if(type==2 && sym !="CaTiO3") 
	    dist <- dist - 0.5
#	  if(type==2 && sym =="CaTiO3") 
#	    dist <- dist - sqrt(3)/2
#	  if(type==3 && sym =="CaTiO3") 
#	    dist <- dist - sqrt(2)/2	
		
	}
    else{
      dist <- nanopNB[i,1]*miller[1]/lattice[1] + nanopNB[i,2]*miller[2]/lattice[1] + nanopNB[i,3]*miller[3]/lattice[1] + nanopNB[i,4]*miller[4]/lattice[2]
	  if(type==2) 
	    dist <- dist - lattice[3]
	}
	
	if(abs(dist) < 0.1) extractPlane <- append(extractPlane, i) 
  
  }
  extractPlane <- extractPlane[-1]
  extractPlane  
  
 }
 