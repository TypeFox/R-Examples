AT.gamma.response <- function( d.Gy,
			gamma.model,
			gamma.parameter,
			lethal.event.mode){

	number.of.doses	<- length(d.Gy)
	gamma.parameter <- c(gamma.parameter, rep(0, 9 - length(gamma.parameter)))

	response <- numeric(number.of.doses)

	res <- .C("AT_gamma_response_R", 
			number.of.doses = as.integer(number.of.doses),
			d.Gy = as.single(d.Gy),
			gamma.model = as.integer(gamma.model),
			gamma.parameter = as.single(gamma.parameter),
			lethal.event.mode = as.integer(lethal.event.mode),
			response = as.single(response),PACKAGE="libamtrack")

	 return.list <- NULL
	 return.list[[1]] <- res$response
	 names(return.list) <- c("response")
	 return(return.list)
}

AT.particle.name.from.particle.no <- function(particle.no){

    .Call( "AT_particle_name_from_particle_no_R",
		  particle.no,
                  PACKAGE        = "libamtrack")
}
     

AT.particle.no.from.particle.name <- function(particle.name){

     n                    <- length(particle.name)
     particle.no          <- numeric(n)
     
     for (i in 1:n){
          cur.particle.no     <- numeric(1)
          res                 <- .C( "AT_particle_no_from_particle_name_R",
		                             particle.name  = as.character(particle.name[i]),  
                                     particle.no    = as.integer(cur.particle.no),
									 PACKAGE        = "libamtrack")
          particle.no[i]     <-     res$particle.no
     }          
     return(particle.no)
}

AT.material.name.from.material.no <- function(material.no){

    .Call( "AT_material_name_from_material_no_R",
		  material.no,
                  PACKAGE        = "libamtrack")

}
     

AT.material.no.from.material.name <- function(material.name){

     n                    <- length(material.name)
     material.no          <- numeric(n)
     
     for (i in 1:n){
          cur.material.no     <- numeric(1)
          res                 <- .C( "AT_material_no_from_material_name_R",
		                             material.name  = as.character(material.name[i]),  
                                     material.no    = as.integer(cur.material.no),
									 PACKAGE        = "libamtrack")
          material.no[i]     <-     res$material.no
     }          
     return(material.no)
}
