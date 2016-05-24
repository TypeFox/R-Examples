##' create particles for povray input
##'
##' writes a list of particles
##' @title particles_povray
##' @param positions matrix of positions
##' @param angles matrix of Euler angles in radians
##' @param sizes matrix of particle sizes
##' @param out output filename
##' @return side-effect only
##' @author baptiste Auguie
##' @export
##' @family user_level povray
particles_povray <- function(positions,
                      angles,
                      sizes, out="positions.pov"){

cat(paste("object{ Particle scale <", apply(round(sizes, 5), 1, paste, collapse=", "), "> EulerRotate(<",
          apply(round(angles*180/pi, 5), 1, paste, collapse=", "), ">) translate <",
          apply(round(positions, 5), 1, paste, collapse=", "), 
          "> }", sep=""), sep="\n",
    file=out, append=FALSE)

}

##' create long string of small particles for povray input
##'
##' writes a list of particles, as well as cylinder
##' @title curve_povray
##' @param positions matrix of positions
##' @param size radius of particles
##' @param radius radius of inner cylinder
##' @param out output filename
##' @return side-effect only (note append=TRUE)
##' @author baptiste Auguie
##' @export
##' @family user_level povray
curve_povray <- function(positions, size=0.005, radius, 
                       out="positions.pov"){

  sizes <- size + 0*positions
cat(paste("object{ Chain scale <", apply(round(sizes, 5), 1, paste, collapse=", "), ">  translate <",
          apply(round(positions, 5), 1, paste, collapse=", "), 
          "> }", sep=""), sep="\n",
    file=out, append=TRUE)

cat(paste("
cylinder
    {
        <0,0,", 1.2*min(positions[, 3]), ">,
       <0,0,", 1.2*max(positions[, 3]), ">,
       ", 
          radius, "
        open
        texture{
                    pigment{color rgbf <0.95, 0.95, 0.95, 0.95>}
                    finish{
                        reflection {0,0.01} 
                        //conserve_energy 
                        ambient 0.1 diffuse 0.9   
                        specular 0.2 roughness 1/100
                    }                                        
                }     
                interior{ior 1.0}// fade_power 2 fade_distance 2}
    }

", sep=""), sep="\n",
    file=out, append=TRUE)
}
