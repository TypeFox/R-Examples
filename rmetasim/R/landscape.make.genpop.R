# Create genpop object from landscape

landscape.make.genpop <- function(Rland)
    {
        genind2genpop(landscape.make.genind(Rland))
    }
