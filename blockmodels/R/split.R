
coordinates_split <- function(coordinates,mass){
    o <- order(-mass)
    io <- 1:length(o)
    io[o] <- 1:length(o)

	C<-.Call( "split",
             list(coordinates=as.matrix(coordinates[o,]),
                  mass=as.matrix(mass[o])),
             PACKAGE = "blockmodels" )

    C[io,]
}

