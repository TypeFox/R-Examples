# Verified 1.3.18
get.shape.venezuela <-
        function(shape.file="venezuela") {
                return(maptools::readShapeSpatial(paste(system.file("shape",package="vetools"), shape.file, sep="/")))
        }