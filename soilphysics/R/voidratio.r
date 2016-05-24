voidratio <- 
function (wetsoil, drysoil, diam.cylinder, height.cylinder, dens.particle, 
    deformation) 
{
    area <- pi * (diam.cylinder/2)^2
    aar <- height.cylinder - deformation
    vra <- aar * area
    drs <- drysoil/vra
    ei <- dens.particle/drs - 1
    return(ei)
}