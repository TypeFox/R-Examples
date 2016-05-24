WGS84GEO_TO_WGS84UTM <-
function(latitude, longitude){
    #### The latitude and longitude should both in decimal format. 
    phi = latitude/(180/pi)
    lambda = longitude/(180/pi)
    #########################
    #### Constants ##########
    ### WGS84LL
    N0 = 0
    E0 = 500000
    phi0 = 0
    
    if(longitude < 114){
         lambda0 = 111/(180/pi) ### (zone 49Q)
         zone = 49
    } else {
    if(longitude >= 114)
         lambda0 = 117/(180/pi) ### (zone 50Q)
         zone = 50
    }
    
    m0 = 0.9996
    M0 = 0
    niu_s = 6381309.467
    rou_s = 6344897.718
    psi_s = 1.005738745
    a = 6378137
    e2 = 6.694379989e-3
    
    A0 = 1 - ((e2)/4) - (3*(e2^2)/64)
    A2 = (3/8)*(e2 + (e2^2)/4)
    A4 = (15/256)*(e2^2)
    M = a*(A0 * phi - A2 * sin(2 * phi) + A4 * sin(4*phi))
    M0 = a*(A0 * phi0 - A2 * sin(2 * phi0) + A4 * sin(4*phi0))
    
    #### Eq. 1
    N = N0 + m0*((M - M0) + niu_s*(sin(phi))*((lambda - lambda0)^2/2)*(cos(phi)))
    #### Eq. 2
    E = E0 + m0*(niu_s*(lambda - lambda0)*cos(phi) + niu_s*((lambda - lambda0)^3/6)*(cos(phi)^3)*(psi_s - tan(phi)^2))
    res <- list(N, E, zone)
    names(res) <- c("N", "E", "zone")
    return(res)
}
