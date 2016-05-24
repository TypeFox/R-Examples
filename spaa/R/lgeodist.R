lgeodist <-
function(L1, phi1, L2, phi2)
{
    Ri = 6371
    d = ( Ri * acos(sin(phi1*(pi/180))*sin(phi2*(pi/180)) + 
         cos(phi1*(pi/180))*cos(phi2*(pi/180))*cos((L1 - L2)*(pi/180))))
    res <- round(d, 3)
    return(res)
}

