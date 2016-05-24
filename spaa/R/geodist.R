geodist <-
function(L1, phi1, L2, phi2){
    a = 6378.14
    f = 1/298.257
    F = (phi1+phi2)/2
    G = (phi1 - phi2)/2
    ramda <- (L1 - L2)/2
    
    S = (sin(G*pi/180)^2)*(cos(ramda*pi/180)^2) + 
        (cos(F*pi/180)^2)*(sin(ramda*pi/180)^2)
    C = (cos(G*pi/180)^2)*(cos(ramda*pi/180)^2) + 
        (sin(F*pi/180)^2)*(sin(ramda*pi/180)^2)
    
    omega = atan(sqrt(S/C))
    R = sqrt(S*C)/omega
    D = 2*omega*a
    
    H1 = (3*R-1)/(2*C)
    H2 = (3*R+1)/(2*S)
    
    res = D*(1 + f*H1*(sin(F*pi/180)^2)*(cos(G*pi/180)^2) - 
      f*H2*(cos(F*pi/180)^2)*(sin(G*pi/180)^2))
    return(round(res,3))
}

