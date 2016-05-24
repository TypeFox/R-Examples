WGS84UTM_TO_WGS84GEO <-
function(N, E, zone = c(49, 50)){
    ### The unit for N and E is meter
    #########################
    #### Constants ##########
    ### WGS84LL
    N0 = 0
    E0 = 500000
    phi0 = 0
    if(zone == 49){
        lambda0 = 111 ## if (zone 49Q)
    } 
    if(zone == 50){
        lambda0 = 117 ## if (zone 50Q)
    } 
    m0 = 0.9996
    M0 = 0
    niu_s = 6381309.467
    rou_s = 6344897.718
    psi_s = 1.0053738745
    a = 6378137
    e2 = 6.694379989e-3

    A0 = 1 - (e2/4) - (3*(e2^2)/64)
    A2 = (3/8)*(e2 + (e2^2)/4)
    A4 = (15/256)*(e2^2)
    delta_N = N - N0
    fm <- function(x){
         return((((delta_N + M0) / m0)/a + A2 * sin (2*x) - A4 *sin (4*x))/A0 )
    }
    #### iterations
    iterate <- function(x, d){  
        a=x;  
        b=fm(a);  
        k=0; #//Count the number of loops  
        while(( (a-b) > d) | ((a-b) < -1 * d)){  
            ### print(a);  
            a=b;  
            b=fm(a);  
            k = k + 1; 
            if(k>1000){  
                stop("The equation does not converge. Stopped computation."); 
                return(0);
            }  
        }  
       return(b);    
    }  
    phi_rou <- iterate(0.5, d = 1e-30)
    ######################################
    ##### Verification that the answer is correct
    #####  a*(A0 * phi_rou - A2 * sin(2 * phi_rou) + A4 * sin(4*phi_rou))
    ##### (delta_N + M0) / m0
    
    t_rou <- tan(phi_rou)
    niu_rou <- a/sqrt(1 - e2*(sin(phi_rou)^2))
    rou_rou <- a*(1 - e2)/((1 - e2*(sin(phi_rou)^2))^(3/2))
    psi_rou <- niu_rou/rou_rou
    
    delta_E = E - E0
    
    #### Equation 4. 
    lambda = (lambda0/(180/pi) + (1/cos(phi_rou))*(delta_E/(m0*niu_rou)) - 
             (1/cos(phi_rou))*(delta_E^3/(6*m0^3*niu_rou^3))*(psi_rou + 2*t_rou^2))
    
    #### Equation 5.
    phi = phi_rou - (t_rou/(m0*rou_rou))*(delta_E^2/(2*m0*(niu_rou)))
    
    res <- list(latitude = phi*(180/pi), longitude = lambda*(180/pi))
    return(res)
}
