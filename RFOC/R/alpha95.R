`alpha95` <-
function(az ,   iang)
  {
    ####  given azimuths and iangles (angle from z=down

   ###  az is measured from North pointing up-wards
    ###  iang is measured from the nadir (z-down)
    DEG2RAD = pi/180
    
    a = iang * DEG2RAD;
    b = (90-az) * DEG2RAD;
    x  = sin(a) * cos(b);
    y  = sin(a) * sin(b);
    z =  cos(a);
    v = cbind(x,y,z);


    ###   this is the matrix from woodcock
    KapT = t(v) %*% v

    B = length(x)*diag( 3) - KapT
    ###   from Davis p 335 
    E1 = eigen(B)
    E = eigen( KapT )

    
    Rn = sum(y )
    Re = sum(x )
    Rd = sum(z )

    N = length(x);

    Ir = 180*atan2( sqrt(Rn^2+Re^2), Rd )/pi;
    Dr = 180*atan2(Re, Rn)/pi;

    R = sqrt(Rn^2+Re^2+Rd^2)

    K = (N-1)/(N-R);
    S = 81/sqrt(K);

    Kappa = log(E$values[1]/E$values[2])/  log(E$values[2]/E$values[3])

    Alpha95 = 180*acos(1- ( (N-R)*((20^(1/(N-1)))-1)/R))/pi;
    return(list(Ir=Ir, Dr=Dr, R=R, K=K, S=S, Alph95=Alpha95, Kappa =Kappa, E=E, MAT=v))
  }

