mogiM<-function(R=1,F=1,A=0.1,P=1e5,E=10e9,nu=0.25)
  {
##################   adapted from the matlab code of Francois Beauducel
######%MOGI   Mogi's model (point source in elastic half-space).
######%        computes radial
######%       and vertical displacements Ur and Uz, ground tilt Dt, radial and
######%       tangential strain Er and Et on surface, at a radial distance R
######%       from the top of the source due to a hydrostatic pressure inside a
######%       sphere of radius A at depth F, in a homogeneous, semi-infinite elastic
######%       body and approximation for A << F (center of dilatation). Formula by
######%       Anderson [1936] and Mogi [1958].
######%
######%       MOGI(R,F,V) and MOGI(R,F,A,mu,P) are also allowed for compatibility
######%       (Mogi's original equation considers an isotropic material with Lame's
######%       constants equal, i.e., lambda = mu, Poisson's ratio = 0.25).
######%
######%       Input variables are:
######%          F: depth of the center of the sphere from the surface,
######%          V: volumetric change of the sphere,
######%          A: radius of the sphere,
######%          P: hydrostatic pressure change in the sphere,
######%          E: elasticity (Young's modulus),
######%         nu: Poisson's ratio,
######%          mu: rigidity (Lame's constant in case of isotropic material).
######%
######%       Notes:
######%               - Equations are all vectorized, so variables R,F,V,A,mu and P are
######%                 scalar but any of them can be vector or matrix, then outputs
######%                 will be vector or matrix of the same size.
######%               - Convention: Uz > 0 = UP, f is depth so in -Z direction.
######%               - Units should be constistent, e.g.: R, F, A, Ur and Uz in m imply
######%                 V in m3; E, mu and P in Pa; Dt in rad, Er, Et and nu dimensionless.
######%

######%
######%       Author: Francois Beauducel <beauducel@ipgp.fr>
######%         Institut de Physique du Globe de Paris
######%       Created: 1997
######%       Updated: 2010-01-05
######%
######%       References:
######%         Anderson, E.M., Dynamics of the formation of cone-sheets, ring-dikes,
######%               and cauldron-subsidences, Proc. R. Soc. Edinburgh, 56, 128-157, 1936.
######%         Mogi, K., Relations between the eruptions of various volcanoes and the
######%               deformations of the ground surfaces around them, Bull. Earthquake Res.
######%               Inst. Univ. Tokyo, 36, 99-134, 1958.

######%       Copyright (c) 1997-2009, Francois Beauducel, covered by BSD License.
######%       All rights reserved.
######%

###########  converted to R by: J. M. Lees March, 2010
    r = R
    f = F
    
    if(missing(P) & missing(E) & missing(nu))
      {  ###  case 3
        v = A;
        nu = 0.25;
        y = v/pi;
        
      }
    if(missing(P) & missing(E) & !missing(nu))
      {  ###  case 4
        v = A
        nu = nu
        y = v/pi;
      }
    if(missing(P) & !missing(E) & !missing(nu))
      {  ###  case 5
        mu = E/(2*(1+nu));
        a =  A
        p = P
        nu = 0.25;
        v = A
        y = (a^3)*p/mu;
        
      }
    

    if(!missing(P) & !missing(E) & !missing(nu))
      {  ###  case 5
        mu = E/(2*(1+nu));
        a =  A
        p = P
        nu = 0.25;
        v = A
        y = (a^3)*p/mu;
        
      }
    
    et = (1-nu)*y/((f^2 + r^2)^1.5);
    ur = r*et;
    uz = f*et;
    dt = 3*et*f*r/(f^2 + r^2);
    er = dt*(f^2 - 2*r^2)/3;
    
###   radial  displacements Ur
###   vertical displacements Uz,
    ###  ground tilt Dt,
    ### radial strain Er 
### 	tangential strain  Et on surface
    return(list(ur=ur, uz=uz, dt=dt, er=er, et=et))
  }
