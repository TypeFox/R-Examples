mogi1 <-
function(d=1, f=1, a=0.1, P=1e5, mu=4e+09, nu=0.25)
  {

    if(missing(d)) { d = 1}
    if(missing(f)) {f=1}
    if(missing(a)) {a=0.1}
    if(missing(P)) {P=1e5}
    if(missing(mu)) {mu=4e+09}
    if(missing(nu)) { nu = 0.25 }
    
  ###  mu = E/(2*(1+nu))

 ###  E, mu  and P in Pa

    
 ###    %	   E: elasticity (Young's modulus),
 ###%	  nu: Poisson's ratio,
 ###%	   mu: rigidity (Lame's constant in case of isotropic material).

###- Units should be constistent, e.g.: R, F, A, Ur and Uz in m imply
###		  V in m3; E, mu and P in Pa; Dt in rad, Er, Et and nu dimensionless.
    

    ##############   five parameters mogi source calculation
    #######################################
    ###  a = radius of sphere injected
    ### P = hydrostatic pressure of injection
    ###  d = distance along surface
    ###  mu =shear modulus
    ###  f = depth to source
    
 ###  DELTAd = radial displacement
 ###  DELTAh = vertical  displacement
 ### %          F: depth of the center of the sphere from the surface,
 ### %          V: volumetric change of the sphere,
 ### %          A: radius of the sphere,
 ### %          P: hydrostatic pressure change in the sphere,
 ### %          E: elasticity (Young's modulus),
 ### %         nu: Poisson's ratio,
 ### %          mu: rigidity (Lame's constant in case of isotropic material).

 ###  %         Mogi, K., Relations between the eruptions of various volcanoes and the
 ### %               deformations of the ground surfaces around them, Bull. Earthquake Res.
 ### %               Inst. Univ. Tokyo, 36, 99-134, 1958.
 ### %MOGI   Mogi's model (point source in elastic half-space).
 ### %       [Ur,Uz,Dt,Er,Et] = MOGI(R,F,V,nu) or MOGI(R,F,A,P,E,nu) computes radial
 ### %       and vertical displacements Ur and Uz, ground tilt Dt, radial and
 ### %       tangential strain Er and Et on surface, at a radial distance R
 ### %       from the top of the source due to a hydrostatic pressure inside a
 ### %       sphere of radius A at depth F, in a homogeneous, semi-infinite elastic
 ### %       body and approximation for A << F (center of dilatation). Formula by
 ### %       Anderson [1936] and Mogi [1958].
 ### %
 ### %       MOGI(R,F,V) and MOGI(R,F,A,mu,P) are also allowed for compatibility
 ### %       (Mogi's original equation considers an isotropic material with Lame's
 ### %       constants equal, i.e., lambda = mu, Poisson's ratio = 0.25).
 ### %%	  Anderson, E.M., Dynamics of the formation of cone-sheets, ring-dikes,
 ### %%		and cauldron-subsidences, Proc. R. Soc. Edinburgh, 56, 128-157,	1936.
 ### %%	  Mogi, K., Relations between the eruptions of various volcanoes and the
 ### %%		deformations of the ground surfaces around them, Bull. Earthquake Res.
 ### %%		Inst. Univ. Tokyo, 36, 99-134, 1958.
   
    denom = (4*mu*(f^2 + d^2)^(1.5))
    DELTAd =  3*a^3*P*d/denom
    DELTAh =  3*a^3*P*f/denom

   return(list(ur=DELTAd, uz=DELTAh))
  }

