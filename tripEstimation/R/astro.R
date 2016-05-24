`astro` <-
function(lon, lat, astro.calc) {
	RA <- astro.calc$RA
	DELTA <- astro.calc$DEC
	MODJD <- astro.calc$MJD

	## beware I'm moving to cartographic longitude
	TAU <- 15.0 * (LMST(MODJD, -lon) - RA)
	EQUHOR(DELTA, TAU, lat)

}

`EQUHOR` <-
function(DEC,TAU,PHI) {

    DEC <- DEC * (pi/180)
    TAU <- TAU * (pi/180)
    PHI <- PHI * (pi/180)


    CS_PHI = cos(PHI); SN_PHI= sin(PHI);
    CS_DEC = cos(DEC); SN_DEC =sin(DEC); CS_TAU=cos(TAU);

    X =CS_DEC*SN_PHI*CS_TAU - SN_DEC*CS_PHI;
    Y =CS_DEC * sin(TAU);
    Z =CS_DEC*CS_PHI*CS_TAU + SN_DEC*SN_PHI;

   # print(list(X, Y, Z))
    HAZ <- POLAR(X,Y,Z)
    #list(H = HAZ$theta, AZ = HAZ$phi)
    HAZ
}

`FRAC` <-
function(x) {
	## return the fractional portion, this is not documented in the reference
	##   just relies on (assuming) how Pascal's TRUNC behaves
	x <- x - trunc(x)
	#if (x < 0) x <- x + 1
	x
}

`julcent` <-
function(time) {
	## julian centuries
	 tm <- as.POSIXlt(time, tz = "GMT")
	 hh <- tm$hour
	 mm <- tm$min
	 ss <- tm$sec
	 jday <- julday(tm) + (hh + (mm + ss/60)/60)/24
	(jday - 2451545)/36525
}

`LMST` <-
function(MJDay, LAMBDA) {
    ## (* LMST: local mean sidereal time
    MJD0 = trunc(MJDay);

    UT = (MJDay-MJD0)*24;
    time =(MJD0-51544.5)/36525.0;
    GMST =6.697374558 + 1.0027379093*UT +
    	(8640184.812866+(0.093104 - 6.2E-6*time)*time)*time/3600.0;
    24.0*FRAC( (GMST-LAMBDA/15.0) / 24.0 );

}

`lunar` <-
function(time) {
   ARC = 206264.8062;
   P2  = 6.283185307;
   COSEPS = 0.91748;
   SINEPS = 0.39778;
   stime <- julcent(time)
#    (* mean elements of lunar orbit *)
    L0 =   FRAC(0.606433+1336.855225*stime); #(* mean longitude Moon (in rev) *)
    L = P2*FRAC(0.374897+1325.552410*stime); #(* mean anomaly of the Moon     *)
    LS =P2*FRAC(0.993133+  99.997361*stime); #(* mean anomaly of the Sun      *)
    dd = P2*FRAC(0.827361+1236.853086*stime); #(* diff. longitude Moon-Sun     *)
    ff = P2*FRAC(0.259086+1342.227825*stime); #(* mean argument of latitude    *)
    DL = 22640*sin(L) - 4586*sin(L-2*dd) + 2370*sin(2*dd) +  769*sin(2*L)  -
          668*sin(LS)- 412*sin(2*ff) - 212*sin(2*L-2*dd) - 206*sin(L+LS-2*dd) +
          192*sin(L+2*dd) - 165*sin(LS-2*dd) - 125*sin(dd) - 110*sin(L+LS) +
          148*sin(L-LS) - 55*sin(2*ff-2*dd);
    S = ff + (DL+412*sin(2*ff)+541*sin(LS)) / ARC;
    H = ff - 2*dd;
    N = -526*sin(H) + 44*sin(L+H) - 31*sin(-L+H) - 23*sin(LS+H) +
          11*sin(-LS+H) -25*sin(-2*L+ff) + 21*sin(-L+ff);
    L_MOON = P2 * FRAC ( L0 + DL/1296E3 ); #(* in rad *)
    B_MOON = ( 18520.0* sin(S) + N ) / ARC; #(* in rad *)
    #(* equatorial coordinates *)
    CB = cos(B_MOON);
    X = CB*cos(L_MOON); V =CB*sin(L_MOON); W =sin(B_MOON);
    Y =COSEPS*V-SINEPS*W; Z=SINEPS*V+COSEPS*W; RHO=sqrt(1.0-Z*Z);
    DEC = (360.0/P2)*atan(Z/RHO);
    RA  = ( 48.0/P2)*atan(Y/(X+RHO)); ifelse(RA<0, RA+24.0, RA)

   MODJD <- MJD(time)
    list(RA = RA, DEC = DEC, MJD = MODJD)
}

`mini.sun` <-
function(time) {

   P2  = 6.283185307;
   COSEPS = 0.91748;
   SINEPS = 0.39778;
    stime <- julcent(time)
    M  = P2*FRAC(0.993133+99.997361*stime);
    DL = 6893.0*sin(M)+72.0*sin(2*M);
    L  = P2*FRAC(0.7859453 + M/P2 + (6191.2*stime+DL)/1296E3);
    SL = sin(L);
    X =cos(L);
    Y = COSEPS*SL;


    Z = SINEPS*SL;
    RHO = sqrt(1.0-Z*Z);
    DEC = (360.0/P2)*atan(Z/RHO);

    RA  = ( 48.0/P2)*atan(Y/(X+RHO));
    RA <- ifelse(RA < 0, RA+24.0, RA);

    MODJD <- MJD(time)
    list(RA = RA, DEC = DEC, MJD = MODJD)
}

`MJD` <-
function(date) {
    ## modified julian date (it's a smaller number)
    date <- as.POSIXlt(date, tz = "GMT")
    YEAR <- as.numeric(format(date, "%Y"))
    MONTH <- as.numeric(format(date, "%m"))
    DAY <- as.numeric(format(date, "%d"))
    HOUR <- date$hour + date$min/60 + date$sec/3600
    A =10000.0*YEAR+100.0*MONTH+DAY;
    YEAR = ifelse(MONTH <= 2, YEAR-1, YEAR)
    MONTH = ifelse(MONTH <= 2, MONTH + 12, MONTH);

    B = ifelse(A<=15821004.1, -2 + trunc((YEAR+4716)/4)-1179,
                                  trunc(YEAR/400)-trunc(YEAR/100)+trunc(YEAR/4))
    A =365.0*YEAR-679004.0;
    A+B+trunc(30.6001*(MONTH+1))+DAY+HOUR/24.0;
}

`POLAR` <-
function(X,Y,Z) {
    RHO =X*X+Y*Y;  R = sqrt(RHO+Z*Z)
    PHI <- atan2(Y,X)*180/pi ; PHI <- ifelse(PHI<0, PHI+360.0, PHI);

    ##Thanks to Nick.Ellis@csiro.au for the report.  2009-05-29
    THETA <- atan2(Z, sqrt(RHO)) * 180/pi # this line had error, was RHO now sqrt(RHO)

    list(r = R, theta= THETA, phi = PHI)
}

