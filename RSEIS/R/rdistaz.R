`rdistaz` <-
  function( olat, olon, tlat, tlon)
{
#### *   This subroutine will compute the distance, azimuth, and back
##### * azimuth (in degrees), given the latitude and longitude (in degrees)
#### * of an origin point and a target point.  (E+,W-; N+,S-)

  
  twopi = 2*pi;
  R.MAPK=6378.2064
  DEG2RAD=pi/180
  RAD2DEG=180/pi
  
  
  OK=0;

  Norg = length(olat)

  BIGL = list()
  

  

  for(i in 1:Norg)
    {
      
      L = list( del=rep(0, length(tlat)) , dist=rep(0, length(tlat)) , az=rep(0, length(tlat)) , baz=rep(0, length(tlat)), err=rep(0, length(tlat)) )
      

      olat[is.na(olat)] = -100
      tlat[is.na(tlat)] = -100

      olat[olat < -90. | olat > 90.] = NA
      
      tlat[tlat < -90. | tlat > 90.] = NA

####   wsame = which[  olat==tlat &  olon ==tlon ]

      wout = (olat==tlat &  olon ==tlon) | is.na(tlat)

      olon = RPMG::fmod(olon,360)
      tlon = RPMG::fmod(tlon,360)

      KLAT = olat
      KLON = olon
      MLAT = tlat[!wout]
      MLON = tlon[!wout]

      L$err[wout] = 0
      L$del[wout] = NA
      L$az[wout]  = NA
      L$baz[wout]  = NA
      L$dist[wout]  = NA


      L$del[olat==tlat &  olon ==tlon] = 0
      L$dist[olat==tlat &  olon ==tlon] = 0
########################################  

      
      L$err[!wout]=1;
      clat = 90. - KLAT;
      clon = KLON ;
####    if(clon < 0.) { clon =clon+ 360.; }
      clon = RPMG::fmod(clon,360)
      clar = DEG2RAD*clat;
      clor = DEG2RAD*clon ;
      stho = sin(clar);
      ctho = cos(clar);
      ctlat = 90. - MLAT;
      ctlon = MLON;
#### if(clon < 0.) ctlon =ctlon+ 360.;
      ctlon = RPMG::fmod(ctlon,360)
      
      ctlar = DEG2RAD*ctlat ;
      ctlor = DEG2RAD*ctlon;
      sth = sin(ctlar);
      cth = cos(ctlar);
      dph = ctlor - clor;
      sdph = sin(dph);
      cdph = cos(dph);

      
      delr = acos(stho * sth * cdph + ctho * cth);
      del = RAD2DEG* delr ;

#### /* compute forward azimuth */

####	if(sth == 0.) { azr = 0.;}
####	else { azr = atan2(sdph,stho*cth/sth-ctho*cdph);}

      azr = rep(0, length(sdph))
      azr[sth!=0.0] = atan2(sdph,stho*cth/sth-ctho*cdph);
      az = RAD2DEG*azr;
      azr = RPMG::fmod(azr,360)

###/* compute back azimuth */
      bazr = rep(0, length(sdph))

      
      bazr[stho!=0.0] = bazr = atan2(-sdph,sth*ctho/stho-cth*cdph);

      baz = RAD2DEG*bazr;
      bazr = RPMG::fmod(bazr,360)

      L$err[!wout] = 1
      L$del[!wout] = del
      L$az[!wout]  =az
      L$baz[!wout]  = baz
      L$dist[!wout]  = L$del[!wout]*2*pi*R.MAPK/360

      BIGL[[i]] = L

    }
  if(Norg==1)
    {
      return(L);
    }
  if(Norg>1)
    {
      return(BIGL); 
    }




  
}

