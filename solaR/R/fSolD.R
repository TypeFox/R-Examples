fSolD<-function(lat, BTd, method='michalsky'){

  if (abs(lat)>90){
    lat=sign(lat)*90
    warning(paste('Latitude outside acceptable values. Set to', lat))
  }
  
  lat=d2r(lat)

  if (missing(BTd)){
    BTd=fBTd(mode='prom')
  } else {
    BTd=unique(truncDay(BTd))
  }
  
  dn <- doy(BTd)                         #día del año
  jd <- as.numeric(julian(BTd + 12*3600, ##12:00:00 UTC
                          origin='2000-01-01 12:00:00 UTC'))
  X = 2*pi*(dn-1)/365

  methods=c('cooper', 'spencer', 'michalsky', 'strous')
  method=match.arg(method, methods)

  ##Declination
  decl=switch(method,
    cooper={
      ##Cooper, P.I., Solar Energy, 12, 3 (1969). "The Absorption of Solar Radiation in Solar Stills"
      decl=23.45*sin(2*pi*(dn+284)/365) 
      decl=d2r(decl)
    },
    spencer={
      ##Spencer, Search 2 (5), 172
      ##http://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html
      decl = 0.006918 - 0.399912*cos(X) + 0.070257*sin(X) -
        0.006758*cos(2*X) + 0.000907*sin(2*X) -
          0.002697*cos(3*X) + 0.001480*sin(3*X) 
    }, 
    strous={
      meanAnomaly = (357.5291 + 0.98560028*jd)%%360
      coefC=c(1.9148, 0.02, 0.0003)
      sinC=sin(outer(1:3, d2r(meanAnomaly), '*'))
      C = colSums(coefC*sinC)
      trueAnomaly=(meanAnomaly + C)%%360
      eclipLong=(trueAnomaly + 282.9372)%%360
      excen=23.435
      sinEclip=sin(d2r(eclipLong))
      sinExcen=sin(d2r(excen))
      decl=asin(sinEclip*sinExcen)
    }, 
    michalsky={
      meanLong=(280.460+0.9856474*jd)%%360
      meanAnomaly=(357.528+0.9856003*jd)%%360
      eclipLong=(meanLong +1.915*sin(d2r(meanAnomaly))+0.02*sin(d2r(2*meanAnomaly)))%%360
      excen=23.439-0.0000004*jd
      sinEclip=sin(d2r(eclipLong))
      sinExcen=sin(d2r(excen))
      decl=(asin(sinEclip*sinExcen))
    }
    )

  ##Distancia sol-tierra, 1/r2
  ##ro=1.496E8                        #distancia media Tierra-Sol (km)
  eo = switch(method,
    cooper = 1 + 0.033*cos(2*pi*dn/365),
    spencer = 1.000110 + 0.034221*cos(X) + 0.001280*sin(X) + 0.000719*cos(2*X) + 0.000077*sin(2*X),
    michalsky = 1.000110 + 0.034221*cos(X) + 0.001280*sin(X) + 0.000719*cos(2*X) + 0.000077*sin(2*X),
    strous = 1.000110 + 0.034221*cos(X) + 0.001280*sin(X) + 0.000719*cos(2*X) + 0.000077*sin(2*X)
    )


  ##Equation of Time, minutes
  ##según Alan M.Whitman "A simple expression for the equation of time"
  ##EoT=ts-t, donde ts es la hora solar real y t es la hora solar media
  ##Valores negativos implican que el sol real se retrasa respecto al medio
  M=2*pi/365.24*dn
  EoT.min=229.18*(-0.0334*sin(M)+0.04184*sin(2*M+3.5884))
  EoT=h2r(EoT.min/60)                   #radianes

  ##daylength
  cosWs=-tan(lat)*tan(decl)
  ws=suppressWarnings(-acos(cosWs)) #sunrise, negative since it is before noon
  polar <- which(is.nan(ws))        ##Polar day/night
  ws[polar] <- -pi*(cosWs[polar] < -1) + 0*(cosWs[polar] >1)
  
  ##Extraterrestial irradiance
  Bo=1367                               #solar constant
  Bo0d=-24/pi*Bo*eo*(ws*sin(lat)*sin(decl)+cos(lat)*cos(decl)*sin(ws)) #el signo negativo se debe a la definición de ws

  result<-zoo(data.frame(decl=decl, eo=eo, ws=ws, Bo0d=Bo0d, EoT=EoT), truncDay(BTd))
  attr(result, 'lat')=r2d(lat)
  result
}
