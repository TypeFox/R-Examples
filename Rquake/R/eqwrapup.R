eqwrapup<-function(Ldat,
                   EQ,
                   vel,
                   distwt=20, lambdareg = 0.0, verbose=FALSE )
{
  #  following an earthquake location, calculate the errors
  ##   (PUT the origin at the earthquake)

  if(missing(distwt)) distwt=20
  if(missing(lambdareg))lambdareg = 0.0
  if(missing(verbose))verbose=FALSE

  
  MLAT = EQ$lat
  MLON = EQ$lon

  proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)

  XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)

###   these should be zero -
###  but I put this in just in case the origin is modified
  b = GEOmap::GLOB.XY(EQ$lat, EQ$lon, proj)
  EQ$x = b$x
  EQ$y = b$y

  delx = EQ$x-XY$x
  dely = EQ$y-XY$y
  
  deltadis =sqrt( (delx)^2 +  (dely)^2)

  neqns = length(XY$x)

  ROWZ = matrix(ncol=4, nrow = neqns)
  
  RHS = rep(0, length=neqns)
  RHSw = rep(0, length=neqns)
  wts = rep(1, length=neqns)

  if(TRUE)
    {
      temp = Ldat$err
      dwt =  (1.0/(1. + ((deltadis^2)/(distwt^2))))
#########   this is the Lquake weighting scheme:
      wts = dwt/sqrt(temp^3);
    }
  G1 = GETpsTT(Ldat$phase, eqz=EQ$z, staz=0, delx=delx, dely=dely,  deltadis=deltadis , vel)
  
  RHS  =      ( Ldat$sec - EQ$t - G1$TT ) 
  neqns = length(RHS)
  ROWZ =  cbind(rep(1, neqns), G1$Derivs)
  

  
###   ROWZw = diag(wts) %*% ROWZ
  RHSw = wts * RHS
  if(TRUE)
    {
      dels = rep(NA, 4)
      
      
      S = svd(ROWZ)
      
      LAM = diag(S$d/(S$d^2+lambdareg^2) )
      Gdagger = S$v %*% LAM %*% t(S$u)
####  covariance of the data:
      covD = diag(Ldat$err)
      
###  covariance of the model parameters
      covB =   Gdagger %*% covD  %*% t(Gdagger)
      
#######  extract the diagonals (variances) and get sqrt
      dels = sqrt(diag(covB))
      
    }
  

  sqres = RHS^2
  ssqres = sum(  sqres )

  rms = sqrt(mean(sqres))
  meanres =  mean(RHS)
  sdres  = sqrt(var(RHS))
  sdmean = sdres/(neqns-1)
  sswres =  sum(RHSw^2)
  ndf = neqns - 2

  colnames(covB) = c("t", "x", "y", "z")

   gap = getGAP(EQ, Ldat)
  herr = sqrt( dels[2]^2 +  dels[3]^2 )

  ez = dels[4]
   distmin = min(deltadis)
  

  Q1 = "D"
  Q2 = "D"
  
  if(rms < .5 & herr <= 5.0) { Q1 = 'C'} 
  if(rms < .3 & herr <= 2.5 & ez <= 5.0)  { Q1 = 'B'}	
  if(rms < .15 & herr <= 1.0 & ez <= 2.0)  { Q1= 'A'}
  
  if(gap <= 180.0 & distmin <= 50.0)  { Q2 = 'C'}
  if(gap <= 135.0 & distmin <= max( c(EQ$z*2.0,10.0)) ) {  Q2= 'B'}
  if(gap <= 90.0 & distmin <= max( c(EQ$z,5.0)) )  { Q2 = 'A'}
  
  
  
  giveback=list(
    res = RHS,
    rms = rms,
    meanres =  meanres,
    sdres  = sdres,
    sdmean = sdmean,
    sswres = sswres,
    ndf = ndf,
    sterrx = dels[2],
    sterry = dels[3],
    sterrz =dels[4],
    sterrt = dels[1],
    cov = covB,
    lam=diag(LAM),
    gap = gap,
    herr=herr,
    distmin = distmin,
    Q1 = Q1,
    Q2 = Q2)
  
  ##  print(data.frame(giveback))
  
  if(verbose)
    {
      cat(paste("rms = ", format(rms)) , sep="\n")
      cat(paste("dx = ", format(dels[2])) , sep="\n")
      cat(paste("dy = ", format(dels[3])) , sep="\n")
      cat(paste("dz = ", format(dels[4])) , sep="\n")
      cat(paste("Quality:", Q1, Q2), sep="\n")
    }

  return(giveback)
}
