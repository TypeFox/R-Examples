`makenet` <-
function()
{
###########   MN = make.net()
                                        #% SMALL RPMG::circles  (Latitudes)
  alons = pi*seq(from=0, to=180, by=5)/180
  alats = seq(from=-80, to=80, by=10)
  x1 =  matrix(nrow= (length(alons)+1),  ncol=(length(alats)))
  y1 =  matrix(nrow= (length(alons)+1),  ncol=(length(alats)))

  lam0 = pi/2
  for(i in 1:length(alats) )
    {
      phi = alats[i]*pi/180
      R = sqrt(2)/2
      kp = sqrt(2/(1+cos(phi)*cos(alons-lam0)))
      x1[,i] =  c(R*kp*cos(phi)*sin(alons-lam0), NA)
      y1[,i] =  c(R * kp*sin(phi), NA)
    }
  alons =seq(from=10, to=170, by=10)
  alats = seq(from=-90, to=90, by=5)*pi/180

  x2 =  matrix(nrow=length(alats)+1 , ncol= length(alons) )
  y2 =  matrix(nrow=length(alats)+1 ,  ncol=length(alons))

  for(i in 1:length(alons) )
    {
      lam = alons[i]*pi/180
      R = sqrt(2)/2
      kp = sqrt(2/(1+cos(alats)*cos(lam-lam0)))
      x2[,i] = c(R*kp*cos(alats)*sin(lam-lam0), NA)
      y2[,i] = c(R * kp*sin(alats), NA)
    }
  return(list(x1=x1, y1=y1,x2=x2, y2=y2 ))
}

