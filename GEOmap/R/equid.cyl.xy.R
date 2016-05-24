`equid.cyl.xy` <-
function(LON0, phi0, phi, lam)
  {
    DEG2RAD=pi/180
    theta = (RPMG::fmod(lam, 360)  - LON0);
    R.MAPK  =MAPconstants()$R.MAPK
    ####  theta = RPMG::fmod(theta, 360) 
    x = DEG2RAD * R.MAPK * theta*cos(DEG2RAD *phi0)
    y =R.MAPK * DEG2RAD * phi 
    return(list(x=x, y=y))
  }

