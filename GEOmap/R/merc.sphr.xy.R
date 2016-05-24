`merc.sphr.xy` <-
function(LON0, phi, lam)
  {

    DEG2RAD =pi/180
    theta = (RPMG::fmod(lam, 360)  - LON0);
    R.MAPK = MAPconstants()$R.MAPK
    ####  theta = RPMG::fmod(theta, 360) 
    x = DEG2RAD * R.MAPK * theta;
    y = R.MAPK * log(tan(DEG2RAD * (45 + phi / 2)));
    return(list(x=x, y=y))
  }

