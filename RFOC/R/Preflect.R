`Preflect` <-
function(az, dip )
{
    #  /* Reflect to lower hemisphere */
    if (dip < 0.0)
      {
	az = RPMG::fmod(az + 180., 360.);
	dip = -dip;
    }
    return(list(az=az, dip=dip));
}

