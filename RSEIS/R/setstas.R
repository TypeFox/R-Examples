`setstas` <-
function(stafile)
{
################  read in a LAT LON Z for a station file - the z is positive km above sealevel (0)
  if(is.null(stafile)) { return(NULL) }
  
  if(!is.list(stafile) & is.character(stafile) )
    {  sta = scan(file=stafile, list(name="", lat=0, lon=0, z=0), quiet=TRUE, flush=TRUE) }
  else
    {
      sta =stafile
    }
  if(any(abs(sta$z)>100)) { print("BEWARE: STATIONS are in meters!, dividing by 1000") ; sta$z = sta$z/1000 }
  
  invisible(sta)
}

