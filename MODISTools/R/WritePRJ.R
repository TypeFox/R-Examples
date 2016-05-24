WritePRJ <-
function(Path)
{
  prj <- paste('PROJCS["Sinusoidal",GEOGCS["GCS_Undefined",DATUM["Undefined",',
               'SPHEROID["User_Defined_Spheroid",6371007.181,0.0]],PRIMEM["Greenwich",0.0],',
               'UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],',
               'PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],',
               'PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]',
               sep = "")
  write(prj, file = Path)
}