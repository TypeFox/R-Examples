`fixUWstasLL` <-
function(STAS, stafile)
{
  if(length(STAS$name)<1) return(NULL) 
  sta =setstas(stafile)
  stmatch = match(STAS$name, sta$name)			
  STAS$lat = sta$lat[stmatch]	
  STAS$lon = sta$lon[stmatch]	
  STAS$z = sta$z[stmatch]	
  
  return(STAS) 
}

