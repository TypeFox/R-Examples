`foc.icolor` <-
function(rake)
  {
#### /* strikeslip */
    if(abs(rake) <= 15.0) return(1);
    
####  /* strikeslip */
    if(abs((180.0 - abs(rake))) <= 15.0) return(1);

####  /* rev-obl strk-slp */
    if((rake >= 15.0 & rake < 45) || (rake >= 135 && rake < 165))
      return(2);

####  /* oblique reverse */
  if((rake >= 45.0 & rake < 75) || (rake >= 105 && rake < 135))
    return(3);

####  /* reverse */
  if(rake >= 75.0 & rake < 105.0) return(4);

####  /* norm-oblq strkslp */
  if((rake < -15.0 & rake >= -45) || (rake < -135 && rake >= -165))
    return(5);

####  /* oblq norm */
  if((rake < -45.0 & rake >= -75) || (rake < -105 && rake >= -135))
    return(6);

####  /* normal */
  if(rake < -75.0 & rake >= -105)  return(7);



  }

