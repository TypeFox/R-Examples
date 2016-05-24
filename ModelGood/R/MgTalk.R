MgTalk <- function(x,y){
  if (y>100){
    if (x %in% seq(0,y,100)) message(x)
  }
  else{
    if (y>10){
      if (x %in% seq(0,y,10)) message(x)
    }
    else
      if (y>1)
        message(x)
  }
}

