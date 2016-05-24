internalTalk <- function(x,y,sign="'"){
  if (y>100){
    if (y<500){
      if (x %in% seq(0,y,10)) message(x)
    }
    else{
      if (x %in% seq(0,y,100)) message(x)
    }
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

