"kossze" <-
function (ik=3,BE=BEM) {

  #############################################################
  # 
  # TITLE:  kossze()
  # AUTHOR: SANDOR KABOS
  # DATE:   2003
  # CALLS:  
  # NEEDS:  
  # NOTES:  
  #
  #############################################################

  MI <- BE[,ik]
  MIKI <- match(kul(MI),MI)
  BE[MIKI,]
}

