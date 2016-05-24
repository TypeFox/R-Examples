"dups" <-
function (intab=KI) { 

  #############################################################
  # 
  # TITLE:  dups()
  # AUTHOR: TARMO REMMEL
  # DATE:   7 DEC, 2004
  # CALLS:  N/A
  # NEEDS:  KI OBJECT FROM HDECO
  # NOTES:  CHECKS FOR DUPLICATE RECORDS IN KI MATRIX FOR HDECO
  #         IF OMITX1 IS SPECIFIED, WHERE DUPLICATE RECORDS MIGHT
  #         EXIST
  #
  #############################################################

  # CREATE A MATRIX WITH ONLY THE UNIQUE ROWS OF intab  
  intabunique <- unique(intab)
    
  # COLLAPSE THE FULL TABLE
  a2 <- apply(intab, 1, paste, collapse=":")
  
  # COLLAPSE THE UNIQUE TABLE
  b2 <- apply(intabunique, 1, paste, collapse=":")
  
  # COUNT UNIQUE ROW OCCURRENCES AND WRITE TO .COUNT
  assign(".COUNT",c(table(c(a2,unique(b2)))[b2] - 1) ,pos=1)
  
  # WRITE THE UNIQUE TABLE TO .KIUNIQUE
  assign(".KIUNIQUE", intabunique, pos=1)   
}

