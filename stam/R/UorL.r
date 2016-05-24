UorL<-function (x,charlower=TRUE) {
  #=======================================================
  #
  #  TITLE:     Converting all the letters of Character variable into uppercase or lowercase
  #  FUNCTION:  uORl()
  #  AUTHOR:    Zhijie Zhang
  #  DATE:      11 JANUARY 2010
  #  CALLS:     
  #  NEEDS:
  #  NOTES: It defaults to change all the letters into lowercase.
  #  x-Character variable whose letters to be changed into uppercase or lowercase.
  #  charlower- Specify whether you want to get the uppercase letters or lowercase letters. 
  #  CHANGE HISTORY:
  #=======================================================
 ifelse(charlower==TRUE,tolower(as.character(x)),toupper(as.character(x)))
}