"ls." <-
function () {

  #############################################################
  # 
  # TITLE:  ls.()
  # AUTHOR: FERENC (FERKO) CSILLAG (MODIFIED BY TARMO REMMEL)
  # DATE:   DECEMBER 3, 2003
  # CALLS:  N/A
  # NEEDS:  
  # NOTES:  LISTS ONLY HIDDEN OBJECTS; THOSE STARTING WITH A .
  #
  #############################################################

  # BUILD AN OBJECT LISTING ALL OBJECTS IN A WORKSPACE
  .AA <- ls(all=TRUE, pos=1)

  # RETAIN ONLY THOSE OBJECTS WHOSE FIRST CHARACTER IS A .
  .AA <- .AA[substring(.AA,1,1)=="."]

  # REMOVE COMMON HIDDEN OBJECTS FROM THE LISTING SUCH THAT
  # ONLY OBJECTS OF INTEREST ARE DISPLAYED
  .AA[.AA != ".Traceback" & .AA != ".Last" & .AA != ".Random.seed"]

}
