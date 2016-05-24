
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
write8lyse <- function(lyse,fields,check)
#TITLE ancillary function 
#DESCRIPTION to display the contents of an object
# (The arguments are not described because it is
#  a masked function.)
#FUTURE improve the formatting of the content
# (removing undue break lines and formatting the lists)
#CREATED 14_07_02
#REVISED 14_07_02
#--------------------------------------------
{
  # checking
  rrrobject9(lyse,"list");
  rrrobject9(fields,"character");
  # displaying
  cat("\n### content of",lyse$name,"###\n\n");
  if (length(lyse)>0) {
    nbf <- 0; chvi <- character(0);
    for (ily in rrrbf(lyse)) {
      quel <- names(lyse)[ily];
      vale <- lyse[[ily]];
      if (quel %in% fields) {
        if (length(vale)==0) {
          chvi <- c(chvi,quel);
        } else {
          nbf <- nbf+1;
          cat("-<",quel,">- : \n",sep="");
          rrrform3parag(lyse[[ily]],titre=-2);
        }
      }
    }
    if (nbf==0) {
      cat("There were fields but no one was required...\n");
    } else {
      if (length(chvi)>0) {
        cat("The following field(s) were empty:",
            paste(chvi,collapse="/"),"\n");
      }
    }
  } else {
    cat("No one field was provided for this documented object!\n");
  }
  # monitoring
  if (check) { rrrpause(paste0("Content of ",lyse$name," right?"));}
  # returning
  invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
