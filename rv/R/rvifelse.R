	

rvifelse <- function (test, yes, no) {
  #   rvifelse - If-Then-Else For Random Vectors
  rvmapply(base:::ifelse, test=test, yes=yes, no=no)
}


