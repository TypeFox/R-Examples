make.cohort <-
function(df,ce){
#  make.cohort()  -  setup a cohort vector
#                 - ce is vector of df col nos for parts of cohort env code to be
#                    pasted together to make a single cohort env code
  ne <- length(ce)
  ni <- length(df$Id)
  cohort <- rep(0,ni)
  for(i in 1:ni) {
    noce <- 0
    for(k in 1:ne) {
      if( is.na(df[[ce[k]]][i])) {
        noce <- 1
      }
    }
    if(noce == 1) {
      cohort[i] <- NA
    }
    else {
      cohort[i] <- paste(df[[ce[1]]][i])
      if(ne > 1) {
        for(k in 2:ne) {
          cohort[i] <- paste(cohort[i],df[[ce[k]]][i],sep="")
        }
      }
    }
  }
  return(cohort) 
}
