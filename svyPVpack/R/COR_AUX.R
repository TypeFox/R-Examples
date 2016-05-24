########################################################################
# COR AUX FUNCTIONS
########################################################################
cor.aux.func <- function(x, y, svydes){
  withReplicates(svydes, function(w,data)
    cov.wt(cbind(x, y), 
           wt=w, cor =TRUE)$cor[1,2])}