`DetAR` <-
function(phi){
z<-ARToPacf(phi)
1/prod((1-z^2)^(1:length(phi)))
}

