`InvertibleQ` <-
function(phi){
    identical(TRUE,try(all(abs(ARToPacf(phi))<1),silent=TRUE))
}
