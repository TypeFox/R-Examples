"bits" <-
function(object,bit) {
  (object %/% (2^bit)) %% 2
}

"bits<-" <-
function(object,bit,value) {
  mask <- 2^bit
  object <- object+(value - ((object %/% mask) %% 2))*mask
  object
}


