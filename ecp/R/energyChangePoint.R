 getWithin <- function (alpha_, X_) 
.Call("getWithin", alpha_, X_, PACKAGE = "ecp")

 getBetween <- function (alpha_, X_, Y_) 
.Call("getBetween", alpha_, X_, Y_, PACKAGE = "ecp")

 splitPointC <- function (s_, e_, D_, min_size_) 
.Call("splitPointC", s_, e_, D_, min_size_, PACKAGE = "ecp")

getBounds <- function(n_, lvl_, eps_)
.Call("getBounds", n_, lvl_, eps_, PACKAGE = "ecp")


eFastC <- function(Z_, K_, delta_, alpha_, eps_, verbose_ )
.Call("eFastC", Z_, K_, delta_, alpha_, eps_, verbose_, PACKAGE = "ecp")
