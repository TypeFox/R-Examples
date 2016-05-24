is.fdata<-function(fdataobj){
if (!inherits(fdataobj, "fdata"))     return(FALSE)
else {
type<-switch(class(fdataobj),
matrix={if (ncol(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)},
data.frame={if (ncol(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)},
numeric={if (length(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)},
integer={if (length(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)})
}
return(TRUE)
}
