disp <-
function (x) {
H <- herf(x)
H.norm <- herf(x, norm=1)
H.eq <- herf.eq(x)
G <- gini(x)
G.norm <- gini(x, norm=1)
dispvalues <- c(H, H.norm, H.eq, G, G.norm)
names(dispvalues) <- c("HHI","HHI*","HHI_eq","GINI","GINI*")
return (dispvalues)
}
