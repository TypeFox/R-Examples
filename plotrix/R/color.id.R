color.id<-function(col) {
 c2 <- col2rgb(col)
 coltab <- col2rgb(colors())
 cdist <- apply(coltab, 2, function(z) sum((z - c2)^2))
 colors()[which(cdist == min(cdist))]
}
