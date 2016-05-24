trueP <-
function(level_1, site.name){
true.p <- table(site.name, level_1)
scaled.true.p <- t(apply(true.p,1,norm<-function(x){return (x/sum(x))}))
scaled.true.p.dist <- apply(true.p,2,norm<-function(x){return (x/sum(x))})
true.p <- list(p = round(scaled.true.p,5), p_dist = round(scaled.true.p,5))
cat("** Warning **:", "the observed relative abundances/species distributions are calculated from the trait data!", "\n", "This may not be the correct value. Please refer to the documentation for details!", "\n")

return(true.p)
}
