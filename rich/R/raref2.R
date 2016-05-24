raref2 <-
function(matrix, dens, tolerance, nrandom=99)
 {
  thinning<-function(comD=comD,d,ub=ub,lb=lb)
  {
    it<-1;sss<-sum(comD)
    co<-comD[d,] 
    while(sss>ub){
    repeat{
    if(is.null(dim(co)[1])==TRUE || sum(co)<ub) {break}
    x<-sample(seq(1:dim(co)[1]),1)
    co<-co[-x,]
    sss<-sum(co)}
    if(sum(co)<lb) {co<-comD;sss<-sum(co)} 
    if(is.null(dim(co)[1])==TRUE) {co<-comD;sss<-sum(co)} 
    if(is.null(dim(co)[1])==FALSE && dim(co)[1]<2) {co<-comD;sss<-sum(co)} 
    if(it>100) {break}#{print("quit cause no solution");break}
    it<-it+1}
    return(specpool(co)$Species)
  }

ef<-as.matrix(matrix);dens<-dens;tolerance<-tolerance
if(is.null(nrandom)==TRUE | nrandom<10) nrandom<-99
if(is.null(dens)==TRUE | dens<=0) stop("invalid dens value")
if (any(is.na(ef))) stop("na entries in table")
if(dens>=sum(ef)) stop("invalid dens value: too high")
seuil<-dens
lb<-seuil-seuil*tolerance
ub<-seuil+seuil*tolerance
cat("computing...", "\n")
ef_prime.boot<-boot(ef, thinning, R=nrandom,ub=ub,lb=lb)
cat("done", "\n")
return(list(mean.boot=mean(ef_prime.boot$t),sd.boot=sd(as.vector(ef_prime.boot$t))))}