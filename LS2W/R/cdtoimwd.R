cdtoimwd <-
function (cddews) 
{
data <- cddews$S
struct <- cddews$structure
nlev <- cddews$nlevels
datadim <- cddews$datadim
family <- cddews$family
filter.number <- cddews$filter.number

zero.mat <- matrix(0,datadim[1],datadim[2])
z.imwd <- imwd(zero.mat, filter.number, family, type="station")

for( i in (nlev-1):0){

tmpH<-lt.to.name(i,"CD")
tmpV<-lt.to.name(i,"DC")
tmpD<-lt.to.name(i,"DD")

z.imwd[[tmpV]] <- as.vector(data[nlev-i,,])
z.imwd[[tmpH]] <- as.vector(data[nlev-i+nlev,,])
z.imwd[[tmpD]] <- as.vector(data[nlev-i+2*nlev,,])
}
return(z.imwd)

}

