designMPRM <-
function(daten){

kateg.zahl <- length(table(daten))
item.zahl <- ncol(daten)
  
des_core <- diag(1, nrow=(item.zahl-1)*(kateg.zahl-1), ncol=(item.zahl-1)*(kateg.zahl-1))

app <- seq((kateg.zahl-1),(item.zahl-1)*(kateg.zahl-1), by=(kateg.zahl-1))

catadd <- rep(0, ncol(des_core))

des_corS <- lapply(split(des_core,app, drop=F), function(x){
  rbind(matrix(x, ncol=(item.zahl-1)*(kateg.zahl-1), byrow=T), catadd, deparse.level=0)
})

des_corSm <- do.call(rbind,des_corS)

fixIt <- matrix(rep(diag(-1, ncol=kateg.zahl-1, nrow=kateg.zahl-1),ncol(des_corSm)/(kateg.zahl-1)), ncol=ncol(des_corSm))

rbind(des_corSm, fixIt,catadd, deparse.level=0)

}
