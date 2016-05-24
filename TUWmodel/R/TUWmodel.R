TUWmodel <- function (prec, airt, ep, area=1, param=c(1.2,1.2,2,-2,0,0.9,100,3.3,0.5,9,105,50,2,10,26.5), incon=c(50,0,2.5,2.5), itsteps=NULL) {
 nzones <- ifelse(is.vector(prec), 1, dim(prec)[2])
 itsteps <- ifelse(is.null(itsteps), length(prec)/nzones, itsteps)
 if (nzones == 1) {
  parametri <- as.matrix(t(param))
  inconditions <- as.matrix(t(incon))
 } else if (nzones < 1) {
  cat("\nFormatting harddisk. Please smile!\n")
 } else {
  if (is.matrix(param)) {
   parametri <- param
  } else if (is.vector(param)) {
   parametri <- matrix(rep(param, nzones), ncol=nzones)
  }
  if (is.matrix(incon)) {
   inconditions <- incon
  } else if (is.vector(incon)) {
   inconditions <- matrix(rep(incon, nzones), ncol=nzones)
  }
 }
 storage.mode(prec) <- "double"
 storage.mode(airt) <- "double"
 storage.mode(ep) <- "double"
 storage.mode(area) <- "double"
 storage.mode(parametri) <- "double"
 storage.mode(inconditions) <- "double"
 output <- array(-777, dim=c(nzones, 20, itsteps))
 storage.mode(output) <- "double"

 dummy <- .Fortran("hbvmodel", itsteps=as.integer(itsteps), nzones=as.integer(nzones), area=area, param=parametri, incon=inconditions,
                    prec=prec, airt=airt, ep=ep, output=output, PACKAGE="TUWmodel") 
 names(dummy$param) <- c("SCF","DDF","Tr","Ts","Tm","LPrat","FC","BETA","k0","k1","k2","lsuz","cperc","bmax","croute")
 names(dummy$incon) <- c("SSM0","SWE0","SUZ0","SLZ0")
 dummy$qzones <- t(dummy$output[,1,])
 if (nzones > 1) {dummy$q <- apply(dummy$qzones,1,sum)} else {dummy$q <- dummy$qzones}
 dummy$swe <- t(dummy$output[,2,])
 dummy$melt <- t(dummy$output[,6,])
 dummy$q0 <- t(dummy$output[,7,])
 dummy$q1 <- t(dummy$output[,8,])
 dummy$q2 <- t(dummy$output[,9,])
 dummy$moist <- t(dummy$output[,3,])
 dummy$rain <- t(dummy$output[,4,])
 dummy$snow <- t(dummy$output[,5,])
 dummy$eta <- t(dummy$output[,10,])
 dummy$suz <- t(dummy$output[,11,])
 dummy$slz <- t(dummy$output[,12,])

 return(dummy)
}







