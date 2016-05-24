require(PET)
set.seed(1)

# generation of the phantom
# generation of the phantom
nP <- readline("Press 'Enter' to generate a Shepp-Logan head phantom ('nP=6'). Otherwise provide a number 'nP' between '1' and '6': ")
if(is.na(as.integer(nP))) nP <- 6 else nP <- as.integer(nP)
if(nP<1 || nP>6){
print("'nP' has to be between '1' and '6'. 'nP=6' is used. \n")
nP <- 6
}
design <- switch(nP, c(0.5, 0, 0, 0.4, 0.6),
       matrix( c(0.6, -0.35, 0, 0.4, 0.6, 0.3,  0.5,  0, 0.2, 0.35), nrow=2, byrow=TRUE),
       matrix( c(0.8, 0, 0, 0.6, 0.6, -0.4, 0, 0, 0.4, 0.2), nrow=2, byrow=TRUE),
       matrix( c(0.3, 0, -0.2, 0.5, 0.5, 0.6, 0,  0.2, 0.3, 0.5), nrow=2, byrow=TRUE),
        "B", "A")
phantomName <- paste("Phantom", nP, ".pet", sep="")
if(!(file.access(phantomName, 0))){
P <- readData(phantomName)$Signal
}else{
if(nP==6) P <- phantom(design=design,addIm="blurred1")
else P <- phantom(design=design)
try(writeData(P, phantomName, DebugLevel="HardCore"))
}
phantomName <- paste("Phantom", nP, sep="")
viewData(P,"Phantom with 257x257 samples")

# generation of the radon data
nS <- readline("Press 'Enter' to simulate Radon data with 'nS*100000' events, 'nS=10'. Otherwise provide a number 'nS' between '1' and '30': ")
if(is.na(as.integer(nS))) nS <- 10 else nS <- as.integer(nS)
if(nS>30){
print("nS larger than 30 takes to much time for a demo, nS=30 is used \n")
nS <- 30
}
namerP <- paste(phantomName, ".nS=", as.character(nS), ".pet", sep="")
nS <- nS*100000
if(!(file.access(namerP, 0))){
rP <- readData(namerP)$Signal
}else{
rP <- markPoisson(P,nSample=nS)
try(writeData(rP$rData, namerP, fileHeader=rP$Header, DebugLevel="HardCore"))
rP <- rP$rData
}
viewData(rP,paste("Simulated PET data with", nS, "counts"))

readline("Press 'Enter' return for iradon 'FB':")
irP.FB <- iradon(rP, 257, 257)
irP.FB1 <- irP.FB
viewData(irP.FB$irData,"Reconstruction with 'FB'")
readline("Press 'Enter' return for iradon 'BF':")
irP.BF <- iradon(rP, 257, 257, mode="BF")
viewData(irP.BF$irData,"Reconstruction with 'BF'")
readline("Press 'Enter' return for iradon 'CNF':")
irP.CNF <- iradon(rP, 257, 257, mode="CNF")
viewData(irP.CNF$irData,"Reconstruction with 'CNF'")
readline("Press 'Enter' return for iradon 'CBF':")
irP.CBF <- iradon(rP, 257, 257, mode="CBF")
viewData(irP.CBF$irData,"Reconstruction with 'CBF'")
readline("Press 'Enter' return for iradon 'CC':")
irP.CC <- iradon(rP, 257, 257, mode="CC")
viewData(irP.CC$irData,"Reconstruction with 'CC'")
readline("Press 'Enter' return for iradon 'CNC':")
irP.CNC <- iradon(rP, 257, 257, mode="CNC")
viewData(irP.CNC$irData,"Reconstruction with 'CNC'")
readline("Press 'Enter' return for iradon 'CBC':")
irP.CBC <- iradon(rP, 257, 257, mode="CBC")
viewData(irP.CBC$irData,"Reconstruction with 'CBC'")

readline("Press 'Enter' to compute and display the RMSE's:")
norm.all <- c(norm(P, irP.FB$irData), norm(P, irP.BF$irData), norm(P, irP.CNF$irData), norm(P, irP.CBF$irData), norm(P, irP.CC$irData), norm(P, irP.CNC$irData), norm(P, irP.CBC$irData))
names(norm.all) <- c("FB", "BF", "CNF", "CBF", "CC", "CNC", "CBC")
barplot(norm.all, col="red3", ylim=c(0,max(norm.all)+0.05), xpd=FALSE, main="Comparsion the RMSE of different reconstruction methods with 'iradon'")

readline("Press 'Enter' to smooth radon data (require 'adimpro' package):")
if(!("adimpro" %in% .packages(all.available=TRUE))){
stop("Package 'adimpro' doesn't exist")
}
require(adimpro)
degree <- readline("Press 'Enter' to smoothing the radon data using a local linear model ('degree=1'). Otherwise provide the number '2' using a lokal quadratic model: ")
if(is.na(as.integer(degree))) degree <- 1 else degree <- as.integer(degree)
if(!(degree %in% c(1,2))){
print("'degree' has to be '1' or '2'. 'degree=1' is used. \n")
nP <- 1
}
hmax <- readline("Press 'Enter' to smoothing the radon data using maximum bandwith 'hmax=9' when 'degree=1' and  'hmax=15' when 'degree=2'. Otherwise provide a number between '2' and '20': ")
if(is.na(as.double(hmax))) hmax <- switch(degree,9,15) else hmax <- as.double(hmax)
if(hmax<2 || hmax>20){
print("'hmax' has to be between '2' and '20'. 'hmax=9' when 'degree=1' and 'hmax=15' when 'degree=2' is used. \n")
nP <- switch(degree,9,15)
}
rP.part <- cutMatrix(rP)
tmp.smooth <- make.image(scaleImage(rP.part$A, mode="max"), gamma=TRUE, compress=FALSE)
tmp.smooth <- awspimage(tmp.smooth, hmax=hmax, degree=degree, varmodel="Linear", compress=FALSE, graph=TRUE)$img
rP.smooth <- matrix(0, nrow=rP.part$dimOrg[1], ncol=rP.part$dimOrg[2])
rP.smooth[, rP.part$pattern[1]:rP.part$pattern[2]] <- tmp.smooth
viewData(rP.smooth,paste("Smoothed sinogram with 'degree=",degree,"' and 'hmax=",hmax,"'.",sep=""))
rm(degree,hmax)

readline("Press 'Enter' return for iradon 'FB':")
irP.FB <- iradon(rP.smooth, 257, 257)
viewData(irP.FB$irData,"Reconstruction of the smoothed sinogram with 'FB'")
readline("Press 'Enter' return for iradon 'BF':")
irP.BF <- iradon(rP.smooth, 257, 257, mode="BF")
viewData(irP.BF$irData,"Reconstruction of the smoothed sinogram with 'BF'")
readline("Press 'Enter' return for iradon 'CNF':")
irP.CNF <- iradon(rP.smooth, 257, 257, mode="CNF")
viewData(irP.CNF$irData,"Reconstruction of the smoothed sinogram with 'CNF'")
readline("Press 'Enter' return for iradon 'CBF':")
irP.CBF <- iradon(rP.smooth, 257, 257, mode="CBF")
viewData(irP.CBF$irData,"Reconstruction of the smoothed sinogram with 'CBF'")
readline("Press 'Enter' return for iradon 'CC':")
irP.CC <- iradon(rP.smooth, 257, 257, mode="CC")
viewData(irP.CC$irData,"Reconstruction of the smoothed sinogram with 'CC'")
readline("Press 'Enter' return for iradon 'CNC':")
irP.CNC <- iradon(rP.smooth, 257, 257, mode="CNC")
viewData(irP.CNC$irData,"Reconstruction of the smoothed sinogram with 'CNC'")
readline("Press 'Enter' return for iradon 'CBC':")
irP.CBC <- iradon(rP.smooth, 257, 257, mode="CBC")
viewData(irP.CBC$irData,"Reconstruction of the smoothed sinogram with 'CBC'")

readline("Press 'Enter' to compute and display the RMSE's:")
norm.all2 <- c(norm(P, irP.FB$irData), norm(P, irP.BF$irData), norm(P, irP.CNF$irData), norm(P, irP.CBF$irData), norm(P, irP.CC$irData), norm(P, irP.CNC$irData), norm(P, irP.CBC$irData))
names(norm.all2) <- c("FB", "BF", "CNF", "CBF", "CC", "CNC", "CBC")

norm.all <- rbind(norm.all, norm.all2)
rownames(norm.all) <- c("Without smoothing", "AWS")
barplot(norm.all, beside=TRUE, col=c("red","mistyrose"), legend=rownames(norm.all), ylim=c(0, max(norm.all)+0.05), xpd=FALSE, main="Comparsion the RMSE of different reconstruction with 'iradon'")

readline("Press 'Enter' to display the 'FB' results in a new window:")
viewData(list(P, irP.FB1$irData, irP.FB$irData), list("Orginal data", "iradon with 'FB'", "AWS + iradon with 'FB'"), curWindow=FALSE)
rm(design, irP.BF, irP.CBC, irP.CBF, irP.CC, irP.CNC, irP.CNF, irP.FB, irP.FB1, nP, phantomName, namerP, norm.all, norm.all2, nS, P, rP, rP.part, rP.smooth, tmp.smooth)
