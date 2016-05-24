require(PET)
set.seed(1)
if(file.access(getwd(),2))
stop("You don't have a write permission in your working directory. It is not possible to create ini-file.'")

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
if(nP==6) P <- phantom(design=design,addIm="blurred1") else P <- phantom(design=design)
writeData(P, phantomName, DebugLevel="HardCore")
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
writeData(rP$rData, namerP, fileHeader=rP$Header, DebugLevel="HardCore")
rP <- rP$rData
}
viewData(rP, paste("Simulated PET data with", nS, "counts"))
readline("Press 'Enter' to create the first ini-file:")
ini <- file("PET.FB.ini", "wb")
writeLines("mode=FB", con=ini)
writeLines(paste("rData=",namerP, sep=""), con=ini)
writeLines("XSamples=257" ,con=ini)
writeLines("YSamples=257" ,con=ini)
writeLines("FilterTyp=Hamming2" ,con=ini)
close(ini)
readline("Press 'Enter' return for iradon (mode='FB'):")
irP.FB <- iradon(iniFile="PET.FB.ini")
writeData(irP.FB$irData, "irP.FB.pet", fileHeader=irP.FB$Header)
readline("Press 'Enter' to create the second ini-file:")
ini <- file("PET.EM.ini", "wb")
writeLines("mode=EM", con=ini)
writeLines(paste("rData=",namerP, sep=""), con=ini)
writeLines("StartImage=irP.FB.pet", con=ini)
writeLines("RadonKernel=RL", con=ini)
writeLines("Iterations=6", con=ini)
writeLines("ConstrainMin=0", con=ini)
writeLines("ConstrainMin=50", con=ini)
close(ini)
readline("Press 'Enter' return for iradonIT (mode='EM'):")
irP.EM <- iradonIT(iniFile="PET.EM.ini")
readline("Press 'Enter' to display the results in a new window:")
viewData(list(irP.FB$irData,irP.EM$irData), 
         list("iradon (mode='FB')","iradon (mode='EM')"),
         curWindow=FALSE)
rm(design,ini,irP.EM,irP.FB,namerP,nP,nS,P,phantomName,rP)
file.remove("irP.FB.pet")
