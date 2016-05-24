require(PET)
set.seed(1)

# generation of the phantom
if(!(file.access("Phantom1.pet", 0))){
P <- readData("Phantom1.pet")$Signal
}else{
P <- phantom()
try(writeData(P, "Phantom1.pet", DebugLevel="HardCore"))
}  
nameP <- "Phantom1"
rP1 <- radon(P)$rData
viewData(rP1,"Radon transformed data")

# generation of the radon data
nS <- readline("Press 'Enter' to simulate Radon data with 'nS*100000' events, 'nS=10'. Otherwise provide a number 'nS' between '1' and '30':")
if(is.na(as.integer(nS))) nS <- 10 else nS <- as.integer(nS)

if(nS>30){
print("nS larger than 30 takes to much time for a demo, nS=30 is used \n")
nS <- 30
}
namerP2 <- paste(nameP, ".nS=", as.character(nS), ".pet", sep="")
nS <- nS*100000

if(!(file.access(namerP2, 0))){
rP2 <- readData(namerP2)$Signal
}else{
rP2 <- markPoisson(P,nSample=nS)
try(writeData(rP2$rData, namerP2, fileHeader=rP2$Header, DebugLevel="HardCore"))
rP2 <- rP2$rData
}
viewData(rP2,"Simulated PET data")

readline("Press 'Enter' return for iradon with a Hamming filter:")
irP1a <- iradon(rP1, 257, 257, FilterTyp="Hamming1")
viewData(irP1a$irData,"Reconstruction of Radon data: FilterTyp='Hamming1'")
readline("Press 'Enter' return for iradon  with a Ramp filter:")
irP1b <- iradon(rP1, 257, 257, FilterTyp="Ramp")
viewData(irP1b$irData,"Reconstruction of Radon data: FilterTyp='Ramp'")
readline("Press 'Enter' return for iradon  with a Hamming filter:")
irP2a <- iradon(rP2, 257, 257, FilterTyp="Hamming1")
viewData(irP2a$irData,"Reconstruction of PET data: FilterTyp='Hamming1'")
readline("Press 'Enter' return for iradon  with a Ramp filter:")
irP2b <- iradon(rP2, 257, 257, FilterTyp="Ramp")
viewData(irP2b$irData,"Reconstruction of PET data: FilterTyp='Ramp'")


viewData(list(irP1a$irData, irP1b$irData, irP2a$irData, irP2b$irData), list("Reconstructed Radon data: FilterTyp='Hamming1'", "Reconstructed Radon data: FilterTyp='Ramp'", "Reconstructed PET data: FilterTyp='Hamming1'", "Reconstructed PET data: FilterTyp='Ramp'"), curWindow=FALSE)
rm(irP1a,irP1b,irP2a,irP2b,nameP,namerP2,nS,P,rP1,rP2)
