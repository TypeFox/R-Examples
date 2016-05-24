LOAD <- function(filename){



liste     <- list.files(filename)
gen       <- vector("list",length(liste))

for(xx in 1:length(liste)){

genX            <- paste(filename,"/",liste[xx],sep="")
tryl            <- try(PopGenread(genX),silent=TRUE)
# tryl            <- try(read.dna(genX,format="fasta",as.character=TRUE))

if(is.matrix(tryl)){
   print(genX)
   gen[[xx]] <- get_data(tryl)
}

}      

}
