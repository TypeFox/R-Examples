`read.Data` <-
function(blocksperarray=4,spotter="arrayjet",writetable=FALSE,printFlags=FALSE,
        fileName="Flagged_spots.csv", remove_flagged=NULL, ...){

    ## read gprfiles
    temp.dat <-  read.gpr(blocksperarray=blocksperarray,spotter=spotter, remove_flagged=remove_flagged, ...)
    ## annotate arrays
    incubation <- read.slidedescription()
    temp.dat.i <- annotate.arrays (temp.dat,incubation)
    ## annotate samples
    sampleID <- read.sampledescription()
    temp.dat.ii <- annotate.samples(temp.dat.i,sampleID)
    ## export txt file
    if(writetable){
        write.Data(temp.dat.ii)
    }
    ## print flagged spots


    if (printFlags){
        for (i in 1:ncol(temp.dat.ii[[1]])){

            flines <- temp.dat.ii[[5]][,i]==-100
            if ( any (flines)){
                tempflags <- temp.dat.ii[[4]][flines,]

                write.csv(tempflags,file=paste(temp.dat.ii[[3]]["target",i],"_",
                                temp.dat.ii[[3]]["AB_ID",i],"_",
                                fileName,sep=""))
            }
        }
    }


    return(temp.dat.ii[1:4])
}

