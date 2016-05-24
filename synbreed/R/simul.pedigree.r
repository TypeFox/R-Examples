simul.pedigree <- function(generations=2,ids=4,animals=FALSE,familySize=1){

        # if only one value is given, set samenumber for all generations
        if(length(ids)==1) ids <- rep(ids,times=generations)

        # initialisation
        gener <- rep(1:generations,times=ids)
        ID <- 1:sum(ids)
        Par1 <- Par2 <- rep(0,length(ID))

        # random mating for plants (inbreeds are likely)
        if(!animals){
          for (i in 2:generations){
            Par1[gener==i] <- ID[sample(ID[gener==i-1],size=ids[i],replace=TRUE)]
            Par2[gener==i] <- ID[sample(ID[gener==i-1],size=ids[i],replace=TRUE)]
          }
          ped <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener-1,sex=NA)
        }
        # define sire and dams for animals (no inbreeds)
        else{
           # define sex for 1st generation
           # 0 = female
           # 1 = male
           sex <- rep(0,length(ID))
           sex[gener==1] <- sample(rep(c(0,1),length=sum(gener==1)),sum(gener==1),replace=FALSE)

           for (i in 2:generations){
            sex[gener==i] <- sample(rep(c(0,1),length=sum(gener==i)),sum(gener==i),replace=FALSE)

            Par1[gener==i] <- ID[sample(ID[(gener==i-1) & sex==1],size=ids[i],replace=TRUE)]
            Par2[gener==i] <- ID[sample(ID[(gener==i-1) & sex==0],size=ids[i],replace=TRUE)]
          }
           ped <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener-1,sex=sex)
        }

        # create full-sib families in the last generation
        pedTemp <- ped[ped$gener==generations-1,]
        ped <- ped[ped$gener<generations-1,]
        pedFamily <- data.frame(ID=rep(1:(nrow(pedTemp)*familySize)+pedTemp$ID[1]-1),
                                Par1=rep(pedTemp$Par1,each=familySize),
                                Par2=rep(pedTemp$Par2,each=familySize),
                                gener=rep(pedTemp$gener,each=familySize),
                                sex=rep(pedTemp$sex,each=familySize))
        ped <- rbind(ped,pedFamily)
        ped$ID <- paste("ID", ped$ID, sep="")
        ped$Par1[ped$gener!=0] <-paste("ID", ped$Par1[ped$gener!=0], sep="")
        ped$Par2[ped$gener!=0] <-paste("ID", ped$Par2[ped$gener!=0], sep="")

        class(ped) <- c("pedigree","data.frame")
        return(ped)
}

simul.pedigree()
