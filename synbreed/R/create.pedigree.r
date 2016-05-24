create.pedigree <- function(ID,Par1,Par2,gener=NULL,sex=NULL,add.ancestors=FALSE,unknown=0){
     ID <- as.character(ID)
     Par1 <- as.character(Par1)
     Par2 <- as.character(Par2)
     n <- length(ID)
     # only unique IDs
     if(length(unique(ID))!=n)  warning("ID is not unique, removing duplicated individuals")

     # further checks
     if(length(Par1)!=n) stop("Par1 must have same length as ID")
     if(length(Par2)!=n) stop("Par2 must have same length as ID")

     # NA for unkonwn pedigree
     if(!(0 %in% ID | "0" %in% ID | any(unknown %in% ID))){
        Par1[Par1=="0" | Par1%in%unknown | Par1==0] <- NA
        Par2[Par2=="0" | Par2%in%unknown | Par2==0] <- NA
     } else if(as.character(ID[ID==0|ID=="0"]) == as.character(Par1[ID==0|ID=="0"]) & as.character(ID[ID==0|ID=="0"]) == as.character(Par2[ID==0|ID=="0"])){
        Par1<- Par1[ID!=0|ID!="0"]
        Par2<- Par2[ID!=0|ID!="0"]
        ID <- ID[ID!=0|ID!="0"]
        Par1[Par1=="0" | Par1%in%unknown | Par1==0] <- NA
        Par2[Par2=="0" | Par2%in%unknown | Par2==0] <- NA
     }
     Pars <- unique(c(Par1[!is.na(Par1)], Par2[!is.na(Par2)]))

     ancestors <- Pars[!Pars%in%ID]
     ID <- c(ancestors,ID)
     Par1 <- c(rep(NA,length(ancestors)),Par1)
     Par2 <- c(rep(NA,length(ancestors)),Par2)
     if (!is.null(gener)) gener <- c(rep(NA,length(ancestors)),gener)
     if (!is.null(sex)) sex <- c(rep(NA,length(ancestors)),sex)
     # update n
     n <- length(ID)

     # add gener if NULL
     if(is.null(gener)){
        #if(any(! Par1 %in% c(0,ID) ) | any(! Par1 %in% c(0,ID) )) stop("parents without pedigree, try to use argument 'add.ancestors=TRUE'")
        generOld <- gener <- rep(n+100,n)
        gener[is.na(Par1) & is.na(Par2)] <- 0
        i<-0
        while(!all(generOld == gener)){
            generOld <- gener
            gener[Par1 %in% ID[gener==i]] <- i+1
            gener[Par2 %in% ID[gener==i]] <- i+1
            i<-i+1
            if(i>n+10) break
        }
     }

     if(add.ancestors) ancestors <- FALSE

     # gener starts from 0
     if(!is.null(sex)) pedigree <- data.frame(ID=ID[!ID%in%ancestors],Par1=Par1[!ID%in%ancestors],Par2=Par2[!ID%in%ancestors],gener=gener[!ID%in%ancestors],sex=sex[!ID%in%ancestors],stringsAsFactors=FALSE)
     else pedigree <- data.frame(ID=ID[!ID%in%ancestors],Par1=Par1[!ID%in%ancestors],Par2=Par2[!ID%in%ancestors],gener=gener[!ID%in%ancestors],stringsAsFactors=FALSE)

     # removing duplicated entries
     pedigree <- pedigree[!duplicated(pedigree), ]

     # sort by generation
     pedigree <- pedigree[order(pedigree$gener,partial=pedigree$ID),]
     class(pedigree) <- c("pedigree","data.frame")

     # missing values as 0
     pedigree[is.na(pedigree)] <- 0

     return(pedigree)
}
