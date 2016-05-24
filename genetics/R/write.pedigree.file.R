# $Id: write.pedigree.file.R 666 2006-03-10 18:19:34Z nj7w $

write.pedigree.file <- function(data,
                                family, pid, father, mother, sex,
                                file="pedigree.txt"
                                )
  {
    # pedigree file format
    # --------------------
    #
    # <family> <pid> <father> <mother> <sex> <genotype_1> ... <genotype_n>
    #
    # <family> is a unique identifier for each family, and within each family
    # <pid> is a unique identifier for an individual.
    # <father> and <mother> identify the individuals father and mother
    #    (if this line refers to a founder, these should be set to
    #     zero).
    #  <sex> denotes the individuals sex, using the convention
    #     1=male, 2=female.
    #
    # Each <genotype> is encoded as two integer allele numbers.

    if(missing(family)) 
      family <- 1:nrow(data)
    if(missing(pid))
      pid <- 1:nrow(data)
    if(missing(father))
      father <- rep(0,nrow(data))
    if(missing(mother))
      mother <- rep(0,nrow(data))
    if(missing(sex))
      sex <- rep(0,nrow(data))

    pedigree <- list()
    pedigree$family <- format(family)
    pedigree$pid    <- format(pid)
    pedigree$father <- format(father)
    pedigree$mother <- format(mother)
    pedigree$sex    <- format(sex)

    which <- sapply(data, is.genotype)
    if(!all(which)) warning("Data contianed non-genotype variables.",
                            " These have been ignored: ",
                            paste(colnames(data)[!which]) )

    data <- data[,which]

    allele.number <- function(g, ind) {
      as.numeric(factor(allele(g, ind),
                        levels = allele.names(g)))
    }

    for( col in names(data) )
      {
        name.1 <- paste(col,".1")
        name.2 <- paste(col,".2")

        

       ## allele.number <- function(g, ind)
       ##   as.numeric(as.factor(allele(g,ind), levels=allele.names(g, ind) ))
        
        pedigree[[name.1]] <- allele.number( data[[col]], 1)
        pedigree[[name.2]] <- allele.number( data[[col]], 2)
      }
    pedigree <- as.data.frame(pedigree)
    
    # NA's are represented as 0
    pedigree[is.na(pedigree)] <- 0

    write.table(pedigree,   file=file, sep=" ", row.names=FALSE,
                col.names=F, quote=F)
  }

write.marker.file<-function(data, location, file="marker.txt")
{
    # marker.map file format
    # --------------------
    #
    # MARKERID   NAME      LOCATION
    # <col#>     <name>    <location>

    which <- sapply(data, is.genotype)
    if(!all(which)) warning("Data contianed non-genotype variables.",
                            " These have been ignored: ",
                            paste(colnames(data)[!which]) )

    data <- data[,which]
    if(missing(location)) location <- 1:ncol(data)

    ## Create marker.map file data frame
    marker.map <- cbind(
                             formatC(1:ncol(data), width=8),
                             formatC(colnames(data), width=8, flag="-"),
                             formatC(location, width=8)
                             )

    marker.map <- rbind(c("MARKERID", "NAME    ","LOCATION"), marker.map)
    
    write.table(marker.map, file=file, sep="   ", row.names=FALSE,
                col.names=F, quote=F)    
}
