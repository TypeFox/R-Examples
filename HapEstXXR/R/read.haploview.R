read.haploview <-
function ( ped.file , map.file ) {

  inped.tab <- read.table ( ped.file , header=F , stringsAsFactors=F)

  famid <- inped.tab[,1]
  patid <- inped.tab[,2]
  dad <- inped.tab[,3]
  mom <- inped.tab[,4]
  sex <- inped.tab[,5]
  trait <- inped.tab[,6]
  geno <- allele2to1(as.matrix(inped.tab[,7:(dim( inped.tab)[2])]))
  storage.mode(geno) <- "integer"


  inmap.tab <-  read.table ( map.file , header=F , stringsAsFactors=F )
  marker.names <- inmap.tab [,1]
  marker.position  <- inmap.tab [,2]

  rownames(geno) <- patid
  colnames(geno) <- marker.names
  
  return (  list (famid=famid , patid=patid , dad=dad , mom=mom , sex=sex ,
       trait=trait , genotypes=geno , marker.names=marker.names ,
       marker.position=marker.position  ) )

}
