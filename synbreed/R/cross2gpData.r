cross2gpData <- function(cross){
    # check for class
    # sometimes multiple classes in package 'qtl'
    if(all(class(cross)!="cross")) stop("object '",substitute(cross),"' not of class 'cross'")  
    # phenotypic data
    pheno <- cross$pheno
    # names of chromosomes
    chr <- names(cross$geno)
    # extract informations from first chromosome
    geno <- cross$geno[[chr[1]]]$data
    map <- data.frame(chr=chr[1],pos=cross$geno[[chr[1]]]$map)
    # loop over chromosomes 2:n
    if (length(chr)>1){
       for(i in chr[-1]){
          # add new genotypes
          geno <- cbind(geno,cross$geno[[i]]$data)
          # add new map info
          map <- rbind(map,data.frame(chr=i,pos=cross$geno[[i]]$map))
       }
    }
    # combine pheno, geno, and map  in class 'gpData'
    ret <- create.gpData(pheno=pheno,geno=geno,map=map)
    # recode
    # coding in cross: 1 = AA, 2 = AB, 3 = BB,
    ret$geno <- ret$geno-1
    return(ret)
} 