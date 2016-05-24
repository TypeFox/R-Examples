#### Function reverse as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

reverse <-
function(dna){
     splitstr <- substring(dna,1:nchar(dna),1:nchar(dna))
     result <- paste(splitstr[length(splitstr):1], collapse = "")
     return(result)
}

