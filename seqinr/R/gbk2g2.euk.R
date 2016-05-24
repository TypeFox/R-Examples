gbk2g2.euk <- function(
  gbkfile = system.file("sequences/ame1.gbk", package = "seqinr"), g2.coord = "g2.coord")
{
  input <- readLines(gbkfile)
  outfile <- file( description = g2.coord, open ="w+")

 
  #
  # Keep lines with annotation flag:
  #

  cds <- which(substring(input,1,8) == "     CDS")

  features <- which(substr(input,6,6)!=" ")
  
  features <- c(features,length(input))

  genes <- which(substr(input,22,26)=="/gene")
  
  genes.cds <- character(length(cds))
  
  for(i in seq_len(length(cds))){
    print(i)
    this.gene=genes[which(genes>cds[i])[1]]
    
    nextfeat=features[which(features>cds[i])[1]]
    
    if(this.gene<nextfeat){
      genes.cds[i]=unlist(strsplit(input[this.gene],split="\""))[2]
      
    }
    else{
      print(paste("no id for cds #",i))
    }

  }

  
  #
  # Extract boundaries strings
  #
  get.boundaries <- function( index.line )
  {
    
    join <- grep("join",input[index.line])
   
    if(length(join)>0){ ## there are introns !!

     
      end.par <- grep("\\)",input[index.line])
     

      if(length(end.par)==0){ ## the exons are written on more than one line

        next.end.par <-grep("\\)",input)
        next.end.par <-next.end.par[which(next.end.par>index.line)[1]] ## we take the first ending parenthesis
        exons <- paste(input[index.line],paste(substr(input[(index.line+1):next.end.par],22,nchar(input[(index.line+1):next.end.par])),sep="",collapse=""),collapse="",sep="")
      }
      else{
        exons <- input[index.line]
      }
     

    

    complement <- grep("complement",exons)

    if(length(complement)==0){

      exons=unlist(strsplit(exons, split="\\("))[2]
      exons=unlist(strsplit(exons, split="\\)"))[1]

      exons=unlist(strsplit(exons,split=","))
      exons.begin=unlist(lapply(exons, function(x) unlist(strsplit(x, split="\\.\\."))[1]))
      exons.end=unlist(lapply(exons, function(x) unlist(strsplit(x, split="\\.\\."))[2]))
      

      
    }

    else{
      exons=unlist(strsplit(exons, split="\\("))[3]
      exons=unlist(strsplit(exons, split="\\)"))[1]

      exons=unlist(strsplit(exons,split=","))
      exons.end=unlist(lapply(exons, function(x) unlist(strsplit(x, split="\\.\\."))[1]))
      exons.begin=unlist(lapply(exons, function(x) unlist(strsplit(x, split="\\.\\."))[2]))
      
      
      
    }
    }
    else{
      complement <- grep("complement",input[index.line])

      exons=unlist(strsplit(input[index.line],split=" "))
      exons=exons[exons!=""]
      exons=exons[2]

      if(length(complement)==0){
        exons.begin=unlist(strsplit(exons, split="\\.\\."))[1]
        exons.end=unlist(strsplit(exons, split="\\.\\."))[2]
        
      }
      else{
        exons=unlist(strsplit(exons, split="\\("))[2]
        exons=unlist(strsplit(exons, split="\\)"))[1]
        exons.begin=unlist(strsplit(exons, split="\\.\\."))[2]
        exons.end=unlist(strsplit(exons, split="\\.\\."))[1]
      }
      
    }
    
  
    return(list(exons.begin,exons.end))

    
  }

   

  
  
  
  already=character(0)
   
  if(length(cds)>0){
    for(i in seq_len(length(cds))){
      
      boundaries=get.boundaries(cds[i])
      exons.begin=boundaries[[1]]
      exons.end=boundaries[[2]]
      
      
      
      for(j in seq_len(length(exons.begin))){
        phrase=paste(genes.cds[i],exons.begin[j],exons.end[j])
        
        if(!phrase%in%already){
          writeLines(paste(genes.cds[i],exons.begin[j],exons.end[j],sep=" "),outfile)
          already=c(already,phrase)
        }
      }
    }
     
   } 
    

 
  close(outfile)
  
}

