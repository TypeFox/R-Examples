SuperMatrix <- function(missing = "-",
                        prefix = "concatenated",
                        save = T){
  # get file names
  file.names <- list.files()
  # read DNA
  DNA <- list()
  for(i in 1:length(file.names)){
    print(paste("Reading alignment", i))
    DNA[[i]] <- read.dna(file=file.names[i], 
                         format = "f", 
                         as.character=T)
  }
  # calculate total alignment length
  total.length <- 0
  for(i in 1:length(DNA)){
    total.length <- total.length + ncol(DNA[[i]])
  }
  # get all taxa names
  taxa <- vector()
  for(i in 1:length(DNA)){
    counter <- length(taxa) + 1
    for(j in 1:nrow(DNA[[i]])){
      taxa <- c(taxa, 
                row.names(DNA[[i]])[!row.names(DNA[[i]]) %in% taxa])
    }
  }
  # create empty alignment
  seqmatrix <- matrix(as.character(missing),
                                    length(taxa),
                                    total.length)
  rownames(seqmatrix) <- taxa
  # create partition record holder
  partitions <- as.data.frame(matrix(,length(DNA),3))
  colnames(partitions) <- c("part", "start", "stop")
  #build the supermatrix
  c.col <- 0
  print("Creating supermatrix")
  for(i in 1:length(DNA)){
    print(paste("Processing alignment", i))
    gene <- DNA[[i]]
    print(i)
    for(j in 1:nrow(gene)){
      c.row <- which(rownames(seqmatrix) == rownames(gene)[j])
      seqmatrix[c.row, (c.col + 1):(c.col + ncol(gene))] <- gene[j, ]
    }
    partitions[i,1:3] <- c(file.names[i], c.col + 1, c.col + ncol(gene))
    c.col <- c.col + ncol(gene)
  }
  results <- list()
  results[[1]] <- partitions
  results[[2]] <- seqmatrix
  if(save == T){
    print("saving files")
    write.dna(seqmatrix, 
              file = paste(prefix, ".fasta", sep = ""), 
              format = "f")
    write.csv(partitions, 
              row.names = F, 
              file = paste(prefix, ".partitions.csv", sep = ""))
  }
  return(results)
}