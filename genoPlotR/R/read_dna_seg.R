################################################################################
# File reading functions: read dna_seg
################################################################################
# SUPPORT FUNCTIONS FOR READ_DNA_SEG_FROM_GENBANK
# Get start and end position for a currentFeature Object
# Please note, the function makes sure that only start and end of the
# following two types are ok:
# 'XXX (numeric)..(numeric)' or 'XXX complement((numeric)..(numeric))',
# where XXX may contain a-Z, 0-9 and
# also ! " # $ % & ' ( ) * + , - . / : ; < = > ? @ [ \ ] ^ _ ` { | } ~.
get_start <- function(line){
  ifelse (length(grep("complement", line)) > 0,
          start <- as.numeric(gsub("_|[[:blank:]]|[[:alpha:]]|\\(|\\)|\\.\\..*", "",
                   grep("^[[:graph:]]+ complement\\([[:digit:]]+\\.\\.[[:digit:]]+\\)$",
                   line, value=TRUE))),
          start <- as.numeric(gsub("_|[[:blank:]]|[[:alpha:]]|\\(|\\)|\\.\\..*", "",
                   grep("^[[:graph:]]+ [[:digit:]]+\\.\\.[[:digit:]]+$", line, value=TRUE))))
  start
}
get_end <- function(line){
  ifelse(length(grep("complement", line)) > 0,
         end <-
         as.numeric(gsub("_|[[:blank:]]|[[:alpha:]]|\\(|\\)|.*\\.\\.", "",
                         grep("^[[:graph:]]+ complement\\([[:digit:]]+\\.\\.[[:digit:]]+\\)$",
                              line, value=TRUE))),
         
         end <- as.numeric(gsub("_|[[:blank:]]|[[:alpha:]]|\\(|\\)|.*\\.\\.", "",
                                grep("^[[:graph:]]+ [[:digit:]]+\\.\\.[[:digit:]]+$",
                                     line, value=TRUE))))
  end
}
# Extracts data from feature lines.
extract_data <- function(extract, cF){
  extract <- gsub(extract, "", grep(paste("^",extract,sep=""), cF, value=TRUE))
  if (length(extract) == 0) { extract <- "NA" }
  extract[1]
}
# Added for support, to run an embl file directly
read_dna_seg_from_embl <- function(file, tagsToParse=c("CDS"), ...){
  read_dna_seg_from_file(file, tagsToParse, fileType="embl", ...)
}
# Added for support, to run an genbank file directly
read_dna_seg_from_genbank <- function(file, tagsToParse=c("CDS"), ...){
  read_dna_seg_from_file(file, tagsToParse, fileType="genbank", ...)
}

# MAIN FUNCTION
# Read genes from a GenBank file
read_dna_seg_from_file <- function(file, tagsToParse=c("CDS"),
                                   fileType="detect", meta_lines=2,
                                   gene_type="auto", header=TRUE,
                                   extra_fields=NULL, ...){
  
  # Import data from file into variable
  importedData <- readLines(file)
  
  # Find file type
  TYPE <- "Unknown"
  if (fileType == "detect" || fileType == "Detect" || fileType == "DETECT") {
    if (length(grep(">", importedData[1]))) {TYPE <- "Fasta"}
    if (length(grep("^ID", importedData))) {TYPE <- "EMBL"}
    if (length(grep("^LOCUS", importedData))) {TYPE <- "Genbank"}
  }
  else if (fileType == "EMBL" || fileType == "embl" || fileType == "Embl") {
    TYPE <- "EMBL"
  }
  else if (fileType == "Genbank" || fileType == "GENBANK" || fileType == "genbank") {
    TYPE <- "Genbank"
  }
  else if (fileType == "Ptt" || fileType == "PTT" || fileType == "ptt") {
    TYPE <- "PTT"
  }
  else if (fileType == "Fasta" || fileType == "FASTA" || fileType == "fasta") {
    TYPE <- "Fasta"
  }
  if (TYPE == "Unknown") {
    stop("fileType has to be either 'detect', 'embl', 'genbank' or 'ptt'. Note if file type is .ptt, please specify this rather than using 'detect'.")
  }
  
  # If type is PTT, call already made function
  if(TYPE == "PTT") {
    dna_seg <- read_dna_seg_from_ptt(file, meta_lines, header, ...)
    return(dna_seg)
  }
  else if (TYPE == "Fasta"){
    dna_seg <- read_dna_seg_from_fasta(file)
    return(dna_seg)
  }
  
  # If type isn't PTT do everything else...
  else {
    
    # Extract and name main segments
    if(TYPE == "Genbank") {
      mainSegments <- grep("^[[:alnum:]]", importedData)
      names(mainSegments) <- gsub("*| .*", "", grep("^[[:alnum:]]",
                                                    importedData, value=TRUE))
    }
    
    # SIMPLE ERROR HANDLING
    if(TYPE == "Genbank") {
      if(length(grep("FEATURES|DEFINITION", names(mainSegments))) < 2) {
        stop("FEATURES or DEFINITION segment missing in GBK File.")
      }
      if(length(grep("LOCUS", names(mainSegments))) != 1) {
        stop("Number of LOCUS should be 1.")
      }
    }
    if(TYPE == "EMBL") {
      if(length(grep("^ID|^FH", importedData)) < 2){
        stop("ID or FH segment missing in EMBL File.")
      }
      if(length(grep("^ID", importedData)) != 1) {
        stop("Number of ID lines should be 1.")
      }
    }
    
    # Extract META data
    if(TYPE == "EMBL") {
      seg_name <- gsub("^DE[[:blank:]]+", "",
                       grep("^DE", importedData, value=T))
    }
    if(TYPE == "Genbank") {
      seg_name <- gsub("DEFINITION {1,}", "",
                       importedData[mainSegments["DEFINITION"]])
    }
    
    # Extract features only, handles whether FEATURES is the last (or not)
    # entry in GBK file
    if(TYPE == "Genbank") {
      ifelse(which(names(mainSegments) == "FEATURES") == length(mainSegments),
             dataFeatures <-
             importedData[mainSegments["FEATURES"]:
                          (length(importedData) - 1)],
             dataFeatures <-
             importedData[mainSegments["FEATURES"]:
                          (mainSegments[which(names(mainSegments) ==
                                              "FEATURES")+1] - 1)])
    }
    if(TYPE == "EMBL") {
      dataFeatures <- grep("^FT", importedData, value=T)
    }
    
    # SIMPLE ERROR HANDLING
    if(TYPE == "Genbank") {
      if(length(dataFeatures) < 2){ stop("No FEATURES in GBK file.") }
    }
    if(TYPE == "EMBL") {
      if(length(dataFeatures) < 1){ stop("No FEATURES in GBK file.") }
    }
    
    
    # Extract each start line for each feature
    if(TYPE == "Genbank") {
      startLineOfFeature <- c(1:length(dataFeatures))[- grep("^ {6,}",
                                                             dataFeatures)]
    }  
    if(TYPE == "EMBL") {
      startLineOfFeature <- grep("FT   [[:alnum:]]", dataFeatures)
    }
    startLineOfFeature <- c(startLineOfFeature, length(dataFeatures)+1)
    
    
    # Define variables for storage
    nF <- length(startLineOfFeature)-1
    name <- character()
    start <- numeric()
    end <- numeric()
    strand <- numeric()
    length <- numeric()
    pid <- character()
    gene <- character()
    synonym <- character()
    product <- character()
    color <- character()
    proteinid <- character()
    feature <- character()
    geneType <- character()
    extra <- list()
    for (tag in extra_fields){
      extra[[tag]] <- character()
    }
    
    # Loop over all features                     
    for(counter in 1:nF){
       
      # Get feature, normally 20ish lines.
      currentFeature <- (dataFeatures[ startLineOfFeature[counter] :
                                      (startLineOfFeature[counter+1]-1) ])
      
      # Clean up feature, decreases number of lines etc.
      if(TYPE == "Genbank") {
        currentFeature <-
          gsub("^ |:|\"| $", "",
               gsub("[[:blank:]]+|[[:space:]]+",
                    " ", strsplit(paste(currentFeature,
                                        collapse=""), "   /")[[1]]))
      }
      if(TYPE == "EMBL") {
        currentFeature <-
          gsub("^ |:|\"| $", "",
               gsub("[[:blank:]]+|[[:space:]]+", " ",
                    strsplit(paste(gsub("FT", "", currentFeature),
                                   collapse=""), "   /")[[1]]))
      }
      
      # If feature is of a type to parse. Default is only CDS tags.
      if(length(grep(gsub(" [[:print:]]+", "", currentFeature[1]),
                     tagsToParse)) > 0){
        # Create list with exons to parse
        tag <- gsub(" [[:graph:]]+", "", currentFeature[1])
        qualif <- gsub("[[:graph:]]+ ", "", currentFeature[1])
        exonVector <- strsplit(gsub("[[:alpha:]]|_| |\\(|\\)|", "",
                                   currentFeature[1]), ",")
        if (length(grep("complement", currentFeature[1])) > 0){
          exonVector <- paste(tag, " complement(", exonVector[[1]], ")",
                              sep="")
        }
        if (length(grep("complement", currentFeature[1])) == 0){
          exonVector <- paste(tag, " ", exonVector[[1]], sep="")
        }

        # If there are more than one exon, insert introns
        if (length(exonVector) > 1) {
          exonVector2 <- 0
          for (i in 1:(length(exonVector)-1)) {
            if (length(grep("complement", currentFeature[1])) > 0) {
              exonVector2 <-
                c(exonVector2,exonVector[i],
                  paste(tag, "_intron complement(", get_end(exonVector[i])+1,
                        "..", get_start(exonVector[i+1])-1, ")", sep=""))
            }
            if (length(grep("complement", currentFeature[1])) == 0) {
              exonVector2 <- c(exonVector2,exonVector[i],
                               paste(tag, "_intron ", get_end(exonVector[i])+1,
                                     "..", get_start(exonVector[i+1])-1,
                                     sep=""))
            }
          }
          exonVector2[length(exonVector)*2] <- exonVector[length(exonVector)]
          exonVector2 <- exonVector2[2:length(exonVector2)]
          exonVector <- exonVector2
        }
        
        # For each exon in currentFeature 
        for (currentExon in exonVector) {
          
          # Set currentExon to currentFeature line 1
          currentFeature[1] <- currentExon
        
          # SIMPLE ERROR HANDLING AND Only continue parsing if start and stop
          # is valid...
          # Extract gene name or ID AND start and end, THEN, check if it's ok.
          ifelse(length(grep("gene=", currentFeature)) > 0,
                 nameTEMP <- extract_data("gene=", currentFeature),
                 nameTEMP <- extract_data("locus_tag=", currentFeature))
          startTEMP <-get_start(currentFeature[1])
          endTEMP <- get_end(currentFeature[1])
          if(length(startTEMP) == 0 || length(endTEMP) == 0) {
            warning(paste("Start and stop position invalid for "),
                    nameTEMP, ". This entry has been excluded.", sep="") }
          
          # Continue if start and end is ok... Otherwise, skip this feature
          if (length(startTEMP) > 0 && length(endTEMP) > 0) {
            
            # Save name, start, end, length
            name <- c(name, nameTEMP)
            start <- c(start, startTEMP)
            end <- c(end, endTEMP)
            length <- c(length, (get_end(currentFeature[1]) -
                                 get_start(currentFeature[1]) + 1)/3 - 1)
            
            # Set strand to 1 or -1
            ifelse (length(grep("complement", currentFeature[1])) > 0,
                    strand <- c(strand, -1), strand <- c(strand, 1))
            
            # Extract PID
            if (TYPE == "Genbank") {
              pid <- c(pid, extract_data("db_xref=GI", currentFeature))
            }
            if (TYPE == "EMBL") {
              pidTEMP <- extract_data("db_xref=UniProtKB/Swiss-Prot",
                                      currentFeature)
              if (pidTEMP == "NA"){
                pidTEMP <- extract_data("db_xref=UniProtKB/TrEMBL",
                                        currentFeature)
              }
              pid <- c(pid, pidTEMP)
            }
            
            # Extract gene
            ifelse(length(grep("gene=", currentFeature)) > 0,
                   gene <- c(gene, extract_data("gene=", currentFeature)),
                   gene <- c(gene, "-"))
            
            # Extract synonym
            synonym <- c(synonym, extract_data("locus_tag=", currentFeature))
            
            # Extract protein ID
            proteinid <- c(proteinid, extract_data("protein_id=",
                                                   currentFeature))
            
            # Extract product
            product <- c(product, extract_data("product=", currentFeature))

            # Extract color
            color <- c(color, extract_data("(color|colour)=", currentFeature))

            # Extract extra
            for (tag in names(extra)){
              extra[[tag]] <- c(extra[[tag]],
                                extract_data(paste(tag, "=", sep=""),
                                             currentFeature))
            }
            
            # Set geneType
            if (length(grep("intron", currentFeature[1])) > 0){
              geneType <- c(geneType, "introns")
            }
            if (length(grep("intron", currentFeature[1])) == 0){
              geneType <- c(geneType, gene_type)
            }
            
            # Return tag feature info, with or without added _pseudo tag
            if (length(grep("^pseudo", currentFeature)) > 0) {
              feature <- c(feature,
                           paste(gsub(" [[:print:]]+", "",
                                      currentFeature[1]), "_pseudo", sep=""))
            }
            if (length(grep("^pseudo", currentFeature)) == 0) {
              feature <- c(feature, gsub(" [[:print:]]+", "",
                                         currentFeature[1]))
            }
            
            
          # end of parse if start and end ok
          }
          
        # End of exon loop
        }
        
      # End of CDS loop
      }
      
    # End of loop over all features
    }
    
    # SIMPLE ERROR HANDLING
    if(is.numeric(start) == FALSE) stop("Start is not numeric.")
    if(is.numeric(end) == FALSE) stop("End is not numeric.")
    if(is.numeric(length) == FALSE) stop("Length is not numeric.")
    if(is.numeric(strand) == FALSE) stop("Strand is not numeric.")
    if(is.character(pid) == FALSE) stop("PID is not character.")
    if(is.character(name) == FALSE) stop("Name is not character.")
    if(is.character(gene) == FALSE) stop("Gene is not character.")
    if(is.character(synonym) == FALSE) {
      stop("Synonym is not character.")
    }
    if(is.character(product) == FALSE) {
      stop("Product is not character.")
    }
    if(is.character(proteinid) == FALSE) {
      stop("Protein ID is not character.")
    }
    if(is.character(feature) == FALSE) {
      stop("Feature is not character.")
    }
    # Check color
    # Eventually, change Artemis colors to their RGB equivalent
    artCol <- artemisColors()
    if (length(color) > 0 && all(color %in% c("NA", artCol$n))) {
      for (i in 1:length(color)){
        if (color[i] != "NA") color[i] <- artCol$colors[artCol$n == color[i]]
      }
    }

    # If gene_type is auto, change form arrows to blocks, otherwise arrows
    if (length(grep("intron", geneType))>=1 && gene_type == "auto") {
      geneType[geneType== "auto"] <- "exons"
    }
    if (length(grep("intron", geneType))==0 && gene_type == "auto") {
      geneType[geneType== "auto"] <- "bars"
    }
    
    # Cut table to include only added features    
    table <- data.frame(name=name, start=start, end=end, strand=strand,
                        length=length, pid=pid, gene=gene, synonym=synonym,
                        product=product, proteinid=proteinid, feature=feature,
                        gene_type=geneType, stringsAsFactors=FALSE)
    if (!all(color == "NA")){
      color[color == "NA"] <- "blue"
      table$col <- color
    }
    ## Adding extra fields
    for (tag in extra_fields){
      table[[tag]] <- extra[[tag]]
    }
    # SIMPLE ERROR HANDLING
    if (dim(table)[1] == 0)
      return (NULL)
      #stop("Nothing to return in table. I.e. no features extracted.")
    
    # Go to next function
    .read_dna_seg(table, seg_name, ...)
    
  # End of not PTT
  }
  
# End of genbank to dna_seg
}
read_dna_seg_from_fasta <- function(file, ...){
  # read data
  data <- readLines(file)
  title <- data[1]
  if (length(grep("^>\\w+", title, perl=TRUE)) < 1){
    stop(paste(file, "does not seem like a valid fasta file"))
  }
  name <- unlist(strsplit(title, " ", fixed=TRUE))[1]
  name <- substr(name, 2, nchar(name))
  seq <- data[-1]
  len <- sum(nchar(seq))
  dna_seg <- as.dna_seg(data.frame(name=name, start=1, end=len, strand=1,
                                   stringsAsFactors=FALSE), ...)
  return(dna_seg)
}
# reading genes from a file. Use source=tab or ptt to specify type
read_dna_seg_from_ptt <- function(file, meta_lines=2, header=TRUE, ...){
  # reads meta info
  seg_name <- readLines(file, n=1)
  seg_name <- strsplit(seg_name, "/,|-/", fixed=TRUE)[[1]][1]
  # reads ptt table
  ptt <- read.table(file, skip=meta_lines, as.is=TRUE, header=header,
                    sep="\t", quote="")
  if (header){
    names(ptt) <- tolower(names(ptt))
  }
  else {
    names(ptt) <- c("location", "strand", "length", "pid", "gene",
                    "synonym", "code", "cog", "product")
  }
  # parse location
  location <- strsplit(ptt$location, "..", fixed=TRUE)
  start <- as.numeric(sapply(location, function(x) x[[1]]))
  end <- as.numeric(sapply(location, function(x) x[[2]]))
  # parse strand
  strand <- ptt$strand
  strand[strand=="-"] <- -1
  strand[strand=="+"] <- 1
  strand <- as.numeric(strand)
  # parse gene name from name or synonym if not present
  name <- ifelse(ptt$gene == "-", ptt$synonym, ptt$gene)
  table <- data.frame(name=name, start=start, end=end, strand=strand,
                      length=ptt$length, pid=ptt$pid, gene=ptt$gene,
                      synonym=ptt$synonym, code=ptt$code, cog=ptt$cog,
                      product=ptt$product,
                      stringsAsFactors=FALSE)
  .read_dna_seg(table, seg_name, ...)
}
read_dna_seg_from_tab <- function(file, header=TRUE, ...) {
  table <- read.table(file, as.is=TRUE, header=header, sep="\t", quote="")
  if (ncol(table) < 4) stop("Insufficent number of columns in table")
  col_names <-  c("name", "start", "end", "strand")
  names(table)[1:length(col_names)] <- col_names
  # parse name from file name by default
  seg_name <- basename(file)
  .read_dna_seg(table, seg_name, ...)
}
.read_dna_seg <- function(table, seg_name, reverse=FALSE, xlim=NULL, ...){
  # check args
  if (ncol(table) < 4) stop("Insufficent number of columns in table")
  if (nrow(table) < 1) stop("No lines in table")
  col_names <-  c("name", "start", "end", "strand")
  if (!all(col_names %in% names(table)))
    stop("Table should contain at least columns name, start, end and strand")
  # make dna_seg object, set seg_name attribute
  dna_seg <- as.dna_seg(table, ...)
  dna_seg <- trim.dna_seg(dna_seg, xlim)
  if (reverse) dna_seg <- reverse.dna_seg(dna_seg)
  attr(dna_seg, "seg_name") <- seg_name
  dna_seg
}
