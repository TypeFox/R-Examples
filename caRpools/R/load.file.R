load.file = function(filename = NULL, header= TRUE, sep="\t", comment.char="", type = NULL)
  
  #miaccs = "MIACCS.xlsx", files = NULL, group=NULL, libraryfasta=NULL, gff=NULL, path.to.scripts="./scripts", path.to.data="./data", extract=TRUE, seq.pattern = "default", machine.pattern = "default", sensitivity = "very-sensitive-local", match = "perfect", reversecomplement = FALSE, threads = 1, bowtieparams = "", header= TRUE, sep="\t", comment.char="", type = NULL, agg.function=sum, extractpattern=expression("^(.+?)_.+"))
{
#   # create caRpools Environment
#   if(!exists("cp", mode="environment"))
#   {
#     cp <- new.env()
#   }
#  
#   if(is.null(libraryfasta))
#   {
#     stop("No library FASTA file provided. Please provide the filename inclduing the .fasta extension of your library.")
#   }
  
  # make loading files more convenient by just telling the function the file names as either
  # a comma separated string
  # or a list
  # that must indicate wether it is
  # belonging to a group (numeric, integer), in this case it is used as a replicate
  # it is a read-count file or a fastq which needs extraction/mapping or a LIBRARY fasta file
  # the name of the file as a description (either comma separated or as a list in the SAME order)
  
# arguments:
  # miaccs.file (MUST BE SET), default=MIACCS.xls
  # files (list or comma-separated)
  # group (by default untreated and treated, which is either untreated(1) or treated(2) by default)
    # group = list(untreated = c("1","2"), treated = c("3","4") )
  # libraryfasta (required for some tools AND data mapping, must be FASTA format)
  # gff E-CRISP GFF file (optional)
  # scriptpath=NULL, datapath=NULL
  # extract = FALSE, pattern = "default"
  # maschinepattern = "default", createindex = FALSE
  # bowtie2file = NULL, mapping = FALSE
  # reversecomplement = FALSE, threads = 1
  # bowtieparams = "", sensitivity = "very-sensitive-local"
  # match = "perfect"
  
  
###### READ MIACCS FILE
  # read XLSX file
#   if(type=="MIACCS")
#   {
#     filemiaccs = xlsx::read.xlsx (miaccs.file, sheetName="MIACCS", header=FALSE, stringsAsFactors=FALSE)
#     
#     # only take those with information
#     filemiaccs = filemiaccs[which(filemiaccs[,1] !=""),]
#     
#     #set identifiers to as rownames
#     rownames(filemiaccs) = filemiaccs[,1]
#     
#     cp$miaccs = filemiaccs
#   }
#   else
# {
# 
# ##### LOAD FILES SETUP
#   
# # check if comma separated or list
#   # store number of files for each treatment group which its name in separate data frame
# if(!is.list(files))
# {
#   # no list, so comma-separated
#   files <- as.list(unlist(strsplit(files, ",")))
# }
#   
# ############ Load files
# ## Eiether with extraction / mapping or direct read-count
# # it is a list, so we go through the list an construct a data frame
# # sgRNA | untreated | untreated | untreated | treated |
# # Read Count or FastQ files which need to be mapped/extracted?
# # First check ending, if Fastq go for extraction and mapping or mapping only, if txt go for direct loading.
# 
#   
#   ######## LOAD FASTA LIBRARY
#   # First load FASTA LIBRARY to create data frame needed
#   # Read fasta reference file for sgRNA sequence list at the end
#   # fasta file has      <IDENTIFIER
#   # next line           SEQUENCE
#   if("seqinr" %in% rownames(installed.packages()) == FALSE) {install.packages("seqinr")}
#   #library("seqinr")
#   file.raw <- seqinr::read.fasta(file = paste(datapath, libraryfasta, sep="/"), 
#                                  seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
# 
#   # make as data frame
#   cp$libFILE = data.frame(
#     design = attributes(file.raw),
#     sequence = unlist(file.raw),
#     stringsAsFactors=FALSE)
#   
#   cp$readcount = data.frame(
#     design = attributes(file.raw),
#     stringsAsFactors=FALSE)
#   colnames(cp$readcount) = c("design")
#  
# ##### Load Files with Data 
# # go through files with lapply
#   files.loaded = lapply(files, function(x)
#   {
#     x <- as.character(x)
#     # check if txt or fastq
#     if(length(unlist(strsplit(x,"[.]"))) < 2)
#     {
#       stop("No file extension provided. Please provide .txt extension for read-count files and .fastq for raw data FASTQ files.")
#     }
#     else if(length(unlist(strsplit(x,"[.]"))) > 2)
#     {
#       stop("The filename contains ., which is not allowed.")
#     }
#     else
#     {
#       # check what kind of files are provided
#       # check if extraction/mapping is necessary
#       # do the mapping etc and load the final read-count file
#       
#       if(unlist(strsplit(x,"[.]"))[2] == "fastq" || unlist(strsplit(x,"[.]"))[2] == "FASTQ")
#       {
#         # FASTQ file will be extracted if necessary and mapped against the fasta library reference
#         fileFASTQ = data.extract(
#           scriptpath=path.to.scripts, datapath=path.to.data,fastqfile=unlist(strsplit(x,"[.]"))[1],
#           extract=extract, pattern=seq.pattern, machinepattern=machine.pattern, createindex=TRUE,
#           referencefile=unlist(strsplit(libraryfasta,"[.]"))[1],
#           mapping=TRUE, reversecomplement=reversecomplement, threads, bowtieparams,
#           sensitivity=sensitivity,match=match) 
#         
#         # Returned value fileFASTQ will provide the file name including the extension
#         fileread = read.table(paste(path.to.data, fileFASTQ, sep="/"), header=header, sep=sep, comment.char=comment.char, stringsAsFactors=FALSE)
#         colnames(fileread) = c("design", "reads")
#       }
#       else if(unlist(strsplit(x,"[.]"))[2] == "txt" || unlist(strsplit(x,"[.]"))[2] == "TXT")
#       {
#         # or directly load the read-count file
#         fileread = read.table(paste(path.to.data, x,sep="/"),header=header, sep=sep, comment.char=comment.char, stringsAsFactors=FALSE)
#         #return(fileread.test)
#         colnames(fileread) = c("design", "reads")
#       }
#       else
#       {
#         # file extension neither fstq nor txt
#         stop("The file extension is neither .txt nor .fastq . Please provide files with appropriate extension.")
#       }
#  
#       #str(fileread)
#       # Merge to data frame
#       cp$readcount <- merge(cp$readcount, fileread, by="design")
#       
#     }
#   }
#   )
#   
# ### check for groups and apply to data frame
# cp$treatmentgroup = group
# 
# colnames(cp$readcount) <- c("design", cp$treatmentgroup$untreated, cp$treatmentgroup$treated)
# rownames(cp$readcount) <- cp$readcount$design
# 
# ###### BACKUP old carpools -> create single data frames for each replicate
# 
# # iterate through treatment group definition
# 
# for (name in names(cp$treatmentgroup)) {
#   print(name)
#   cp[name] = list()
#   str(cp)
#   for(i in 1:length(cp$treatmentgroup[[name]]))
#   {
#     print(as.numeric(cp$treatmentgroup[[name]][i])+1)
#     cp$files = c(cp$files, data.frame(design = cp$readcount$design,
#                               count = cp$readcount[,as.numeric(cp$treatmentgroup[[name]][i])+1],
#                               stringsAsFactors = FALSE))
#   }
#   #str(cp$files)
#   #print(length(cp$treatmentgroup[[name]]))
#   #print(cp$treatmentgroup[[name]])
# }
# paste(name,cp$treatmentgroup[[name]][i], sep="") = 
# # iterate through length
# #for(i in 1: ncol(cp$readcount))
# #{
# #  
# #}
# cp$files = list()

###### Load E-CRISP / CLD GFF file
# cp$gff


###### Aggregate to gene level and store in seperate variable
# cp$readcount.gene



###### Normalize Read Count data (normalize = TRUE, norm.function = median)
# cp$readcount.norm -> will be used in functions if normalize == TRUE



  
  
  if(identical(type,"fastalib"))
  {
    # Read fasta reference file for sgRNA sequence list at the end
    # fasta file has      <IDENTIFIER
    # next line           SEQUENCE
    if("seqinr" %in% rownames(installed.packages()) == FALSE) {install.packages("seqinr")}
    #library("seqinr")
    file.raw <- seqinr::read.fasta(file = filename, 
                           seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)

    
    #str(attributes(file.raw))
    # make as data frame
    file = data.frame(
      ID = attributes(file.raw),
      sequence = unlist(file.raw),
      stringsAsFactors=FALSE)
  }
  
# read XLSX file
else if(identical(type,"xlsx"))
{
  #library("xlsx")
  file = xlsx::read.xlsx (filename, sheetName="MIACCS", header=FALSE, stringsAsFactors=FALSE)
  
  # only take those with information
  file = file[which(file[,1] !=""),]
  
  #set identifiers to as rownames
  rownames(file) = file[,1]
}

  
  else
  {
   # OLD for single file loading
    file = read.table(filename,header=header, sep=sep, comment.char=comment.char)
    
    # Multiple file loading with list of file name
  }
  
  return(file)
}
#}