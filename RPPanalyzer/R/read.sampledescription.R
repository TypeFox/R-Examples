`read.sampledescription` <-
function (){
    ## check for sampledescription.txt in working directory
    allfiles <- dir()
    if ( ! "sampledescription.txt" %in% allfiles) {
         stop("can not find sampledescription.txt in current working directory.")
    }
    
    ## read in data frame
    sampleID<-read.delim("sampledescription.txt",header=T)
    
    ## fill empty cells with NAs
    #for (i in 1:ncol(sampleID)){
     #        empty.cells <- grep("^$",sampleID[,i])
     #        sampleID[empty.cells,i]= "not_available"
     #        }
    
                    
      ## check header for required columns
      reqCols <- c("plate","row","column","sample_type","sample","concentration")
      if ( !all( reqCols %in% colnames(sampleID))){
          stop("sampledescription file: columns are missing or header incorrect!")
      }
      
      ## check source plate description
      if ( mode(sampleID[,"plate"])!="numeric" | mode(sampleID[,"column"])!="numeric" | mode(sampleID[,"row"])!="numeric"){
          stop ("sampledescription file: data format in columns plate,column and row have to be numeric!") 
      }
      
      ## check sample types      
      reqtypes <- c("measurement","control","blank","neg_control")
      types <- as.character(levels(sampleID[,"sample_type"]))
      if ( !all(types %in% reqtypes)){
          stop("sampledescription file: allowed sample_types are: measurement, control, blank, neg_control")
      }
    
    ## generate sample identifier from source plate location  
    ID<-paste(sampleID[,"plate"],sampleID[,"row"],sampleID[,"column"],sep="")
    sampleID<-cbind(ID,sampleID)
    sampleID<-sampleID[,-c(2:4)]
    
    return(sampleID)
}

