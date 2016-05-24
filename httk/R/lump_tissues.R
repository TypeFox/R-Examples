# This function returns the flows, volumes, and partition coefficients for the
# lumped tissues specified in tissue list
# Ktissue2plasma -- tissue to free plasma concentration partition coefficients
#                   for every tissue specified by Schmitt (2008) (the tissue.data table)
# tissuelist -- a list of character vectors, the name of each entry in the list
#               is a lumped tissue, the words in the vector are the Schmitt (2008)
#               tissues that are to be lumped, for example:
#               tissuelist<-list(Rapid=c("Brain","Kidney"))
# species specifies the flow.col and vol.col in the tissuedata.table
lump_tissues <- function(Ktissue2pu.in,
                            tissuelist=NULL,
                            species="Human")
{
  physiology.data <- physiology.data
  tissue.data <- tissue.data
  Parameter <- NULL
  if(length(Ktissue2pu.in) != length(tissue.data[,1]) | !all(tissue.data[,1] %in% c(substr(names(Ktissue2pu.in),2,nchar(names(Ktissue2pu.in))-3)
  [!substr(names(Ktissue2pu.in),2,nchar(names(Ktissue2pu.in))-3) %in% 'rbc'],'red blood cells'))) stop(paste('Ktissue2pu.in must contain the tissues from tissue.data:',paste(tissue.data[,1],collapse=', ')))
  if (!(species %in% colnames(physiology.data)))
  {
    if (toupper(species) %in% toupper(colnames(physiology.data)))
    {
      species <- colnames(physiology.data)[toupper(colnames(physiology.data))==toupper(species)]
    } else stop(paste("Tissue data for",species,"not found."))
  }
  vol.col <- paste(species,"Vol (L/kg)",sep=" ")  
  flow.col <- paste(species,"Flow (mL/min/kg^(3/4))",sep=" ")

# Initialize the output lists:
	vol <- list()
	flow <- list()
	Ktissue2pu.out <- list()
 
# The vector all.tissues indicates whether each tissue in tissue.data has been lumped yet (TRUE/FALSE)
	all.tissues <- rep(FALSE,length(tissue.data$Tissue))
	names(all.tissues) <- tissue.data$Tissue
  #Renames pcs to match tissue names
  names(Ktissue2pu.in) <- substr(names(Ktissue2pu.in),2,nchar(names(Ktissue2pu.in))-3)
  names(Ktissue2pu.in)[names(Ktissue2pu.in) == 'rbc'] <- 'red blood cells'
# Blood cells only need a partioncoefficient:
  Ktissue2pu.out[["red blood cells"]] <- Ktissue2pu.in[["red blood cells"]]	
  all.tissues["red blood cells"] <- T
 
# This loop adds up the volumes and flows for the tissues within each lumped tissue
# as well as Red blood cells
	for (this.lumped.tissue in c(names(tissuelist),"cleanup"))
	{
# Anything that has not yet been lumped is added to the lumped tissue "Rest"
		if (this.lumped.tissue == "cleanup")
		{
			this.lumped.tissue <- "rest"
# First check to see if rest has been created and create it if it is missing:
			if (!("rest" %in% names(vol)))
			{
				vol[["rest"]] <- 0
				flow[["rest"]] <- 0
        Ktissue2pu.out[["rest"]] <- 0
			}
# Every tissue not already lumped gets added to "Rest"
			these.lumped.tissues <- tissue.data$Tissue[!all.tissues]
		}	else{
			vol[[this.lumped.tissue]] <- 0
			flow[[this.lumped.tissue]] <- 0
			Ktissue2pu.out[[this.lumped.tissue]] <- 0
			these.lumped.tissues <- tissuelist[[this.lumped.tissue]]
		}
# Loop over every tissue that is lumped into the tissue:   
		for (this.tissue in these.lumped.tissues)
		{
			if (!(this.tissue %in% tissue.data$Tissue))
				stop(paste(this.tissue,"not in list:",paste(tissue.data$Tissue,collapse=', ')))
			if (all.tissues[[this.tissue]] & this.tissue !="rest")
				stop(paste(this.tissue,"assigned to multiple lumped tissues"))

# Mark that this tissue has been lumped:
			all.tissues[[this.tissue]] <- TRUE
# Find the row in the tissue.data table that corresponds to this tissue: 
			this.row <- tissue.data$Tissue==this.tissue
			
#Add the volume for this tissue to the lumped tissue:
			vol[[this.lumped.tissue]] <- vol[[this.lumped.tissue]] + as.numeric(tissue.data[this.row,vol.col])
#Add the flow for this tissue to the lumped tissue:                             
			flow[[this.lumped.tissue]] <- flow[[this.lumped.tissue]] + as.numeric(tissue.data[this.row,flow.col])
			 
#Add a contribution to the partition coefficient weighted by the volume of this tissue:

			Ktissue2pu.out[[this.lumped.tissue]] <- Ktissue2pu.out[[this.lumped.tissue]] + as.numeric(tissue.data[this.row,vol.col])*Ktissue2pu.in[[this.tissue]]
		}
#Calculate the average parition coefficient by dividing by the total volume of
#the lumped tissue:
		Ktissue2pu.out[[this.lumped.tissue]] <- Ktissue2pu.out[[this.lumped.tissue]]/vol[[this.lumped.tissue]]
	}

  # Must have tissue-specific flows for these tissues (even if lumped) in order
  # to calculate other quantities (e.g. rate of metabolism, renal clearance):
  for (this.tissue in c("liver","gut","kidney"))
    if (is.null(flow[[this.tissue]])) 
    {
      if(this.tissue %in% tissue.data[,"Tissue"]) flow[[this.tissue]] <- as.numeric(tissue.data[tissue.data[,"Tissue"]==this.tissue,flow.col])
      else if (paste(this.tissue,"s",sep="") %in% tissue.data[,"Tissue"]) flow[[this.tissue]] <- as.numeric(tissue.data[tissue.data[,"Tissue"]==paste(this.tissue,"s",sep=""),flow.col])
      else stop(paste("Tissue",this.tissue,"not found in tissue.data table."))            
    }

  # Must have tissue-specific volumes for these tissues (even if lumped) in order
  # to calculate other quantities (e.g. rate of metabolism):
    for (this.tissue in c("liver"))
    if (is.null(vol[[this.tissue]])) 
    {
      if (this.tissue %in% tissue.data[,"Tissue"]) vol[[this.tissue]] <- as.numeric(tissue.data[tissue.data[,"Tissue"]==this.tissue,vol.col])
      else if (paste(this.tissue,"s",sep="") %in% tissue.data[,"Tissue"]) vol[[this.tissue]] <- as.numeric(tissue.data[tissue.data[,"Tissue"]==paste(this.tissue,"s",sep=""),vol.col])
      else stop(paste("Tissue",this.tissue,"not found in tissue.data table."))
    }
    names(Ktissue2pu.out)[names(Ktissue2pu.out) == 'red blood cells'] <- 'rbc'
    names(Ktissue2pu.out) <- paste("K",names(Ktissue2pu.out),"2pu",sep='')
    names(vol) <- paste('V',names(vol),'c',sep='')
    names(flow)[names(flow) == 'liver'] <- 'total.liver'
    flow <- subset(unlist(flow), names(flow) != 'rest') / subset(physiology.data,Parameter=='Cardiac Output')[[species]]
    names(flow) <- paste('Q',names(flow),'f',sep='')
    
 	return(c(Ktissue2pu.out,vol,flow))
}