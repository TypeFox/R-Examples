# This function displays the information specified in "info=" for all chemicals for which the PBPK model can be paramterized.
get_cheminfo <- function(info="CAS",species="Human",exclude.fub.zero=NA,fub.lod.default=0.005,model='3compartmentss',default.to.human=F)
{ 
  chem.physical_and_invitro.data <- chem.physical_and_invitro.data
  if(tolower(species) == 'human') species <- 'Human' 
  else if(tolower(species) == 'rat') species <- 'Rat'
  else if(tolower(species) == 'dog') species <- 'Dog'
  else if(tolower(species) == 'rabbit') species <- 'Rabbit'
  else if(tolower(species) == 'mouse') species <- 'Mouse'
  else stop("Only species of human, rat, mouse, rabbit, and dog accepted.")
  
  if(default.to.human==T) species <- 'Human'
  
  model <- tolower(model)
  if (model == "pbtk" | model == "3compartment" | model == "1compartment")
  {
    necessary.params <- c(paste(species,"Clint",sep="."),paste(species,"Funbound.plasma",sep="."),"MW","logP")
    if (is.na(exclude.fub.zero)) exclude.fub.zero=T
  }
  else if (model == "3compartmentss")
  {
    necessary.params <- c(paste(species,"Clint",sep="."),paste(species,"Funbound.plasma",sep="."),"MW","logP")

    if (is.na(exclude.fub.zero)) exclude.fub.zero <- F

  }else if(model == 'schmitt'){  
    
    necessary.params <- c(paste(species,"Funbound.plasma",sep="."),"logP")
    if (is.na(exclude.fub.zero)) exclude.fub.zero=T
    
  }else stop("Valid models are currently only: pbtk, 1compartment, 3compartment, schmitt, and 3compartmentss.")
  
  good.chemicals.index <- apply(chem.physical_and_invitro.data[,necessary.params],1,function(x) all(!is.na(x)))
  if (exclude.fub.zero) good.chemicals.index <- good.chemicals.index & (chem.physical_and_invitro.data[,paste(species,"Funbound.plasma",sep=".")]>0) 
  good.chemical.data <- chem.physical_and_invitro.data[good.chemicals.index,] 

    if('mw' %in% tolower(info)) info <- c('MW',info[tolower(info) != 'mw'])
    if('pka_accept' %in% tolower(info)) info <- c('pKa_Accept',info[tolower(info) != 'pka_accept'])
    if('pka_donor' %in% tolower(info)) info <- c('pKa_Donor',info[tolower(info) != 'pka_donor'])
    if('logp' %in% tolower(info)) info <- c('logP',info[tolower(info) != 'logp'])
    if('compound' %in% tolower(info)) info <- c('Compound',info[tolower(info) != 'compound'])
    if('cas' %in% tolower(info)) info <- c('CAS',info[tolower(info) != 'cas'])
  
  valid.info <- c("Compound","CAS","logP","pKa_Accept","pKa_Donor","MW","Clint","Clint.pValue","Funbound.plasma")

  if (any(toupper(info)=="ALL")) info <- valid.info
  
  if (any(!(toupper(info) %in% toupper(valid.info)))) stop(paste("Data on",info[!(info %in% valid.info)],"not available. Valid options are:",paste(valid.info,collapse=" ")))

  if (toupper("Clint") %in% toupper(info)) info[toupper(info)==toupper("Clint")] <- paste(species,"Clint",sep=".")
  if (toupper("Clint.pValue") %in% toupper(info)) info[toupper(info)==toupper("Clint.pValue")] <- paste(species,"Clint.pValue",sep=".")
  if (toupper("Funbound.plasma") %in% toupper(info)) info[toupper(info)==toupper("Funbound.plasma")] <- paste(species,"Funbound.plasma",sep=".")
 
 #if (!toupper(paste(species,"Clint",sep=".")) %in% toupper(colnames(chem.physical_and_invitro.data))) stop(paste("Species",species,"not found."))
    
  columns <- colnames(chem.physical_and_invitro.data)
#    c("Compound","CAS","MW",paste(species,"Clint",sep="."),paste(species,"Fub",sep="."))
  this.subset <- good.chemical.data[,toupper(colnames(chem.physical_and_invitro.data))%in%toupper(columns)]
  
  if('CAS' %in% info) rownames(this.subset) <- NULL 
  
 # this.subset <- suppressWarnings(this.subset[!is.na(as.numeric(this.subset[,paste(species,"Clint",sep=".")])),])
  #this.subset <- suppressWarnings(this.subset[!is.na(as.numeric(this.subset[,paste(species,"Fub",sep=".")])),])
  if (exclude.fub.zero) this.subset <- suppressWarnings(this.subset[as.numeric(this.subset[,paste(species,"Funbound.plasma",sep=".")])!=0,])  
  else this.subset[suppressWarnings(as.numeric(this.subset[,paste(species,"Funbound.plasma",sep=".")]) == 0),paste(species,"Funbound.plasma",sep=".")] <- fub.lod.default
                              
  data.table <- this.subset
  
   
    
  return(data.table[,info])
}