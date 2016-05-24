### Functions for creating encounter history files
# Landscape wrapper
createReplicates<-function(n_runs, map, Parameters, ... ){
  # 0. Match argument list 
  additional.args<-list(...)
    folder.dir<-setDefault(additional.args$folder.dir, getwd())    
    run.label<- setDefault(additional.args$run.label, 'rSPACE_X')
    base.name<- setDefault(additional.args$base.name, 'rSPACEx')
    filter.map<-additional.args$filter.map
    printN<-setDefault(additional.args$printN, 1)
    saveParameters<-setDefault(additional.args$saveParameters, 1)
    saveGrid<-setDefault(additional.args$saveGrid, 1)
    skipConfirm<-setDefault(additional.args$skipConfirm, F)
    add<-setDefault(additional.args$add, F)
    overwrite<-setDefault(additional.args$overwrite, F)

  # 0. Set up files
    if(!skipConfirm){
      askConfirm<-("" ==readline(prompt="\n rSPACE creates text files.
        If you're ok with this, press ENTER to continue.
        Typing anything else will exit.\n"))
      if(!askConfirm){
        message('Exiting function')
        return(0)
      }
   }

   folder.dir <- paste(folder.dir, run.label, sep='/')
   if(!file.exists(folder.dir)) dir.create(folder.dir)


   output.dir <- paste(folder.dir, 'output', sep='/')
   if(!file.exists(output.dir)) dir.create(output.dir)
   if(printN){
     printN<-paste0(output.dir,'/N_final.txt')
     if(overwrite){
      message(paste('Restarting', printN))
      file.remove(printN)
     }}

   prevPList<-paste0(output.dir,'/Parameters.rdata')
  
  # 1. Enter parameters
  if(missing(Parameters)) {
    if(file.exists(prevPList)){
      message(paste0('Using existing parameters list for scenario from: ', prevPList))
      load(prevPList)
    } else {
      if(add==T){
        stop('No parameter list available') 
      } else { Parameters<-enter.parameters()}
    }
  } else {
    if(add==T){
      if(file.exists(prevPList)){
        warning(paste0('Using existing parameters list for scenario from: ', prevPList))
        rm('Parameters')
        load(prevPList)
      }}
  }
  
  Parameters<-checkParameters(Parameters, additional.args)

  # 2. Set up map + grid layer
  if(missing(map)) stop("Missing habitat layer")
  map<-checkMap(map, filter.map)

  grid_layer<-createGrid(map, Parameters, filter.map)
  gridIDs<-unique(grid_layer)[unique(grid_layer)>0]

  # 2.5 Shift start for indices if needed
  n.prevFiles<-length(dir(folder.dir, pattern=base.name))
  rn.start<-ifelse(add, n.prevFiles, 0)
  if(add==F & n.prevFiles>0 & overwrite==F) 
    stop(paste('\nExisting rSPACE runs found in', folder.dir, '\n Use "overwrite=T" to replace or "add=T" to add to existing folder'))
    
  # 3. Simulate encounter histories loop ##
  for(rn in (1:n_runs)+rn.start){
    cat(rn,'\n');flush.console()
    ch<-encounter.history(map, Parameters, grid_layer=grid_layer, n_cells=length(gridIDs), printN=printN)
  
    # 4. Output encounter history
    output_file<-paste(folder.dir,'/',base.name,rn,".txt",sep='')
    cat(paste("/*", gridIDs, "*/", ch, "1;"), sep="\n",file=output_file) 
  } # End runs loop
  
  if(saveParameters)
    save(Parameters, file=paste0(output.dir,'/Parameters.Rdata'))
  if(saveGrid)
    writeRaster(setValues(map, grid_layer), filename=paste0(output.dir,'/Grid.tif'), overwrite=T)
    
  return(list(DIR = folder.dir, 
          filenames=paste0(base.name,(1:n_runs)+rn.start,'.txt')))
} # End function
create.landscapes<-function(...){
  .Deprecated("createReplicates")
  return(0)
  }