VariableSelection <- function(maxent,outdir,gridfolder,occurrencelocations,backgroundlocations,additionalargs,contributionthreshold,correlationthreshold,betamultiplier){
    dir.create(outdir)
    occurrencesites <- occurrencelocations
    backgroundsites <- backgroundlocations
    
                                        # the first three columns in the input tables must have the header 'species','longitude', and 'latitude'
    
                                        # Extracting the set of variables of the original csv table
    beginningvariableset <- colnames(read.csv(occurrencesites,header=TRUE))
                                        # Excluding the first three columns which contain no variables, but
                                        # the species name, longitudes and latitudes of occurrence sites
    beginningvariableset <- beginningvariableset[4:length(beginningvariableset)]
    beginning.variable.set <- beginningvariableset # the very first set of variables
        

        
                                        # Start a file that shows the single steps in the variable selection process
    cat(c("Test","Model","betamultiplier",beginningvariableset,"\n"), file = paste(outdir,"/VariableSelectionProcess.txt",sep=""), sep = "\t", fill = FALSE, labels = NULL, append = FALSE)
        
    cat(c("Model","betamultiplier","variables","samples","parameters","loglikelihood","AIC","AICc","BIC","AUC.Test","AUC.Train","AUC.Diff","\n"), file = paste(outdir,"/ModelPerformance.txt",sep=""), sep = "\t", fill = FALSE, labels = NULL, append = FALSE)

                                        # Number the models
    modelnumber <- 1

    for (b in betamultiplier){
        cat("----------------------- ","\n")
        cat("Choosing betamultiplier ", b,"\n")
        occurrencesites <- occurrencelocations
        backgroundsites <- backgroundlocations

        
                                        # Extract only the name and strip off the filepath from the occurrence and background site files
        occurrencesitefilename <- gsub(".*/","",occurrencesites)
        backgroundsitefilename <- gsub(".*/","",backgroundsites)
        
                                        #  Make a copy of the occurrence and backgroundsites in order not to overwrite the original files
        file.copy(from=occurrencesites,to= paste(outdir,"/",occurrencesitefilename,"_VariableSubset.csv",sep=""),overwrite=TRUE)
        file.copy(from=backgroundsites,to= paste(outdir,"/",backgroundsitefilename,"_VariableSubset.csv",sep=""),overwrite=TRUE)
        
        occurrencesites <- paste(outdir,"/",occurrencesitefilename,"_VariableSubset.csv",sep="")
        backgroundsites <- paste(outdir,"/",backgroundsitefilename,"_VariableSubset.csv",sep="")
        
        uncorrelated.variables <- character() # No variable testing was
                                        # yet done
        
        already.tested.variables <- character() # Creating an empty
                                        # character vector that
                                        # will be filled with the
                                        # names of variables for
                                        # which correlation to
                                        # other variables has
                                        # already been tested.
        
                                        #    while (length(beginning.variable.set)>length(uncorrelated.variables)){
        
                                        # As long as this is true, variables had been removed in the
                                        # last step and thus, more variables might be redundant.
        variablenames <- beginningvariableset        
        
        
        while(length(already.tested.variables)<length(variablenames)){
            cat("\n","Number of remaining variables",length(variablenames),"\n")
            cat("Testing variable contributions...","\n")
                                        # As long as there are still
                                        # important variables for
                                        # which correlation was not
                                        # yet tested, continue with the variable reduction steps
            
                                        # Maxentrun for AUC test value
                                        # calculation. The preset
                                        # options keep the
                                        # number of created files to a
                                        # minimum (only the
                                        # maxentResult.csv file is
                                        # needed) and to create with
                                        # every run a new set of
                                        # random test points.
            cat("Calculating average AUC values from 10 maxent models...","\n")
            Maxentrun(maxent = maxent,
                      outdir = outdir,
                      gridfolder = gridfolder,
                      occurrencesites = occurrencesites,
                      backgroundsites = backgroundsites,
                      additionalargs = paste("plots=false writeplotdata=false visible=false autorun=true randomseed=true writebackgroundpredictions=false replicates=10 replicatetype=subsample randomtestpoints=50 redoifexists writemess=false writeclampgrid=false askoverwrite=false pictures=false outputgrids=false -b ",b," ",
                          additionalargs,sep="")
                      )
            AUCs <- MaxentAUC(outdir)

            variablecontributions <- Subsetselection(outdir)
            # Asses the average contribution of each environmental variable


                                        #Maxentrun for information
                                        #criteria calculation. The
                                        #preset options reduce the
                                        #number of unnecessary output
                                        #files. Information criteria
                                        #calculation requires only raw
                                        #output grids, a lambdas file
                                        #and the latitude and
                                        #longitude data of occurrence
                                        #sites.
            cat("Calculating Information criteria...","\n")            
            Maxentrun(maxent = maxent,
                      outdir = outdir,
                      gridfolder = gridfolder,
                      occurrencesites = occurrencesites,
                      backgroundsites = backgroundsites,
                      additionalargs = paste("plots=false writeplotdata=false visible=false autorun=true randomseed=true writebackgroundpredictions=false redoifexists writemess=false writeclampgrid=false askoverwrite=false pictures=false outputformat=raw -b ",b," ",
                          additionalargs,sep="")
                      )
            
            occurrencesites.data <- read.csv(occurrencesites,header=TRUE)
            species <- as.character(occurrencesites.data$species[1])
            lambdas <- paste(outdir,"/",species,".lambdas",sep="")
            gridfolder.without.path <- gsub(".*/","",gridfolder)
            pred.raw <- paste(outdir,"/",species,"_",gridfolder.without.path,".asc",sep="")
           
            ICs <- MaxentIC( occurrencesites, pred.raw, lambdas)

            beginning.variable.set <- names(variablecontributions) # The set of
                                        # variables before
                                        # any is removed in this round
            
            cat(c(modelnumber,b,length(beginning.variable.set),as.vector(ICs),AUCs,abs(diff(AUCs)),"\n"), file = paste(outdir,"/ModelPerformance.txt",sep=""), sep = "\t", fill = FALSE, labels = NULL, append = TRUE)
            

            
            
            
                                        # Extracting the variable contributions for the set of all environmental variables
            matching.variables <- numeric(length(beginningvariableset))
            matching.variables[matching.variables==0] <- NA
            matching.variables[match(names(variablecontributions),beginningvariableset)] <- variablecontributions
            matching.variables <- unlist(matching.variables)
            
            cat(c("Contributions",modelnumber,b,matching.variables,"\n"), file = paste(outdir,"/VariableSelectionProcess.txt",sep=""), sep = "\t", fill = FALSE, labels = NULL, append = TRUE)
            
                                        # selecting the set of variables that exceed a user-defined contribution threshold
            selected.variables <- variablecontributions[,variablecontributions>contributionthreshold]
            
                                        # the names of selected variables
            variablenames <- colnames(selected.variables)
                                        # The first variable name belongs to the most important variable
                                        # because they were sorted by variable contribution
            
                                        # Exclude the variables for which correlation was already tested
            remaining.variables <- variablenames[variablenames%in%already.tested.variables==FALSE]
            
            important.notyettested.variable <- remaining.variables[1]
                                        # As the variablenames were sorted by importance, the first
                                        # remaining one has highest contribution to the model
            
                                        # add this variable to the set
                                        # of variables for which
                                        # correlation to other
                                        # variables was already thested
            already.tested.variables <- c(already.tested.variables,important.notyettested.variable)
                                        # Getting a list with  correlation coefficients of all selected variables and the names of uncorrelated variables
            cat("Testing correlations between relevant variables...","\n")
            correlations <- Correlations(important.notyettested.variable,variablenames,backgroundsites,correlationthreshold)
            
                                        # Extracting the correlation coefficients for the set of all environmental variables
            matching.correlations <- numeric(length(beginningvariableset))
            matching.correlations[matching.correlations==0] <- NA
            matching.correlations[match(names(correlations[[1]]),beginningvariableset)] <- correlations[[1]]
            matching.correlations <- unlist(matching.correlations)
            
            
            cat(c("Correlation",modelnumber,b,matching.correlations,"\n"), file = paste(outdir,"/VariableSelectionProcess.txt",sep=""), sep = "\t", fill = FALSE, labels = NULL, append = TRUE)
                                        # The names of uncorrelated
                                        # variables. This is the set
                                        # of remaining variables
            uncorrelated.variables <- correlations[[2]]
            
                                        # Creating new csv files that contain only the uncorrelated variables
            ExtractVariables(uncorrelated.variables,occurrencesites,backgroundsites)
            
                                        # Increase the model number by 1
            modelnumber <- modelnumber+1
            
        }
    }    
                                        # Transform the table that lists Variable selection processes
    SelectionProcessTable <- read.table(paste(outdir,"/VariableSelectionProcess.txt",sep=""),header=TRUE)
    
    SelectionProcessTable.transposed <- t(SelectionProcessTable)
    
    write.table(SelectionProcessTable.transposed, file = paste(outdir,"/VariableSelectionProcess.txt",sep="") , append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = FALSE)
    
    ModelWithBestPerformance(outdir=outdir) # Extract the models of
                                            # highest performance (and
                                            # the associated
                                            # contribution and
                                            # correlation tests) based
                                            # on AICc (minimum) and
                                            # AUC.Test (maximum)
                                            # values

    ModelSelectionVisualization(outdir=outdir) # Create plots that
                                               # show the range of
                                               # AICc and AUC.Test
                                               # values for a range of
                                               # betamultipliers and
                                               # number of Variables

                                        # List all files that were created in the outdir folder.
    all.outfiles <- list.files(outdir)
                                        # Remove the filenames from this list that shall remain
    all.outfiles <- all.outfiles[all.outfiles!="VariableSelectionProcess.txt"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelPerformance.txt"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelWithMaxAUCTest.txt"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelWithMinAICc.txt"]                  
    all.outfiles <- all.outfiles[all.outfiles!="VariableSelectionMaxAUCTest.txt"]       
    all.outfiles <- all.outfiles[all.outfiles!="VariableSelectionMinAICc.txt"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelSelectionAICc_MarkedMinAICc.png"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelSelectionAUCTest_MarkedMinAICc.png"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelSelectionAICc_MarkedMaxAUCTest.png"]
    all.outfiles <- all.outfiles[all.outfiles!="ModelSelectionAUCTest_MarkedMaxAUCTest.png"]
                                        # remove all the unneccessary output files
    setwd(outdir)
    file.remove(all.outfiles)
}
