dataReading <-
function(IPfile, controlFile, yeastSpecies = NULL){

#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######

    print("**********************************************")
    print("Reading IP and control datasets... ")
 
    IPdata = read.table(IPfile)
    controlData = read.table(controlFile)

    print("... done")

    cdsPositions = yeastSpecies
    
    print("**********************************************")
    print("")

    return(list(IPdata = IPdata, controlData = controlData, cdsPositions = cdsPositions))

# end of function dataReading()
}
