#For analysing single food webs

analyse.single <- function(filename, 
                           omn=c("FALSE", "TRUE"), 
                           cann=c("FALSE", "TRUE"),
                           positions=c("FALSE", "TRUE"),
                           matrix=c("FALSE", "TRUE"),
                           sp.names=c("FALSE", "TRUE"),
                           maxlevels=8) {
    
    omn <- match.arg(omn)
    cann <- match.arg(cann)
    matrix <- match.arg(matrix)
    sp.names <- match.arg(sp.names)
    positions <- match.arg(positions)
    
    #If there is a foodweb matrix loaded, erase it
    if (exists("foodweb")) {rm(foodweb, envir=globalenv())}
    
    
    if (is.numeric(maxlevels)==FALSE) {
      stop("Error: the argument 'maxlevels' should be an integer")
    }   
          
          #Upload the foodweb file
          foodweb <<- read.table(filename, header=FALSE, sep=',', colClasses=numeric()) 
          if (exists("foodweb")) {
            
            #Remove species names if they were supplied
            if (sp.names=="TRUE") {foodweb <<- foodweb[-(1), -(1)]}
            
            #Verify that the matrix contains only ones and zeroes
            if (any(foodweb[,]!=0 & foodweb[,]!=1)) {stop("Error: The food web matrix contains values other than 0 or 1")}
            
            #If the matrix is asymmetrical, create the symmetrical matrix
            if (ncol(foodweb)!=nrow(foodweb)) {
              if (sp.names==FALSE) stop("Error: You've provided an asymmetrical matrix without supplying species names or indicating sp.names=TRUE")
              asym2sym(foodweb=foodweb, name=filename, single=TRUE)
            } 
            
            #Run through the calculations
            name <<- filename
            maxlevels <<- maxlevels
            calculate.metrics(foodweb=foodweb, maxlevels=maxlevels, name=name)
            write.table(indices, file = paste("Results-", filename, sep=""), append=FALSE, quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)
            
            #Create Trophic positions file, if requested.   
            if (positions==TRUE) {
              write.table(levels.l, file = "Trophic positions.csv", append=FALSE, quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)  
            }
      
            #Create file with matrix used for calculations, if requested.
            if (matrix==TRUE) {
              write.table(foodweb[1:S,1:S], file = "Symmetrical matrix.csv", append=FALSE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)  
            }
      
            #Create file with list of omnivores, if requested
            if (omn==TRUE) {
              write.table(omn.l, file = "Omnivores.csv", append=FALSE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)  
            }
            
            #Create file with list of cannibals, if requested
            if (cann==TRUE) {
              write.table(omn.l, file = "Cannibals.csv", append=FALSE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)  
            }
            print(paste("Check your current working directory for the output file(s). Network parameters are in a file with the name Results-", filename, sep=""))
          } else {print("Problem uploading food web matrix file")}
}
  

#For analysing multiple food webs
analyse.seq <- function(prefix, suffix, first, last, maxlevels=8, 
                        sp.names=c("FALSE", "TRUE"),
                        sym=c("FALSE", "TRUE"),
                        output="Results.csv", 
                        problem="Problem.csv",
                        separator=",") {
    prefix <<- prefix
    suffix <<- suffix
    first <<- first
    last <<- last
    maxlevels <<- maxlevels
    output <<- output
    problem <<- problem
    separator <<- separator
    
    if (missing(sym)) sym <<- FALSE else sym <<- match.arg(sym)
    if (missing(sp.names)) sp.names <<- FALSE else sp.names <<- match.arg(sp.names)
    
    #Verify all necessary information was supplied
    if (exists("first") & exists("last") & exists("prefix") & exists("suffix")) {} else {
      stop("first, last, prefix and suffix must be provided in order to perform the calculations")}
    
    #Verify that "first" is an integer
    if (floor(first) - ceiling(first)!=0) {
        stop("Error: 'first' must be an integer")} 
    
    #Verify that "last" is an integer
    if (floor(last) - ceiling(last)!=0) {
        stop("Error: 'last' must be an integer")}
    
    #Verify that "last" is larger than "first"
    if (last<first) {
      stop("Error: 'Last' must be an integer larger than 'first'.")
    }
    
   

    #Perform the calculations iteratively and output
    for (i in first:last) {
        foodweb <- read.table(paste(prefix, i, suffix, sep=""), header=FALSE, sep=separator, colClasses=numeric())
               
        #Verify that the matrix is symmetrical
        if (sp.names=="TRUE") {
          asym2sym(foodweb=foodweb, problem=problem)
          if (sym==TRUE) write.table(foodweb, file = paste(prefix, i, "-sym.csv"), 
                                append=FALSE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE) 
          }
        

        name<- paste(prefix, i, suffix, sep="")
        
        #Verify that matrix has only zeroes and ones 
        if (any(foodweb[,]!=0 & foodweb[,]!=1)) {
           write.table(cbind(name, "Values other than 0 or 1 in matrix") , file = problem, append=TRUE, quote=FALSE, sep=",", 
            col.names=FALSE, row.names=FALSE) 
           next
        }
        
        calculate.metrics(foodweb=foodweb, maxlevels=maxlevels, name=name)
              if (i==first) {
                write.table(indices, file = output, append=FALSE, quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)
              } else {
                write.table(indices, file = output, append=TRUE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
              }
        rm(list = ls())
    }
    print(paste("Check your current working directory for the output file(s). If problems were dected in the matrices, a file called ", problem, " will have been created. For non-problematic files, the network parameters will be in a file named ", output, ".", sep=""))
}


# For analysing a list of food webs

analyse.list <- function(list, maxlevels=8, 
                         sym=c("FALSE", "TRUE"), 
                         output="Results.csv", 
                         problem="Problem.csv",
                         separator=",", 
                         sp.names=c("FALSE", "TRUE")) {
  
    if (missing(sp.names)) sp.names <<- FALSE else sp.names <<- match.arg(sp.names)
    
    if (missing(sym)) sym <<- FALSE else sym <<- match.arg(sym)
    
    files <<- read.table(list, header=FALSE, sep=',')
    separator <<- separator
    maxlevels <<- maxlevels
    problem <<- problem
    output <<- output
    
    write.table(cbind("Web name", "Problem"), file = problem, append=FALSE, quote=FALSE, sep=",", 
                col.names=FALSE, row.names=FALSE)
    
    for (i in 1:nrow(files)) {
        foodweb <<- read.table(as.character(files[i,]), sep=separator, header=FALSE)
        
        name <<- as.character(files[i,])
        
        if (sp.names==TRUE) {
          asym2sym(foodweb=foodweb, problem=problem, name=name)
          if (sym==TRUE) write.table(foodweb, file = paste(unlist(strsplit(name, split=".", fixed=TRUE)), "-sym.csv"), 
                    append=FALSE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE) 
          
        }
        
        #Verify that matrix has only zeroes and ones 
        if (any(foodweb[,]!=0 & foodweb[,]!=1)) {
           write.table(cbind(name, "Values other than 0 or 1 in matrix"), file = problem, append=TRUE, quote=FALSE, sep=",", 
            col.names=FALSE, row.names=FALSE)
           next
        }
        
        calculate.metrics(foodweb=foodweb, maxlevels=maxlevels, name=name)
              if (i==1) {
                write.table(indices, file = output, append=FALSE, quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)
              } else {
                write.table(indices, file = output, append=TRUE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
              }
        rm(list = ls())
    }
    
    print(paste("Check your current working directory for the output file(s). The network parameters are in a file named ", output, ". If problems were dected in the matrices, a file called ", problem, " will have been created.", sep=""))
}


#For plotting a network

plotweb <- function(cols, radii) {  

  levels <- length(unique(t((foodweb[S+1,]))))

  if (exists("foodweb") & exists("S")) {
       if (length(cols)==levels || missing(cols)) {
        if (length(radii)==levels || missing(radii)) {
          create.plot(foodweb=foodweb, radii=radii, cols=cols, levels=levels)} else {
           print(paste("Error: you have not supplied the adequate number of radii for plotting. There are ",  levels, " in your food web."))
              }
      } else {print(paste("Error: you have not supplied the adequate number of colours for plotting. There are ", levels, "trophic levels in your food web."))
              }
    } else {print("Error: there is no food web information available to create the plot. Use the analyse.single() function and try again.")}
  }

# For creating a list of trophic interactors from a food web matrix
mat.2.list <- function(filename, output) {

matrix <- read.table(filename, header=TRUE, sep=',')
Spred <- ncol(matrix)
Sprey <- nrow(matrix)
foodweb <- rbind("Foodweb")
link <- ""

for (PredatorID in 1:Spred) { 
    for (PreyID in 1:Sprey) { 
      if (matrix[PreyID, PredatorID] == "1"){
      link <- paste(PredatorID, PreyID)
      foodweb <- rbind(foodweb, link)
      }
    }
}

write(foodweb[2:nrow(foodweb),],file = output, append=FALSE)

}
