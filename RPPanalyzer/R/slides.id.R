`slides.id` <-
function(x){

## x = slidedescription data.frame from read.slidedescription()
            ## look for gpr files described in slidedescription table
               slides  <- unique(x[,"gpr"])
            ## look for gprfiles in current working directory   
               gprfiles <- dir(pattern="gpr$")
            ## check if all gprs found are described in slidedescription table   
                     proove <- c(gprfiles %in% slides)
                     
               if (all(proove)){
               print(paste("found",length(slides),"slides"))
               }
               else{
               print(paste("identified only",length(slides),"of",
               length(gprfiles),"gpr files saved in current folder"))
               }
               
            ## return vector with gpr filenames from slidedescription table
               return(as.character(slides))
 }

