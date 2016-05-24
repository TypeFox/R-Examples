`read.slidedescription` <-
function (){

    ## check for slidedescription.txt in working directory
    allfiles <- dir()
    if ( ! "slidedescription.txt" %in% allfiles) {
        stop("can not find slidedescription.txt in current working directory.")
    }

    ## read in data frame
    incubation<-read.delim("slidedescription.txt",header=T)

    ## check header for required columns
    reqCols <- c("gpr","pad","slide","incubation_run","spotting_run","target")
    
    if(xor("foreground" %in% colnames(incubation), "background" %in% colnames(incubation))) {
        stop("A 'foreground' and 'background' column have to be specified together.")
    }

    if ( !all( reqCols %in% colnames(incubation))){
        stop("slidedescription file: columns are missing or header incorrect!")
    }
    ## check source plate description
    if ( mode(incubation[,"pad"])!="numeric" 
            | mode(incubation[,"slide"])!="numeric" 
            | mode(incubation[,"incubation_run"])!="numeric"
            | mode(incubation[,"spotting_run"])!="numeric"){
        stop ("slidedescription file: data format in columns pad, slide, incubation_run and spotting_run have to be numeric!") 
    }

    ## check if there is only one fcolumn and bcolumn given per gpr file
    gprFiles <- factor(as.character(incubation$gpr), levels=unique(as.character(incubation$gpr)))


    ## perform the check for the foreground and background check if given
    ## we ensure that there is only one column given for all pads of one slide
    if(("foreground" %in% colnames(incubation)) & ("background" %in% colnames(incubation))) {
        foreGroundCheck <- as.logical(tapply(incubation$foreground, gprFiles, function(xx) {

                            length(unique(as.character(xx))) !=1

                        }))
        if(any(foreGroundCheck)) {
            files <- levels(gprFiles)[foreGroundCheck]
            stop("For the gpr file(s) ", paste(files, collapse=","), " there are more than one foreground column specified.")
        }

        backGroundCheck <- as.logical(tapply(incubation$background, gprFiles, function(xx) {

                            length(unique(as.character(xx))) !=1

                        }))
        if(any(backGroundCheck)) {
            files <- levels(gprFiles)[backGroundCheck]
            stop("For the gpr file(s) ", paste(files, collapse=","), " there are more than one background column specified.")
        }

    }

    return(incubation)
}

