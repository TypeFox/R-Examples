
##########################
##########################
##common.corrections
##########################
##########################

#kr

#description
##########################
#functions to apply corrections
#so take an input, modify and 
#(unless told not to) 
#save it as the same thing


#includes 
##########################
#correctInput
#correctPitotDrift
#


#
#most urgent
######################################
#

#urgent
##########################
#

#to do
##########################
#

#comments
##########################
#



##########################
##########################
##correctInput
##########################
##########################

#kr 23/11/2013 v 0.0.1

#what it does
##########################
#takes an input
#applies a correction
#function
#return that as input 
#assuming overwrite = TRUE



###############################
###############################
##calcPack2
###############################
###############################

calcPack2 <- function(input, ..., settings = NULL, data = NULL){

    #run checks
    extra.args <- listUpdate(list(this.call=match.call(), fun.name = "calcPack2", 
                                  overwrite=TRUE), 
                             list(...))
    att <- attributes(input)

    if(settings$output=="input")
        return(input)

    data <- checkPEMS(data, output="pems")

    if(settings$overwrite==FALSE){
        temp <- make.names(c(names(data), attributes(input)$name), unique=TRUE)
        attributes(input)$name <- temp[length(temp)]
    } else {
        attributes(input)$name <- make.names(attributes(input)$name)
    }
    att$name <- attributes(input)$name

    old.class <- class(data)
    class(data) <- "not.pems"

    if(length(input)<nrow(data$data)){
        input <- c(input, rep(NA, nrow(data$data)-length(input)))
        attributes(input) <- att
    }

    data$data[1:length(input), attributes(input)$name]<-input

    data$units[1, attributes(input)$name] <- if("units" %in% names(attributes(input)))
                                                attributes(input)$units else NA
   
    if("history"  %in% names(data))
        data$history <- c(data$history, extra.args$this.call)
    
    class(data) <- old.class

    if(settings$output=="data.frame")
        return(pemsData(data)) else return(data)
    
}


################################
################################
##correctInput
################################
################################


correctInput <- function(input = NULL, ..., data = NULL,
         correction = NULL){

    #run checks
    extra.args <- listUpdate(list(this.call=match.call(), fun.name = "correctInput", 
                                  overwrite=TRUE), 
                             list(...))
    settings <- do.call(calcChecks, listUpdate(list(data=data), extra.args))

#this could be first check in
#checkInput? then hijack would be 
#hidden

    if(!"hijack" %in% names(extra.args) || extra.args$hijack == FALSE){   
        input <- checkInput(input, data=data, if.missing = "stop", fun.name = extra.args$fun.name)  
    }
    att <- attributes(input)
    temp <- try(names(formals(correction)), silent=TRUE)
    if(class(temp)[1]=="try-error") 
        stop("In ", extra.args$fun.name, ": problem with correction", 
             call. = FALSE)
    ignore <- if("..." %in% temp)
                  NULL else names(extra.args)[!names(extra.args) %in% temp]
    temp2 <- list(input=input)
    if(!"input" %in% temp)
        names(temp2)[1] <- temp[1] 
    ans <- try(do.call(correction, listUpdate(temp2, extra.args, ignore=ignore)), silent=TRUE)
    if(class(ans)[1]=="try-error")
        stop("In ", extra.args$fun.name, ": problem with correction", 
             call. = FALSE)
    attributes(ans) <- att

    calcPack2(input=ans, settings=settings, data=data, 
              this.call=extra.args$this.call, fun.name=extra.args$fun.name)
 
}





################################
################################
##zeroNetagive
################################
################################



zeroNegatives <- function(input = NULL, ..., data = NULL, screen = FALSE){

    #run checks
    extra.args <- listUpdate(list(this.call=match.call(), fun.name = "zeroNegatives", 
                                  overwrite=TRUE), 
                             list(...))
    settings <- do.call(calcChecks, listUpdate(list(data=data), extra.args))
    if(!"hijack" %in% names(extra.args) || extra.args$hijack == FALSE){   
        input <- checkInput(input, data=data, if.missing = "stop", fun.name = extra.args$fun.name)  
    }

    #get current attributes
    att <- attributes(input)

    #main calculation
    ans <- ifelse(input<0, 0, input)

    if(screen){

#proposed plot

#loaPlot(a~1:1000*pems.1$conc.co, panel=panel.compareZcases, 
#scheme="kr.blues", col.regions="Reds", line.col="darkblue")

#



        index <- 1:length(input)
        temp.panel <- function(x=x, y=y, z=z,...){
                            panel.xyplot(x=x, y=y, col="black", type="l", ...)
                            panel.xyplot(x=x, y=z, col="red", type="l", ...)
                            panel.xyplot(x=x, y=y, col=ifelse(y==z, NA, "red"),...)
                      }
        plot.list <- list(x=ans~index*input, grid=TRUE, key=FALSE,
                          panel=temp.panel)
        plot.list <- listUpdate(plot.list, list(...))
        print(do.call(loaPlot, plot.list))

        #accept, discard or rework plot option
        #accept send you on
        #discard stops here
        #rework allows you to redo

    }



    #to use other function to do something similar
    #ans <- correctInput(input = input, correction = function(x) x <- ifelse(x<0,0,x), 
    #                    hijack=TRUE, this.call=extra.args$this.call, 
    #                    fun.name=extra.args$fun.name)
    
    #note: we pass on fun.name, this.call and hijack

    #transfer attributes
    attributes(ans) <- att

    calcPack2(input=ans, settings=settings, data=data, 
              this.call=extra.args$this.call, fun.name=extra.args$fun.name)
 
}


