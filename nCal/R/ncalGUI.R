# Author: Krisztian Sebestyen

file.ext <- function(x)rev(strsplit(x,"\\.")[[1]])[1]
file.stem <- function(x)strsplit(x,"\\.")[[1]][1]

## toolkits: could not get a constintent behaviour across different drivers so settled on tcltk
## tcltk  : intuitive behavior, mature

## RGtk2  : nice but widgets placed in sub-containers do not behave as intuitively as tcltk
## rJava  : nice but widgets placed in sub-containers do not behave as intuitively as tcltk
## arguments 'horizontal' and 'anchor' need to be specified in ggroup() ?

## RwxWidgets : wxWidgets is 116 Mb! (http://www.wxwidgets.org/downloads/#latest_stable)
## Qt     : failed to build (https://r-forge.r-project.org/R/?group_id=454)

## other toolkits have a nice look but do not behave as intuitively as tcltk
## for both RGtk2,rJava more care is needed to place 'state widgets'


    
## uses the gedit as a status bar - no extra window
ncalGUI <- function(verbose = FALSE){
    
#    x <- 1
#    env0 <- env <- "ncal.gui.env"
#    while(1){
#        if(inherits(try(get(env,envir = .GlobalEnv),silent = TRUE),"try-error")) break
#        x <- x + 1
#        env <- env0 %+% x
#    }
#    
#    assign(env,environment(),envir = .GlobalEnv)
#    if(verbose){
#        cat("local gui environment can be found in variable ",env," in global environment and it is: \n");
#        print(environment())
#    }
    
    #require(gWidgetstcltk) || stop('Please install the gWidgetstcltk package')
    options(guiToolkit = "tcltk")
    

    
    state <- NULL
    state$labels <- list(dependent = "             Response variable                   ",
                        predictor  = "             Concentration/dose variable ",
                        bcrm.fit = c("drm","robust Bayesian hierarchical model"),
                        input.data.type = c("Luminex raw output .xls","comma-separated file","tab-delimited file"),
                        input.data = "             Data file",
                        output.data = "             Save estimated concentrations to",
                        output.graph ="             Save calibration plots to                ")
                         
                         
    state$args <- list( dependent = "log(fi)",
                        predictor = "expected_conc",
                        bcrm.fit = FALSE,
                        input.data = "Select a file...",
                        output.data = path.expand("Estimated_concentrations.csv"),
                        output.graph = path.expand("Calibration_plots.pdf"))
                        
    # the handler list is the list of 'handlers' to access a widget's value
    state$handlers <- NULL  

    state$handlers$w <- gwindow("nCal")

    # the container list is a list of 'handlers' of a widget associated with where the widget is contained
    # it should rather be called 'location.handler' or sg like that 
    # initially I intended to keep them around in case extra info needs to be shown for a group of widgets
    # but probably not necessary
    state$containers <- NULL

    glabel("Input:",container = ggroup(container = state$handlers$w))

    # name of input dataset
    state$containers$input.data <- ggroup(container = state$handlers$w)
    glabel(state$labels$input.data,container = state$containers$input.data)
    state$handlers$input.data <- gfilebrowse (text = state$args$input.data, type = "open", quote = TRUE,container = state$containers$input.data, width=60)
    #svalue(state$handlers$input.data) <- state$args$input.data #initial value
    
#    glabel("Input data types are .xls .csv .txt:" ,container = state$handlers$w)    
#    glabel("Luminex raw output (.xls)",container = state$handlers$w)
#    glabel("Comma-separated file containing table only (.csv)",container = state$handlers$w)
#    glabel("Tab-separated file containing table only (.txt)",container = state$handlers$w)
#    glabel("",container = state$handlers$w)

    # glabel("Input data types are .xls .csv .tsv:" ,container = state$containers$input.data)   
    # glabel("Luminex raw output (.xls)",container = state$containers$input.data)
    # glabel("Comma-separated file containing table only (.csv)",container = state$containers$input.data)
    # glabel("Tab-separated file containing table only (.tsv)",container = state$containers$input.data)
    # glabel("",container = state$containers$input.data)

    # radio input dataset type
    state$containers$file.type <- ggroup(container = state$handlers$w)
    glabel("             Data file format    ",container = state$containers$file.type)
#        glabel("                 " ,container = state$containers$file.type)    
    state$handlers$input.data.type<-gradio(state$labels$input.data.type, 
           selected = 1, horizontal = TRUE,container = state$containers$file.type)           
    #svalue(state$handlers$input.data.type) <- state$labels$input.data.type[1] #initial value

    state$containers$dependent <- ggroup(container = state$handlers$w)
    glabel(state$labels$dependent,container = state$containers$dependent)
    state$handlers$dependent <- gedit(text = state$args$dependent, width = 32,initial.msg=state$args$dependent, container = state$containers$dependent) 
    svalue(state$handlers$dependent) <- state$args$dependent #initial value
    
    state$containers$predictor <- ggroup(container = state$handlers$w)
    glabel(state$labels$predictor,container = state$containers$predictor)
    state$handlers$predictor <- gedit(text = state$args$predictor, width = 32,initial.msg=state$args$predictor, container = state$containers$predictor) 
    svalue(state$handlers$predictor) <- state$args$predictor #initial value
     
    glabel("Output:",container = ggroup(container = state$handlers$w))

    # name of output dataset
    state$containers$output.data <- ggroup(container = state$handlers$w)
    glabel(state$labels$output.data,container = state$containers$output.data)
    state$handlers$output.data <- gfilebrowse (text = state$args$output.data, type = "save", quote = TRUE, width=40,
    container = state$containers$output.data,handler = function(h,...)svalue(h$obj) <- file.path(svalue(h$obj)))
#   svalue(state$handlers$output.data) <- state$args$output.data #initial value

    # name of output plots      
    state$containers$output.graph <- ggroup(container = state$handlers$w)
    glabel(state$labels$output.graph,container = state$containers$output.graph)
    state$handlers$output.graph <- gfilebrowse (text = state$args$output.graph, type = "save",quote = TRUE, width=40,
    container = state$containers$output.graph,handler = function(h,...)svalue(h$obj) <- file.path(svalue(h$obj)))
#   svalue(state$handlers$output.graph) <- state$args$output.graph #initial value

    # bcrm radio button
    glabel("Curve fitting methdod:",container = ggroup(container = state$handlers$w))
    state$containers$fit <- ggroup(container = state$handlers$w)
    glabel("            " ,container = state$containers$fit)    
    state$handlers$bcrm.fit<-gradio(state$labels$bcrm.fit, 
           selected = as.numeric(state$args$bcrm.fit) + 1, horizontal = TRUE,container = state$containers$fit)         
#   svalue(state$handlers$bcrm.fit) <- state$labels$bcrm.fit[as.numeric(state$args$bcrm.fit) + 1] #initial value

    call.getConc.fun <- function(h,...){
        state <- get("state",envir = h$action)
        
        ncal.fit <- try(get("gui.current.fit",envir = h$action))
        if(inherits(ncal.fit,"try-error")){
            enabled(state$handlers$sb) <- TRUE
            svalue(state$handlers$sb) <- c("No models were found, please fit a model first.")
            enabled(state$handlers$sb) <- FALSE 
            dispose(state$handlers$data.getConc)            
            return
        }
        ncal.fit <- attr(ncal.fit,"fits")[[1]]

        new.y <- svalue(state$handlers$data.getConc)
        new.y <- strsplit(new.y , ",")
        new.y <- paste(new.y[[1]],collapse = " ")
        new.y <- strsplit(new.y , " ")
        new.y <- as.numeric(new.y[[1]])     
        new.y <- new.y[!is.na(new.y)]
        ncal.pred <- try(do.call("getConc",list(ncal.fit,new.y)))

        if(inherits(ncal.pred,"try-error")){
            enabled(state$handlers$sb) <- TRUE
            svalue(state$handlers$sb) <- ncal.pred
            enabled(state$handlers$sb) <- FALSE 
            #dispose(state$handlers$data.getConc)           
            return
        }

        
        # add pop-up window here. note that window is created invisible until it is filled
        ncal.pred <- as.matrix(ncal.pred) 
        state$handlers$ncal.pred.window <- gwindow(visible = FALSE)
        state$handlers$ncal.pred.window.text <- gtext(container = state$handlers$ncal.pred.window,wrap = FALSE) 
        insert(state$handlers$ncal.pred.window.text,capture.output(ncal.pred),wrap = FALSE)
        visible(state$handlers$ncal.pred.window) <- TRUE
        #size(state$handlers$ncal.pred.window) <- as.integer(size(state$handlers$ncal.pred.window) * c(1.5,1) )
        
        #dispose(state$handlers$data.getConc) #clears text in gtext()
        assign("state",state,envir = h$action) #save state change
        return
    }
    

    update.status.fun <- function(h,...){
        state <- get("state",envir = h$action)
        ####################################################################################
        ## Do not proceed if user hasn't selected a file
        ## But gfile() forcing a selection is annoying
        status <- svalue(state$handlers$input.data)
        if(is.na(status) || (status == state$args$input.data)){
            enabled(state$handlers$sb) <- TRUE      
            svalue(state$handlers$sb) <- "Error: Please select an input file." 
            enabled(state$handlers$sb) <- FALSE         
            tcltk::tcl("update") 
            return()
        }
        ####################################################################################
        
        ####################################################################################        
        # update status bar - how do I trigger an update here(!) without returning ? 
        enabled(state$handlers$sb) <- TRUE      
        svalue(state$handlers$sb) <- "R: performing tasks ..."
        enabled(state$handlers$sb) <- FALSE     
        tcltk::tcl("update") 
        ####################################################################################
                
                                
        ####################################################################################
        ## set up ncal.arg
        # index <- match(intersect(names(formals(ncal.character)),names(state$handlers)),names(state$handlers),nomatch=0)                       
        # ncal.arg <- lapply(state$handlers[index],svalue)

        ncal.arg <- NULL
        ncal.arg$file <- svalue(state$handlers$input.data)

        
       #ncal.arg$is.luminex.xls <- (file.ext(ncal.arg$file) == "xls")
        # # data.type.index is in {1,2,3} for gradio if .xls,.csv,.tsv or '0' if something that should not have happened
        data.type.index <- match(svalue(state$handlers$input.data.type),state$labels$input.data.type,nomatch = 0)
        ncal.arg$is.luminex.xls <- (data.type.index == 1)
        ncal.arg$formula <- as.formula(paste(svalue(state$handlers$dependent),"~",svalue(state$handlers$predictor)))
        ncal.arg$bcrm.fit <- as.logical(match(svalue(state$handlers$bcrm.fit),state$labels$bcrm.fit,nomatch = 1) - 1)
        ncal.arg$return.fits <- TRUE
        ####################################################################################        
        mypdf(mfrow=c(3,4), file = file.stem(svalue(state$handlers$output.graph)))
        result <- try(do.call("ncal.character",ncal.arg))
        if(!inherits(result,'try-error')){
            assign("gui.current.fit",result,envir = h$action)
            write.csv(result,file = svalue(state$handlers$output.data),row.names = FALSE)
            enabled(state$handlers$sb) <- TRUE          
            svalue(state$handlers$sb) <- "Done. Waiting for user input ..."
            enabled(state$handlers$sb) <- FALSE         
            enabled(state$handlers$apply.getConc) <- TRUE #enable the APPLY-button          
        }else{
            enabled(state$handlers$sb) <- TRUE
            svalue(state$handlers$sb) <- result
            enabled(state$handlers$sb) <- FALSE
        }
        dev.off()
        ####################################################################################
                
    }

    # 'OK' or 'EXIT'
    ok_exit <- ggroup(container = state$handlers$w)
    addSpring(ok_exit)
    state$handlers$ok <- gbutton("GO", container = ok_exit, handler = update.status.fun,action = environment())
    addSpace(ok_exit, 40)
    #state$handlers$exit <- gbutton("Close", container = ok_exit,
    #                    handler = function(h, ...){dispose(state$handlers$w);try(dispose(state$handlers$ncal.pred.window));quit(save = "no")})
    glabel("  ",container = ggroup(container = state$handlers$w))

    
    glabel("Apply the last fitted curve to the new data below:",container = ggroup(container=state$handlers$w))
    
    state$containers$getConc.2 <- ggroup(container = state$handlers$w)
    addSpace(state$containers$getConc.2, 15)
    state$handlers$data.getConc <- gedit(text = "", width = 80,initial.msg="", container = state$containers$getConc.2 )

    state$containers$getConc <- ggroup(container = state$handlers$w)
    addSpring(state$containers$getConc)
    state$handlers$apply.getConc <- gbutton("APPLY", container = state$containers$getConc, handler = call.getConc.fun, action = environment())
    addSpace(state$containers$getConc, 40)
    enabled(state$handlers$apply.getConc) <- FALSE
    
    # unlike status bar, height is fixed here !                 
    state$handlers$sb <- gedit("Waiting for user input ...",container = state$handlers$w)
    enabled(state$handlers$sb) <- FALSE # user cannot type 
}


# deprecated
# state of gui is kept by a separate window
ncalGUI.2 <- function(toolkit = c("tcltk","RGtk2","rJava","RwxWidgets","Qt"),visible.state = FALSE){
    
    #require(gWidgets) || stop('Please install the gWidgets package')
    #require(nCal) || stop('Please install the nCal package')

    options(guiToolkit = toolkit[1])
    
    state <- NULL
    state$labels <- list(dependent = "response variable",predictor = "dose/concentration variable",
                        bcrm.fit = c("Fit drm","Fit robust Bayesian hierarchical model"),
                        input.data.type = c("Luminex raw output xls file","csv file containing table only","tab-delimited file containing table only"),
                        input.data = "Input data file",
                        output.data = "Save estimated concentrations to",
                        output.graph = "Save calibration plots to")
                         
                         
    state$args <- list( dependent = "log(fi)",
                        predictor = "expected_conc",
                        bcrm.fit = FALSE,
                        input.data = "Select a file...",
                        output.data = path.expand("Estimated_concentrations.csv"),
                        output.graph = path.expand("Calibration_plots.pdf"))
                        
    # the handler list is the list of 'handlers' to access a widget's value
    state$handlers <- NULL      
    state$handlers$w <- gwindow("nCal")
    state$handlers$sb <- gstatusbar("Waiting for user input ...",container = state$handlers$w)
    


    # the container list is a list of 'handlers' of a widget associated with where the widget is contained
    # it should rather be called 'location.handler' or sg like that 
    # initially I intended to keep them around in case extra info needs to be shown for a group of widgets
    # but probably not necessary
    state$containers <- NULL

    # name of input dataset
    state$containers$input.data <- ggroup(container = state$handlers$w)
    glabel(state$labels$input.data,container = state$containers$input.data)
    state$handlers$input.data <- gfilebrowse (text = state$args$input.data, type = "open", quote = TRUE,container = state$containers$input.data)
    #svalue(state$handlers$input.data) <- state$args$input.data #initial value

    # YOUYI COMMENT
    glabel("Input data types are .xls .csv .txt:" ,container = state$handlers$w)    
    glabel("Luminex raw output (.xls)",container = state$handlers$w)
    glabel("Comma-separated file containing table only (.csv)",container = state$handlers$w)
    glabel("Tab-separated file containing table only (.txt)",container = state$handlers$w)
    glabel("",container = state$handlers$w)
    # glabel("Input data types are .xls .csv .tsv:" ,container = state$containers$input.data)   
    # glabel("Luminex raw output (.xls)",container = state$containers$input.data)
    # glabel("Comma-separated file containing table only (.csv)",container = state$containers$input.data)
    # glabel("Tab-separated file containing table only (.tsv)",container = state$containers$input.data)
    # glabel("",container = state$containers$input.data)

    # YOUYI UNCOMMENT
    # radio input dataset type
    # state$handlers$input.data.type<-gradio(state$labels$input.data.type, 
           # selected = 1, horizontal = TRUE,container = ggroup(container = state$handlers$w))         
    # svalue(state$handlers$input.data.type) <- state$labels$input.data.type[1] #initial value

    state$containers$dependent <- ggroup(container = state$handlers$w)
    glabel(state$labels$dependent,container = state$containers$dependent)
    state$handlers$dependent <- gedit(text = state$args$dependent, width = 32,initial.msg=state$args$dependent, container = state$containers$dependent) 
#   svalue(state$handlers$dependent) <- state$args$dependent #initial value
    
    state$containers$predictor <- ggroup(container = state$handlers$w)
    glabel(state$labels$predictor,container = state$containers$predictor)
    state$handlers$predictor <- gedit(text = state$args$predictor, width = 32,initial.msg=state$args$predictor, container = state$containers$predictor) 
#   svalue(state$handlers$predictor) <- state$args$predictor #initial value
    
    # bcrm radio button
    state$handlers$bcrm.fit<-gradio(state$labels$bcrm.fit, 
           selected = as.numeric(state$args$bcrm.fit) + 1, horizontal = TRUE,container = ggroup(container = state$handlers$w))         
#   svalue(state$handlers$bcrm.fit) <- state$labels$bcrm.fit[as.numeric(state$args$bcrm.fit) + 1] #initial value
        
    # name of output dataset
    state$containers$output.data <- ggroup(container = state$handlers$w)
    glabel(state$labels$output.data,container = state$containers$output.data)
    state$handlers$output.data <- gfilebrowse (text = state$args$output.data, type = "save", quote = TRUE,container = state$containers$output.data)
#   svalue(state$handlers$output.data) <- state$args$output.data #initial value

    # name of output plots      
    state$containers$output.graph <- ggroup(container = state$handlers$w)
    glabel(state$labels$output.graph,container = state$containers$output.graph)
    state$handlers$output.graph <- gfilebrowse (text = state$args$output.graph, type = "save",quote = TRUE,container = state$containers$output.graph)
#   svalue(state$handlers$output.graph) <- state$args$output.graph #initial value


    update.status.fun <- function(h,...){
        tcltk::tcl("update") 

        ####################################################################################
        ## Do not proceed if user hasn't selected a file
        ## But gfile() forcing a selection is annoying
        status <- svalue(state$handlers$input.data)
        if(is.na(status) || (status == state$args$input.data)){
            svalue(state$handlers$sb) <- "Error: Please select an input file." 
            svalue(statew.handler) <- "wait" 
            return()
        }
        ####################################################################################
        
        # update status bar - how do I trigger an update here(!) without returning ? 
        svalue(state$handlers$sb) <- "R: performing tasks ..."
        svalue(statew.handler) <- "process"
        return()
        ####################################################################################
                
    }

    # calls ncal()
    process.fun <- function(h,...){

        tcltk::tcl("update") 

        # when switching from "process" -> "wait" the call to process.fun is still triggered
        # if(svalue(state$handlers$status) == "wait") return()
        if(svalue(get('statew.handler',envir = h$action)) == 'wait') return()
        
        # cat("R process starts \n");op <- options(digits.secs = 6);Sys.time();options(op)
        
        ####################################################################################
        ## Process dataset
        
        ####################################################################################
        ## set up ncal.arg
        # index <- match(intersect(names(formals(ncal.character)),names(state$handlers)),names(state$handlers),nomatch=0)                       
        # ncal.arg <- lapply(state$handlers[index],svalue)

        ncal.arg <- NULL
        ncal.arg$file <- svalue(state$handlers$input.data)
        
        # YOUYI COMMENT
        ncal.arg$is.luminex.xls <- (file.ext(ncal.arg$file) == "xls")
        
        # YOUYI UNCOMMENT
        # # data.type.index is in {1,2,3} for gradio if .xls,.csv,.tsv or '0' if something that should not have happened
        # data.type.index <- match(svalue(state$handlers$input.data.type),state$labels$input.data.type,nomatch = 0)
                                    
        ncal.arg$formula <- as.formula(paste(svalue(state$handlers$dependent),"~",svalue(state$handlers$predictor)))
        ncal.arg$bcrm.fit <- as.logical(match(svalue(state$handlers$bcrm.fit),state$labels$bcrm.fit,nomatch = 1) - 1)

        ####################################################################################        
        mypdf(mfrow=c(3,4), file = file.stem(svalue(state$handlers$output.graph)))
        result <- try(do.call("ncal.character",ncal.arg))
        if(!inherits(result,'try-error')){
            write.csv(result,file = svalue(state$handlers$output.data),row.names = FALSE)
            svalue(state$handlers$sb) <- "Done. Waiting for user input ..."
        }else{
            svalue(state$handlers$sb) <- result
        }
        dev.off()
        ####################################################################################
        # cat("R process ends \n");op <- options(digits.secs = 6);Sys.time();options(op)
        
        # does not get updated so instead destroy and recreate:
        # svalue(state$handlers$status) <- "wait"               
        dispose(get('statew',envir = h$action))
        statew <- gwindow(visible = visible.state)
        assign('statew',statew,envir = h$action)
        statew.handler <- gedit("wait",container = statew)
        assign('statew.handler',statew.handler,envir = h$action)
        svalue(statew.handler ) <- "wait"                   
        addHandlerChanged(statew.handler, handler=process.fun,action = h$action)

    }

    # trigger state update change via 'gedit'
    statew <- gwindow(visible = visible.state)
    statew.handler <- gedit("wait",container = statew)
    svalue(statew.handler) <- "wait"                    
    
    # env <- environment()
    # assign('statew',statew,envir = env)
    # assign('statew.handler',statew.handler,envir = env)
    addHandlerChanged(statew.handler, handler=process.fun,action = environment())


    
    # 'OK' or 'EXIT'
    ok_exit <- ggroup(container = state$handlers$w)
    addSpring(ok_exit)
    state$handlers$ok <- gbutton("Ok", container = ok_exit, handler = update.status.fun,action = environment())
    addSpace(ok_exit, 20)
    state$handlers$exit <- gbutton("Exit", container = ok_exit,
                        handler = function(h, ...){dispose(state$handlers$w);dispose(statew);quit(save = "no")})
}
    
