#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: siena01.r
# *
# * Description: This module contains the code for the gui for creation of a
# * Siena data object.
# *****************************************************************************/
##@installGui Miscellaneous
installGui <- function()
{
    if (.Platform$OS.type =="windows")
    {
        message("The standalone gui is no longer available. ",
                "To run the gui within R:\n",
                "Start R, load RSiena using Packages menu ",
                "then type siena01Gui() (and <ENTER>)")
    }
    else
    {
        stop("Gui only needs installing on Windows: on Linux use sienascript")
    }
}

##@siena01Gui siena01
siena01Gui <- function(getDocumentation=FALSE)
{
   ##DONE (FALSE) ## this is so we can exit cleanly, but seems redundant here
    maxDegree <- NULL
    nMaxDegree <- NULL
    resultsFileID <-  NULL
    estimAns <- NULL
    modelName <- NULL
   ## effEdit <- NULL
    noFiles <- 0
    filename <- NA
    files <- NULL
    fileEditFlag <- FALSE
    depvarnames <- NULL
    ndepvars <- 0
    nettypes <- NULL
    estimVar <- NULL
    effectsVar <-  NULL
    condVar <- NULL
    gainVar  <-  NULL
    stdstartVar <- NULL
    ph2spinVar  <-  NULL
    rsspinVar <- NULL
    rsVar <- NULL
    clustVar  <-  NULL
    clustspinVar <- NULL
    derivVar <- NULL
    ph3spinVar <- NULL
    maxdfVar <- NULL
    session <- NULL
    defaults <- c("Group1","Var1", "","", " ", "Actors",
                  "","Yes", "NA", "1", " ")
    mydata <- NULL
    myeff <- NULL
    mymodel <-  NULL
    ##@addFile internal siena01Gui
    addFile<- function()
    {
        noFiles <<- noFiles+1
        addFiletypes <- paste("{{Data Files} .dat}",
                              "{{Pajek network files} .net}",
                              "{{All files} *}")
        filename[noFiles] <<-
            basename(tclvalue(tkgetOpenFile(filetypes=addFiletypes,
                                            initialdir=initialDir)))
        if (filename[noFiles] == "")
        {
            noFiles <<- noFiles - 1
            return()
        }
        if (tableRows < noFiles)
        {
            addTableRow(tableRows+1)
            tableRows <<- tableRows + 1
        }
        mytkarray[[noFiles, 3]] <<- as.tclObj(format(filename[noFiles],
                                                     length=50), drop=TRUE)
        mytkarray[[noFiles, 6]] <<- "Actors"
        mytkarray[[noFiles, 8]] <<- "Yes"
        mytkarray[[noFiles, 5]] <<- noFiles ## period
        if (substring(filename[noFiles], nchar(filename[noFiles]) - 3,
                      nchar(filename[noFiles]))=='.net')
            tkset(formatspins[[noFiles]], "pajek net")
        tcl(table1, "selection", "clear", "all") ## unselect everything
        tcl(table1, "selection", 'set', paste(noFiles,',3', sep=''))
        tcl(table1, "yview", noFiles)
        ## find the directory we are in and use it next time
        initialDir <<- dirname(filename[noFiles])
    }

    ##@addTableRow internal siena01Gui
    addTableRow <- function(i)
    {
        tkinsert(table1,"rows","end","1")
        ##create spinbox for format
        formatspins[[i]] <<- tkwidget(table1, 'spinbox', state='readonly',
                                      width=20, values=ff, cursor="arrow")
        mypos <- paste(i, ',',4, sep='')
        tkwindow.configure(table1, mypos, window=formatspins[[i]])
        tkbind(formatspins[[i]], "<FocusIn>",
               function(x)tcl(table1,"activate",'1,4'))
        ##create spinbox for type
        typespins[[i]] <<- tkwidget(table1, 'spinbox', state='readonly',
                                    width=25, values=typelist, cursor="arrow")
        mypos <- paste(i, ',',7, sep='')
        tkwindow.configure(table1, mypos, window=typespins[[i]])
        tkbind(typespins[[i]], "<FocusIn>",
               function(x) tcl(table1,"activate",'1,4'))
        tkconfigure(table1, height=i+1)
    }

    ##@applyFn internal siena01Gui
    applyFn <- function() ## prompt to save, then try to create data, then
        ## sienaModelOptions
    {
        if (noFiles == 0)
        {
            tkmessageBox(message='No data to apply', icon='error')
            return()
        }
        ans <- tkmessageBox(message='Do you want to save the session?',
                            type='yesno', icon='question')
        if (tclvalue(ans)=='yes')
        {
            saveFn()
        }
        else
        {
            sessionFromTcl()
            if (is.null(modelName))
            {
                modelName <<- "Siena"
            }
        }
        if (inherits(resp <-
                     try(sienaDataCreateFromSession(session=session,
                                                    modelName=modelName,
                                                    edited=fileEditFlag,
                                                    files=files),
                         silent=TRUE), "try-error"))
        {
            tkmessageBox(message=resp, icon='error')
        }
        else
        {
            mydata <<- resp$mydata
            myeff <<- resp$myeff
            mymodel <<- sienaAlgorithmCreate()
            savedObjectName <- paste(modelName, ".Rdata", sep="")
            save(mydata, myeff, mymodel, file=savedObjectName)
            sienaModelOptions()
        }
    }
    ##@clearFn internal siena01Gui
    clearFn <- function()
    {
        noFiles <<- 0
        filename <<- NULL
        for (i in 1:tableRows)
            for (j in 1:11)
                mytkarray[[i,j]] <<- NULL
        lapply(typespins, function(x) tkset(x, 'network'))
        lapply(formatspins, function(x) tkset(x,'matrix' ))

    }
    ##@deleteTableRow internal siena01Gui
    deleteTableRow <- function(i)
    {
        ## unmap the format window from the table
        mypos <- paste(i,',4', sep='')
        tkwindow.configure(table1, mypos, window="")
        ## delete the window
        tcl(table1, 'window', 'delete', mypos)
        ## remove the tcl variable behind it
        formatspins <<- formatspins[-i]
        ## unmap the type window from the table
        mypos <- paste(i,',7', sep='')
        tkwindow.configure(table1, mypos, window="")
        ## delete the window
        tcl(table1, 'window', 'delete', mypos)
        ## remove the tcl variable behind it
        typespins <<- typespins[-i]
        ## delete the row from the table
        tkdelete(table1, 'rows', i, '1')
        tableRows <<- tableRows - 1
        ##  sessionFromTcl()
        if (noFiles > 0)
        {
            noFiles <<- noFiles - 1
            files <<- files[-i]
        }
        else
            files <<- NULL
    }
    ##@editFile internal siena01Gui
    editFile<- function()
    {
        ##try: may be nothing selected or a box beneath spinbox
        selcursor <- tclvalue(tcl(table1, 'curselection'))
        if (selcursor == "")
        {
            tkmessageBox(message="No file selected")
            return()
        }
        else
        {
            fileno <- as.numeric(strsplit(selcursor, ',')[[1]][1])
            sessionFromTcl()
            files <<- readInFiles(session, fileEditFlag, files)
            tmpfile <- files[[fileno]]
            files[[fileno]] <<- edit(tmpfile) ## may need to undo
            fileEditFlag[fileno] <<-  TRUE
        }
        tkfocus(tt)
        ## put on top globally temporarily
        tcl('wm', 'attributes', tt, '-topmost', 1)
        Sys.sleep(0.1)
        tcl('wm', 'attributes', tt, '-topmost', 0)
        invisible()
    }
    ##@fromFileFn internal siena01Gui
    fromFileFn <- function()
    {

        sessionFiletypes <- paste("{{Text Files} {.txt .csv .prn}}",
                                  ## " {{Excel files} .xls}",
                                  "{{All files} *}")
        loadfilename <- tclvalue(tkgetOpenFile(filetypes =
                                               sessionFiletypes))
        ## browser()
        if (loadfilename == "")
        {
            return(FALSE)
        }
        modelName <<- basename(loadfilename)
        ipos <- max(c(0, gregexpr('.', modelName, fixed=TRUE)[[1]]))
        if (ipos > 1)
        {
            modelName <<- substring(modelName, 1, (ipos - 1))
        }
        session <<- sessionFromFile(loadfilename, tk=TRUE)
        procSession()
        TRUE
    }
    ##@fromFileContFn internal siena01Gui
    fromFileContFn <- function()
    {
        OK <- fromFileFn()
        if (OK)
        {
            ## try to read in the project object
            savedModelName <- paste(modelName, ".Rdata", sep='')
                                        #browser()
            if (inherits(try(load(savedModelName), silent=TRUE), "try-error"))
            {
                tkmessageBox(message="Unable to load saved model", icon="error")
            }
            else
            {
                mydata <<- mydata
                mymodel <<- mymodel
                myeff <<- myeff
                sienaModelOptions()
            }
        }
    }
    ##@helpFn internal siena01Gui
    helpFn <- function() ## display the manual
    {
        RShowDoc("RSiena_Manual", package=pkgname)
    }
    ##@myStop internal siena01Gui
    myStop<- function()
    {
        if (!DONE() && exists("mydata") && exists("myeff") &&
            exists("mymodel") && !is.null(mydata) && !is.null(myeff) &&
            !is.null(mymodel))
        {
            ans <- tkmessageBox(message='Do you want to save the model?',
                                type='yesno', icon='question')
            if (tclvalue(ans)=='yes')
            {
                savefileFn()
            }
        }
        tkdestroy(tt)
        DONE(TRUE)
    }
    ##@procSession internal siena01Gui
    procSession <- function(replace=FALSE) ##
    {
        if (replace)
        {
            if(tableRows != nrow(session))
                browser()
        }
        if (!replace)
        {
            if (tableRows < nrow(session))
                for (i in (tableRows + 1) :(nrow(session)))
                    addTableRow(i)
            else if (tableRows > nrow(session))
                for (i in (nrow(session) + 1) : tableRows)
                    deleteTableRow(i)
        }
        for (i in 1:nrow(session))
        {
            for (j in 1: ncol(session))
                mytkarray[[i, j]] <<- as.tclObj(session[i,j], drop=TRUE)
            tkset(formatspins[[i]], session[i,4])
            tkset(typespins[[i]], session[i, 7])
            filename[[i]] <<- session[i, 3]
        }
        tableRows <<- nrow(session)
        noFiles <<- tableRows
        tcl(table1, "selection", "clear", "all") ## unselect everything
        tcl(table1, "activate", "1, 3")
        tcl(table1, "selection", 'set', paste('1', ',3', sep=''))
    }
    ##@removeFile internal siena01Gui
    removeFile <- function()
    {
        selcursor <- tclvalue(tcl(table1, 'curselection'))
        fileno <- as.numeric(strsplit(selcursor, ',')[[1]][1])
        if (is.na(fileno) || !is.numeric(fileno))
        {
            tkmessageBox(message='No file selected to remove')
            return()
        }
        session <<- NULL
        deleteTableRow(fileno)
        tcl(table1, "selection", "clear", "all") ## unselect everything
        tcl(table1, "activate", '1, 4')
    }
    ##@saveFn internal siena01Gui
    saveFn <- function() ## saves session file
    {
        if (noFiles == 0)
        {
            tkmessageBox(message='No data to save')
            return()
        }
        sessionFromTcl()
        sessionFiletypes <- "{{csv file} *.csv}"
        if (!is.null(modelName))
        {
            init <- modelName
        }
        else
        {
            init <- "Siena"
        }
        savefilename <- tclvalue(tkgetSaveFile(filetypes=sessionFiletypes,
                                               defaultextension='.csv',
                                               initialfile=init))
        if (savefilename != "")
        {
            write.table(session, file=savefilename, sep=',', row.names=FALSE)
        }
        modelName <<- basename(savefilename)
        ipos <- max(c(0, gregexpr('.', modelName, fixed=TRUE)[[1]]))
        if (ipos > 1)
        {
            modelName <<- substring(modelName, 1, (ipos - 1))
        }
    }
    ##@savefileFn internal siena01Gui
    savefileFn <- function() ## saves data and model
    {
        mymodel <<- modelFromTcl()
        modelFiletypes <- "{{R object} *.Rdata}"
        if (!is.null(modelName))
        {
            init <- modelName
        }
        else
        {
            init <- "Siena"
        }
        savefilename <- tclvalue(tkgetSaveFile(filetypes=modelFiletypes,
                                               defaultextension='.Rdata',
                                               initialfile=init))
        if (savefilename != "")
            save(mymodel, mydata, myeff, file=savefilename)
    }
    ##@sessionFromTcl internal siena01Gui
    sessionFromTcl <- function()
    {
        rows <- as.numeric(strsplit(tclvalue(tkconfigure(table1,  '-rows')),
                                    " ")[[1]][5])
        ##height <- as.numeric(strsplit(tclvalue(tkconfigure(table1,  '-height')),
        ##                              " ")[[1]][5])
        if (tableRows != (rows-1))
            browser()
        if (is.null(session))
        {
            session <<- data.frame(Group = 1, Name ="",
                                   Filename = "",
                                   Format = "Matrix",
                                   Period = "1",
                                   ActorSet = "Actors",
                                   Type = "network",
                                   Selected = "Yes",
                                   MissingValues = "NA",
                                   NonZeroCode = "1",
                                   NbrOfActors = "",
                                   stringsAsFactors = FALSE)

            session <<- session[rep(1, noFiles),]
            row.names(session) <<- 1:noFiles
        }
        for (i in 1:noFiles)
        {
            for (j in c(1,2,3,5,6,8,9,10,11))
            {
                if (is.null( mytkarray[[i,j]]) ||
                    tclvalue(mytkarray[[i,j]]) =="")
                {
                    mytkarray[[i,j]] <<- as.tclObj(defaults[j], drop=TRUE)
                }
                session[i, j] <<- trim.blanks(tclvalue(mytkarray[[i,j]]))
            }
            session[i, 4] <<- tclvalue(tkget(formatspins[[i]]))
            session[i, 7] <<- tclvalue(tkget(typespins[[i]]))
        }
        ##one day we will validate too!
    }

    ##@modelFromTcl internal siena01Gui
    modelFromTcl <- function()
    {
       # model <- NULL
        if (!is.null(modelName))
        {
            projname <- modelName
        }
        else
        {
            projname <- "Siena"
        }
        cond <- tclvalue(estimVar) ==
            '1. conditional Method of Moments'
        firstg <- as.numeric(tclvalue(gainVar))
        useStdInits <- tclvalue(stdstartVar) == '1'
        nsub <- as.numeric(tclvalue(ph2spinVar))
        if (tclvalue(rsVar) == '0')
        {
            seed <- NULL
        }
        else
        {
            seed <- as.numeric(tclvalue(rsspinVar))
        }
        FinDiff.method <- tclvalue(derivVar) == '0. crude Monte Carlo'
        n3 <- as.numeric(tclvalue(ph3spinVar))
        degs <- rep(0, nMaxDegree)
        for (i in 1:nMaxDegree)
        {
            degs[i] <- as.integer(tclvalue(maxdfVar[[i, 2]]))
        }
        names(degs) <- depvarnames[maxDegree]
        condvarno <- 0
        condname <- ""
        if (cond)
        {
            if (ndepvars == 1)
            {
               condvarno <- 1
               condname <- ""
            }
            else
            {
                condname <- tclvalue(condVar)
            }
        }
        sienaAlgorithmCreate(projname=projname, useStdInits=useStdInits,
                         cond=cond, firstg=firstg, seed=seed,
                         nsub=nsub, n3=n3, findiff=FinDiff.method,
                         MaxDegree=degs, condvarno=condvarno, condname=condname)
    }
    ##@sienaModelOptions internal siena01Gui
    sienaModelOptions <- function()
    {
        ##@editFn internal siena01Gui
        editFn <- function()
        {
            ## split effects if a variable is selected
            theseEffects <- tclvalue(effectsVar)
            myeffcopy <- myeff
            if (theseEffects != "")
            {
                myeffcopy <- myeff[myeff$name == theseEffects, ]
            }
            if (is.null(myeffcopy$effectNumber))
            {
                myeffcopy <- cbind(effectNumber=1:nrow(myeff), myeff,
                               effect1=rep(0, nrow(myeff)),
                               effect2=rep(0, nrow(myeff)),
                               effect3=rep(0,nrow(myeff)))
            }
             if (is.null(myeffcopy$timeDummy))
            {
                myeffcopy$timeDummy <- rep(",", nrow(myeff))
            }
            editCols <- c("name", "effectName", "type", "include", "fix",
                          "test", "initialValue", "parm", "effectNumber",
                          "effect1", "effect2", "effect3", "timeDummy")
            effEdit <- myeffcopy[, editCols]
            for (i in c("include", "fix", "test"))
            {
                effEdit[,i] <- as.numeric(effEdit[,i])
            }
            effEdit <- utils:::edit.data.frame(effEdit, edit.row.names=FALSE)
            for (i in c("include", "fix", "test"))
            {
                effEdit[,i] <- as.logical(effEdit[,i])
            }
            myeffcopy[, editCols] <- effEdit
            if (theseEffects != "")
            {
            myeff[myeff$name == theseEffects, ] <<- myeffcopy
            }
            else
            {
                myeff <<- myeffcopy
            }
            ##  browser()
            ## make sure this window is top with a global grab,
            ##but only for a second
            tcl('wm', 'attributes', tt, '-topmost', 1)
            Sys.sleep(0.1)
            tcl('wm', 'attributes', tt, '-topmost', 0)
                                        #      tkfocus(tt)
        }
        ##@estimateFn internal siena01Gui
        estimateFn <- function()
        {
            ##create mymodel
            mymodel <<- modelFromTcl()
            if (tclvalue(clustVar) == '0')
            {
                nbrNodes <- 1
            }
            else
            {
                nbrNodes <- as.numeric(tclvalue(clustspinVar))
            }
            if (nbrNodes > 1)
            {
                resp <- try(siena07(mymodel, data=mydata, effects=myeff,
                                    useCluster=TRUE, initC=TRUE,
                                    nbrNodes=nbrNodes),
                            silent=TRUE)
            }
            else
            {
                resp <- try(siena07(mymodel, data=mydata, effects=myeff),
                            silent=TRUE)
            }
            if (inherits(resp, "try-error"))
            {
                tkmessageBox(message=resp, icon="error")
            }
            else ## update the thetas to use next time, if run not interrupted
            {
                estimAns <<- resp
                if (estimAns$cconditional)
                {
                    ## z$condvar has the subscripts of included parameters that
                    ## correspond to the conditional variable
                    if (!is.null(estimAns$rate))
                    {
                        efflist <- apply(myeff[myeff$include, ], 1, function(x)
                                         paste(x[c("name", "shortName",
                                                   "type", "groupName",
                                                   "interaction1",
                                                   "interaction2", "period")],
                                               collapse="|"))
                        condeff <- attr(estimAns$f, "condEffects")
                        condeff$initialValue <- estimAns$rate
                        estimAns$effects$initialValue <- estimAns$theta
                        neweff <- rbind(estimAns$effects, condeff)

                        newlist <- apply(neweff, 1, function(x)
                                         paste(x[c("name", "shortName",
                                                   "type", "groupName",
                                                   "interaction1",
                                                   "interaction2", "period")],
                                               collapse="|"))
                        use <- which(myeff$include)
                        initValues <- rep(0, length(use))
                        initValues <- neweff$initialValue[match(efflist,
                                                                newlist)]
                        myeff$initialValue[myeff$include] <<- initValues
                    }
                }
                else
                {
                    if (!estimAns$termination == "UserInterrupt")
                    {
                        efflist <- apply(myeff[myeff$include, ], 1, function(x)
                                         paste(x[c("name", "shortName",
                                                   "type", "groupName",
                                                   "interaction1",
                                                   "interaction2", "period")],
                                               collapse="|"))
                        newlist <- apply(estimAns$effects, 1, function(x)
                                         paste(x[c("name", "shortName",
                                                   "type", "groupName",
                                                   "interaction1",
                                                   "interaction2", "period")],
                                               collapse="|"))
                        subs <- match(efflist, newlist)
                        myeff$initialValue[myeff$include] <<-
                            estimAns$theta[subs]
                    }
                }
                wasopen <- FALSE
                if (resultsOpen)
                {
                    resultsFn()
                    wasopen <- TRUE
                }
                resultsFn()
                tcl("wm", "attributes", tt, "-topmost", 1)
                Sys.sleep(0.1)
                tcl("wm", "attributes", tt, "-topmost", 0)
                tkfocus(tt)
                if (wasopen)
                {
                    tmp <- tksearch(pt, "-backwards", "New Analysis", "end")
                    tkyview(pt, "-pickplace", tmp)
                }
                tkfocus(resultsframe)
            }
        }
        ##@modelhelpFn internal siena01Gui
        modelhelpFn <- function()
        {
			RShowDoc("RSiena_Manual", package=pkgname)
        }
        ##@randomseedFn internal siena01Gui
        randomseedFn <- function()
        {
            val <- as.numeric(tclvalue(rsVar))
            if (val == 0)
            {
                tkgrid.forget(rsspin)
                tclvalue(rsspinVar) <<- 0
            }
            else
            {
                tkgrid(rsspin, row=3, column=1)
            }
        }
        ##@clustersFn internal siena01Gui
        clustersFn <- function()
        {
            val <- as.numeric(tclvalue(clustVar))
            if (val == 0)
            {
                tkgrid.forget(clustspin)
                tclvalue(clustspinVar) <<- 0
            }
            else
            {
                tkgrid(clustspin, row=4, column=1)
            }
        }
        ##@returnFn internal siena01Gui
        returnFn <- function()
        {
            ans <- tkmessageBox(message="Do you want to save the model?",
                                type="yesno", icon="question")
            if (tclvalue(ans)=="yes")
            {
                savefileFn()
            }
            ## null the objects because
            ##they may not match the session any more
            mydata <<- NULL
            myeff  <<-  NULL
            mymodel <<- NULL

            tkpack.forget(f1)
            tkpack.forget(resultsframe)
            tkpack(frame1, side='top', padx=5)
            tkpack(f0, side="bottom")
        }
        ##@showFn internal siena01Gui
        showFn <- function()
        {
            if (is.null(myeff$effectNumber))
            {
                myeff <- cbind(effectNumber=1:nrow(myeff), myeff,
                                   effect1=rep(0, nrow(myeff)),
                                   effect2=rep(0, nrow(myeff)),
                                   effect3=rep(0,nrow(myeff)))
            }
            if (is.null(myeff$timeDummy))
            {
                myeff$timeDummy <- rep(",", nrow(myeff))
            }
            editCols <- c("name", "effectName", "type", "include", "fix",
                          "test", "initialValue", "parm", "effectNumber",
                          "effect1", "effect2", "effect3", "timeDummy")
            effEdit <- myeff[myeff$include, editCols]
            for (i in c("include", "fix", "test"))
            {
                effEdit[,i] <- as.numeric(effEdit[,i])
            }
            utils:::edit.data.frame(effEdit, edit.row.names=FALSE)
            ##  tkfocus(tt)
            ## make sure this window is top with a global grab,
            ## but only for a second
            tcl('wm', 'attributes', tt, '-topmost', 1)
            Sys.sleep(0.1)
            tcl('wm', 'attributes', tt, '-topmost', 0)
        }
        ##@resultsFn internal siena01Gui
        resultsFn <- function()
        {
            if (resultsOpen)
            {
                resultsOpen <<- FALSE
                tclclose(resultsFileID)
                tkpack.forget(resultsframe)
            }
            else
            {
                resultsFile <- paste(projname, ".out", sep="")
                resultsFileID <<- tclopen(resultsFile)
                tkpack(resultsframe, side="bottom")
                tkconfigure(pt, state="normal")
                tkgrid(pt, yscr)
                tkgrid(xscr)
                tkgrid.configure(yscr,sticky="ns")
                tkgrid.configure(xscr,sticky="ew")
                tkdelete(pt, "1.0", "end")
                tkinsert(pt, "end", tclread(resultsFileID))
                tkconfigure(pt, state="disabled")
                tkfocus(resultsframe)
                resultsOpen <<- TRUE
            }
        }
        ##@saveresultsFn internal siena01Gui
        saveresultsFn <- function()
        {
            if (is.null(estimAns))
            {
                tkmessageBox(message="No results to save", icon="error")
                return()
            }
            mymodel <<- modelFromTcl()
            modelFiletypes <- "{{R object} *.Rdata}"
            if (!is.null(modelName))
            {
                init <- paste(modelName, "Results", sep="")
            }
            else
            {
                init <- "sienaResults"
            }
            savefilename <- tclvalue(tkgetSaveFile(filetypes=modelFiletypes,
                                                   defaultextension='.Rdata',
                                                   initialfile=init))
            if (savefilename != "")
                save(estimAns, file=savefilename)
        }
        ##@screenFromModel internal siena01Gui update screen from saved model
        screenFromModel <- function()
        {
            if (exists("mymodel"))
            {
                if (!is.na(mymodel$cconditional))
                {
                    if (mymodel$cconditional)
                    {
                        tclvalue(estimVar) <<-
                            '1. conditional Method of Moments'
                        if (ndepvars > 1)
                        {
                            tclvalue(condVar) <<- mymodel$condvarno
                        }
                    }
                    else
                    {
                        tclvalue(estimVar) <<-
                            '0. unconditional Method of Moments'
                    }
                }
                tclvalue(gainVar) <<- mymodel$firstg
                tclvalue(stdstartVar) <<- as.numeric(mymodel$useStdInits)
                tclvalue(ph2spinVar) <<- mymodel$nsub
                if (!is.null(mymodel$randomSeed) && mymodel$randomSeed != 0)
                {
                    tclvalue(rsVar) <<- '1'
                    tclvalue(rsspinVar) <<- mymodel$randomSeed
                    randomseedFn()
                }
                if (mymodel$FinDiff.method)
                {
                    tclvalue(derivVar) <<- '0. crude Monte Carlo'
               }
                tclvalue(ph3spinVar) <<- mymodel$n3
                if (any(mymodel$MaxDegree != 0))
                {
                    for (i in 1:nMaxDegree)
                    {
                        maxdfVar[[i, 2]] <<- mymodel$MaxDegree[depvarnames[i]]
                    }
               }
            }
        }
        ##*######################################
        ## start of sienaModelOptions function proper
        ##*#####################################
        resultsOpen <- FALSE
        if (inherits(mydata, "sienaGroup"))
        {
            depvarnames <<- attr(mydata, "netnames")
            ndepvars <<- length(depvarnames)
            nettypes <<- attr(mydata, "types")
        }
        else
        {
            depvarnames <<- names(mydata$depvars)
            ndepvars <<- length(depvarnames)
            nettypes <<- sapply(mydata$depvars, function(x)attr(x, "type"))
        }
        if (!is.null(modelName))
        {
            projname <- modelName
        }
        else
        {
            projname <- "Siena"
        }
        ## get rid of the previous windows
        tkpack.forget(frame1)
        tkpack.forget(f0)

        ## create and display the outer frame of model option screen
        f1 <- tkframe(tt, relief="ridge", width=1000, height=100)
        tkpack(f1)

        ## create and display the labelframe
        optiontt <- tkwidget(f1, 'labelframe',
                             text=paste(" Model Options: ", projname, " "),
                             width=1000, height=300)
        tkpack(optiontt, padx=30, pady=10)

        ## create the results frame to display later
        resultsframe <- tkframe(f1)
        xscr <- tkscrollbar(resultsframe, orient='horizontal',
                            command=function(...)tkxview(pt, ...))
        yscr <- tkscrollbar(resultsframe,
                            command=function(...)tkyview(pt, ...))
        pt <- tktext(resultsframe, width=90, bg="white", font="courier",
                     wrap='none', height=20,
                     xscrollcommand=function(...)tkset(xscr, ...),
                     yscrollcommand=function(...)tkset(yscr, ...))

        ## create and display top inner frame
        optf <- tkwidget(optiontt, 'frame', borderwidth=2, relief="groove")
        tkgrid.configure(optf, padx=5, pady=5)#), columnspan=2)
        maxdf <- tkframe(optf, borderwidth=2, relief='groove')
        tkgrid.configure(maxdf, padx=5, pady=5, rowspan=4, column=4, row=0)

        ## create and display estimation method option box and its label
        estimlist <- c('0. unconditional Method of Moments',
                       '1. conditional Method of Moments')
        estimVar <<- tclVar()
        if (ndepvars == 1)
        {
            tclvalue(estimVar) <<- '1. conditional Method of Moments'
        }
        else
        {
            tclvalue(estimVar) <<- '0. unconditional Method of Moments'
        }
        estim <- ttkcombobox(optf, values=estimlist, state="readonly",
                             textvariable=estimVar, width=30)
        estimlab <- tklabel(optf, text='Estim. method')
        tkgrid(estimlab, estim)

        tkgrid.configure(estimlab, column=0, row=0)
        tkgrid.configure(estim, column=1, row=0)

        ## create options box to select a conditional variable
        if (ndepvars > 1)
        {
            condVar <<- tclVar()
            condvarl <- ttkcombobox(optf, values= depvarnames, state='readonly',
                                    textvariable=condVar, width=10)
            condlab <- tklabel(optf, text=' Conditioning Variable ')
            tkgrid(condlab,  row=1, column=0)
            tkgrid( condvarl, row=1, column=1)
            tclvalue(condVar) <<- depvarnames[1]
        }
        ## create and display initial value of gain parameter box
        gainlab <- tklabel(optf, text=' Initial value of gain parameter ')
        gainVar <<- tclVar()
        tclvalue(gainVar) <<- '0.2'
        gain <- tkentry(optf, width=10, textvariable=gainVar, cursor="arrow")
        tkgrid(gainlab, gain, padx=5)

        ## tidy up the screen
        tkgrid.configure(estimlab, sticky='w')
        tkgrid.configure(gainlab, row=0, column=2, sticky='w')
        tkgrid.configure(gain, row=0, column=3, sticky='w')

        ## create and display box for standard starting values
        stdstartVar <<- tclVar()
        tclvalue(stdstartVar) <<- '0'
        stdstart <- tkcheckbutton(optf, text=' Standard starting value ',
                                  variable=stdstartVar)
        tkgrid(stdstart, row=2, sticky='w', columnspan=2)

        ## create and display box for number of phase 2 subphases
        ph2lab <- tklabel(optf, text=' Number of phase 2 subphases ')
        ph2spinVar <<- tclVar()
        tclvalue(ph2spinVar) <<- 4
        ph2spin <-  tkwidget(optf, 'spinbox', from=0, to=10, width=10,
                             textvariable=ph2spinVar, cursor="arrow")
        tkgrid(ph2lab, row=1, column=2, padx=5)
        tkgrid(ph2spin, row=1, column=3, sticky='w', padx=5)

        ##create and display fields for random number entry
        rsVar <<- tclVar()
        tclvalue(rsVar) <<- '0'
        rs <- tkcheckbutton(optf, text=' Specify random seed: ', variable=rsVar,
                            command=randomseedFn)
        rsspinVar <<- tclVar()
        rsspin <-  tkwidget(optf, 'spinbox', from=0, to=1000000, width=10,
                            textvariable=rsspinVar, cursor="arrow")
        tkgrid(rs, row=3, sticky='w', columnspan=2)

        ##create and display fields for number of processors entry
        clustVar <<- tclVar()
        tclvalue(clustVar) <<- '0'
        clust <- tkcheckbutton(optf,
                               text=' Multiple processors: ',
                               variable=clustVar,
                               command=clustersFn)
        clustspinVar <<- tclVar()
        clustspin <-  tkwidget(optf, 'spinbox', from=2, to=1000, width=10,
                            textvariable=clustspinVar, cursor="arrow")
        tkgrid(clust, row=4, sticky='w', columnspan=2)

        ##create and display field for derivative method
        derivlab <- tklabel(optf, text=' Derivative method ')
        derivlist <- c('0. crude Monte Carlo',
                       '1. score function')
        derivVar <<- tclVar()
        tclvalue(derivVar) <<- '1. score function'
        derivw <- ttkcombobox(optf, values=derivlist, state="readonly",
                              textvariable=derivVar, width=18)
        tkgrid(derivlab,  row=2, column=2, sticky='w', padx=5)
        tkgrid(derivw,  row=2, column=3, sticky='w', padx=5)

        ##create and display field for number of phase 3 iterations
        ph3lab <- tklabel(optf, text=' Number of phase 3 iterations ')
        ph3spinVar <<- tclVar()
        tclvalue(ph3spinVar) <<- 1000
        ph3spin <-  tkwidget(optf, 'spinbox', from=10, to=1000000, width=10,
                             textvariable=ph3spinVar, cursor="arrow")
        tkgrid(ph3lab, row=3, column = 2, sticky='w', padx=5)
        tkgrid(ph3spin, row=3, column=3, sticky='w', padx=5)

        ##create and display field for restricting degree of model
        maxdfVar <<- tclArray()
        xscr2 <- tkscrollbar(maxdf, orient="horizontal",
                             command=function(...)tkxview(table2,...))
        yscr2 <- tkscrollbar(maxdf, command=function(...)tkyview(table2,...))
        ## find out how many are networks
        maxDegree <<- nettypes != "behavior"
        nMaxDegree <<- sum(maxDegree)
        table2 <- tkwidget(maxdf, 'table', cols=2, rows=ndepvars+1,
                           height=min(10, nMaxDegree + 1), variable = maxdfVar,
                           titlerows=1, roworigin=0, cursor="arrow",
                           colorigin=1, anchor="w", background='white',
                           xscrollcommand=function(...) tkset(xscr2,...),
                           yscrollcommand=function(...) tkset(yscr2,...),
                           multiline=0, rowheight=-25, font='Helvetica 10',
                           resizeborders='col', maxwidth=300,
                           selectmode='single', sparsearray=0)
        if(.Platform$OS.type == "windows") ##  these colours only exist there
        {
            tcl(table2, 'tag','configure','active',fg='SystemHighlightText',
                bg='SystemHighlight')
            tcl(table2, 'tag','configure','title',fg='SystemHighlightText',
                bg='SystemHighlight')
        }
        maxdfVar[[0,1]] <<- as.tclObj("NetWork", drop=TRUE)
        maxdfVar[[0,2]] <<- as.tclObj('Max Degree', drop=TRUE)
        for (i in 1:nMaxDegree)
        {
            maxdfVar[[i, 1]] <<- depvarnames[maxDegree][i]
            maxdfVar[[i, 2]] <<- 0
        }
        tkpack(yscr2, fill="y", side="right")
        tkpack(xscr2, fill="x", side="bottom")
        tkpack(table2)
        tcl(table2, 'width', 1, 15)
        tcl(table2, 'width', 2, 10)

        ## create and display the frame for the buttons.
        comf <- tkframe(optiontt,  borderwidth=2, relief='groove')
        tkgrid(comf, row=1, column=0)

        ## create the buttons
        ## create options box to select an effects list
        ## if (ndepvars > 1)
        ## {
        effectsVar <<- tclVar()
        effectsvarl <- ttkcombobox(comf,
                                   values= depvarnames, state='readonly',
                                   textvariable=effectsVar, width=10)
        effectslab <- tklabel(comf, text=' Effects dependent variable ')
        ##}
        editbut <- tkbutton(comf, command=editFn,
                            text=' Edit effects (selected variable) ',
                            width=27)
        showbut <- tkbutton(comf, command=showFn,
                            text=' Show included effects (all) ', width=22)
        applybut <- tkbutton(comf, command=estimateFn, text=' Estimate ',
                             width=22)
        saveresultsbut <- tkbutton(comf, command=saveresultsFn,
                                   text=' Save results ', width=22)
        savefilebut <- tkbutton(comf, command=savefileFn,
                                text=' Save to file ', width=22)
        resultsbut <- tkbutton(comf, command=resultsFn,
                               text=' Display Results ', width=22)
        returnbut <- tkbutton(comf, command=returnFn,
                              text=' Exit Model Options ', width=22)
        helpbut <- tkbutton(comf, command=modelhelpFn, text=' Help ',
                            width=22)
        ## if (ndepvars > 1)
        ## {
        tkgrid(effectslab, effectsvarl,
               applybut,  savefilebut, returnbut, padx=5, pady=5)
        ##}
        ##else
        ## {
        ##   tkgrid(editbut, showbut, applybut,  savefilebut, padx=5, pady=5)
        ## }
        tkgrid(editbut, showbut,  resultsbut, saveresultsbut,
                helpbut, padx=5, pady=5)
        tkgrid(effectsvarl, sticky="w")
        screenFromModel()
       ## make sure this window is top with a global grab, bu only for a second
        tcl('wm', 'attributes', tt, '-topmost', 1)
        Sys.sleep(0.1)
        tcl('wm', 'attributes', tt, '-topmost', 0)
    }
    ## ##############################################
    ## start of siena01Gui function proper
    ## ##############################################
    if (getDocumentation)
    {
        return(getInternals())
    }
    ## check we have the right libraries
    library(tcltk)
    if (!inherits(tclRequire("Tktable"), "tclObj"))
        stop("This function needs the tcl/tk package TkTable: install it, ",
             "or use an alternative data entry method: see RSiena help page")
    ## directory for startup
    initialDir <- getwd()

    ## create toplevel window
    tt <- tktoplevel(class="mytoplevel")
    ## configure it
    tkwm.resizable(tt, FALSE, FALSE)
    tkwm.title(tt, "Siena01")
    tkbind("mytoplevel", "<Destroy>", myStop)

    ## create menu
    topMenu <- tkmenu(tt)
    tkconfigure(tt, menu=topMenu)
    fileMenu <- tkmenu( topMenu, tearoff=FALSE)
    tkadd(fileMenu, "command", label = "Quit",
          command= function() tkdestroy(tt))
    tkadd(topMenu, "cascade", label = "File",
          menu = fileMenu)

    ## create and display topframe
    frame1 <- tkframe(tt,height=300,relief='ridge')
    tkpack(frame1,side='top',padx=5)

    ## create and display buttons on top frame
    fromFile <- tkbutton(frame1, command=fromFileFn, width=20,
                         text=' Load new session from file ')
    fromFileCont <- tkbutton(frame1, command=fromFileContFn, width=20,
                             text=' Continue session from file ')
    tkgrid( fromFile, fromFileCont, padx=50, pady=10)

    ## create and display bottom frame
    f0 <- tkframe(tt, relief='ridge',height=400)
    tkpack(f0, side="bottom")

    ## create and display labelframe Data Definition
    df <- tkwidget(f0, 'labelframe', text='Data Definition',
                   width=2000, height=300)
    tkgrid(df, columnspan=4)

    ## create and display bottom buttons
    save1 <- tkbutton(f0, command = saveFn, width=10, text=' Save to file ' )
    apply1 <- tkbutton(f0, command = applyFn, width=10, text=' Apply ')
    clear1 <- tkbutton(f0, command = clearFn, width=10, text=' Clear ' )
    help1 <- tkbutton(f0, command = helpFn, width=10, text=' Help ')
    tkgrid(save1, apply1,  clear1, help1, pady=10)

    ## create and display session panel
    f1 <- tkframe(df,  height=300)
    tkpack(f1)
    mytkarray <- tclArray() ## variable for the table
    xscr <- tkscrollbar(f1, orient="horizontal",
                        command=function(...)tkxview(table1, ...))
    yscr <- tkscrollbar(f1, command=function(...)tkyview(table1, ...))
    tableRows <- 2
    table1 <- tkwidget(f1, 'table', cols=11, rows=tableRows+1, height=3,
                       variable= mytkarray, titlerows=1, roworigin=0,
                       colorigin=1, anchor="w", cursor="arrow",
                       xscrollcommand=function(...) tkset(xscr,...),
                       yscrollcommand=function(...) tkset(yscr,...),
                       multiline=0, rowheight=-25, maxwidth=900,
                       resizeborders="col",
                       selectmode='single', sparsearray=0)
    ## set the colours on Windows, otherwise defaults
    if (.Platform$OS.type == 'windows')
    {
        tkconfigure(table1, background='white')
        ## final configuration colour uses the active tag, rather than selected
        ## so need to set this although don't use it elsewhere!
        tcl(table1, 'tag','configure','active',fg='SystemHighlightText',
            bg='SystemHighlight')
        tcl(table1, 'tag','configure','title',fg='SystemHighlightText',
            bg='SystemHighlight')
    }

    ## set up headings for the table
    mytkarray[[0,1]] <- 'Group'
    mytkarray[[0,2]] <- 'Name'
    mytkarray[[0,3]] <- 'Filename'
    mytkarray[[0,4]] <- 'Format'
    mytkarray[[0,5]] <- 'Period(s)'
    mytkarray[[0,6]] <- 'ActorSet'
    mytkarray[[0,7]] <- 'Type'
    mytkarray[[0,8]] <- 'Selected'
    mytkarray[[0,9]] <- as.tclObj("MissingValues", drop=TRUE)
    mytkarray[[0,10]] <- as.tclObj("NonZeroCode", drop=TRUE)
    mytkarray[[0,11]] <- "NbrOfActors"

    ##create spinboxes for format
    ff <- c("matrix {pajek net} {Siena net}")
    formatspins <- lapply(1:tableRows, function(x)
                          tkwidget(table1, 'spinbox', state='readonly',
                                   width=20, values=ff, cursor="arrow"))
    ## insert them in the table, and set focus bindings?
    lapply(1:tableRows, function(x)
       {
           mypos <- paste(x, ',',4, sep='')
           tkwindow.configure(table1, mypos, window=formatspins[[x]])
           tkbind(formatspins[[x]], "<FocusIn>",
                  function(y){tcl(table1,"activate", "1,4")
                              tkXselection.own(selection=y)})
       }
           )

    ##create spinboxes for type
    typelist <- c("network", "bipartite", "behavior", "constant covariate",
                  "changing covariate", "constant dyadic covariate",
                  "changing dyadic covariate", "exogenous event")
    typespins <- lapply(1:tableRows, function(x)
                        tkwidget(table1, 'spinbox', state='readonly',
                                 width=25, values=typelist, cursor="arrow"))
    lapply(1:tableRows, function(x) {
        mypos <- paste(x, ',',7, sep='')
        tkwindow.configure(table1, mypos, window=typespins[[x]])
        tkbind(typespins[[x]], "<FocusIn>",
               function(y){ tcl(table1,"activate", "1,4")
                            tkXselection.own(selection=y)})
    })

    ## display the table
    tkpack(yscr, fill='y', side='right')
    tkpack(xscr, fill="x", side='bottom')
    tkpack(table1)

    ## configure the table column widths
    tcl(table1, 'width', 1, 6)
    tcl(table1, 'width', 3, 30)
    tcl(table1, 'width', 5, 8)
    tcl(table1, 'width', 6, 8)
    tcl(table1, 'width', 7, 24)
    tcl(table1, 'width', 8, 8)
    tcl(table1, 'width', 9, 12)
    tcl(table1, 'width', 10, 12)
    tcl(table1, 'width', 11, 10)

    ## create and display the bottom bottons
    f2 <- tkframe(df)
    tkpack(f2, side='bottom')

    add1 <- tkbutton(f2, command = addFile, width=10, text=' Add ')
    rem1 <- tkbutton(f2, command = removeFile, width=10, text=' Remove ' )
    edit1 <- tkbutton(f2, command = editFile, width=10, text=' Edit ')

    tkgrid(add1,  rem1, edit1, pady=10)

    ## to help in debugging, make a copy of the relevant function
    ## then use eg environment(applyFn)$modelName to examine the values!
    ## applyFn <<- applyFn
    ## alternatively, put a browser in a call back, then use
    ## ls(parent.env(environment())) to find the variables

    ## sort out focus
    tkfocus(tt)
    ## put on top globally temporarily
    tcl('wm', 'attributes', tt, '-topmost', 1)
    Sys.sleep(0.1)
    tcl('wm', 'attributes', tt, '-topmost', 0)
    invisible()
}

