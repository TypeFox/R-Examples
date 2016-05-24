readTable  =
function ()
{
	#temp <- ""
    #if (exists("readClipboard")) readClipboard()->temp
   temp<- try(scan("clipboard", what="character", sep="\n", quiet=TRUE), silent=TRUE)

    tf <- tktoplevel()
    tkwm.title(tf, "Choose a file")
   # size=c(320,420,0,132)
   # geo <- paste(size[1], "x", size[2], "+", size[3],"+", size[4], sep = "")
   # tkwm.geometry(tf, geo)
    tkwm.resizable(tf,0,0) 
    tkwm.geometry(tf, paste("+0+",.heigth, sep=""))
    tkwm.deiconify(tf)
    tkgrab.set(tf)
    tkfocus(tf) 

    done <- tclVar(0)
    filenamevar <- tclVar("clipboard")
    tabnamevar <- tclVar("")
    fictrt <- tclVar()
    varnames <- tclVar(1)
    rowNames<- tclVar(1)#######################
    show.data <- tclVar(0)
    sepVar <- tclVar(1)
    otherSepVar <- tclVar("")
    decSepVar <- tclVar(1)
    colClassVar <- tclVar(4)
    
    choosefile <- function() {
        fictrt <- tkgetOpenFile()
        fpath <- tclvalue(fictrt)
         tkfocus(tf)
        if(fpath!=""){
        dat = file(fpath, "r")
         dat1 = readLines(dat, warn=FALSE)
           close(dat)
           
         tkconfigure(TXT, state="normal")
         tkdelete(TXT, "0.0", "end")
         for (i in 1:length(dat1)){        
         tkinsert(TXT,"end",paste(dat1[i]," ", "\n", sep=""))}
       
         tkconfigure(TXT, state="disabled")
            
        tkdelete(file.entry, 0, "end")
        tkinsert(file.entry, "end", fpath)
        tkdelete(tab.entry, 0, "end")
        tkinsert(tab.entry, "end", basename(fpath))
        }
     }
     
    choosefic <- function() {
        if (tclvalue(tabnamevar) != "") {
            tabname <- parse(text = tclvalue(tabnamevar))[[1]]
        }
        else tabname <- "untitled"
        if (tclvalue(filenamevar) != "") {
            filename <- tclvalue(filenamevar)
        }
        else return()
        varn <- as.logical(tclObj(varnames))
        sep <- tclvalue(sepVar)
        if (sep == 1) sepch <- ""
        if (sep == 2) sepch <- ","
        if (sep == 3) sepch <- ";"
        if (sep == 4) sepch <- "\t"
        if (sep == 5) {
            if (tclvalue(otherSepVar) != "") {
                otherSep <- tclvalue(otherSepVar)
            }
            else otherSep <- ""
            sepch <- otherSep
        }
        decSep <- tclvalue(decSepVar)
        if (decSep == 1)
            decsepch <- "."
        if (decSep == 2)
            decsepch <- ","
            rdcom <- paste(tabname, " <<- read.table(file='",
                filename, "', header=", varn, ", sep='", sepch,
                "', dec='", decsepch, "')", sep = "")

        try(eval(parse(text = rdcom)), silent=TRUE)-> flagError
        if (class(flagError) == "try-error") {return(invisible(tk_messageBox(type = "ok", "Something went wrong!!!", caption = "")))}
        tkdestroy(tf)
        
        rowNames.flag <- as.logical(tclObj(rowNames))
        
        if (rowNames.flag) {
        eval(parse(text=paste("rownames(",tabname,")<<-", tabname,"[,1]", sep="")))
        eval(parse(text=paste(tabname,"<<-", tabname,"[,-1]", sep="")))
        }
        
        show.data.flag <- as.logical(tclObj(show.data))
         
    if (show.data.flag)    eval(parse(text = paste("edit(", tabname, ")", sep = "")))
    }
    
    frame1.a <- tkframe(tf, relief = "groove")
    frame1 <- tkframe(tf, relief = "groove")
    frame.preview <- tkframe(tf, relief = "groove")
    
    tkgrid(tklabel(frame1.a, text = "File options:", foreground = "blue"))
    tkpack(frame1.a, fill = "x")
      
    tab.entry <- tkentry(frame1, textvariable = tabnamevar)
    file.entry <- tkentry(frame1, textvariable = filenamevar)
    separator<-tklabel(frame1,text="")
    choosefile.but <- tkbutton(frame1, text = "...", command = function() choosefile())
    tkgrid(tklabel(frame1, text = "Select a file to read: "), file.entry, separator, choosefile.but,  sticky = "w")
    tkgrid(tklabel(frame1, text = "Enter name for data set: "), tab.entry,  sticky = "w")
    tkpack(frame1, fill = "x")
    
    
    
    src <- tkscrollbar(frame.preview, repeatinterval=5, command=function(...)tkyview(TXT,...))
    src1 <- tkscrollbar(frame.preview, repeatinterval=5,orient="horizontal", command=function(...)tkxview(TXT,...))
TXT<-tktext(frame.preview,height=5,width= 37,  yscrollcommand=function(...)tkset(src,...),xscrollcommand=function(...)tkset(src1,...),background="grey", wrap="none")
tkgrid(tklabel(frame.preview,text="Preview:", foreground = "blue"),  sticky = "w")
tkgrid(TXT,src,sticky="w")
tkgrid(src1,sticky="w")
tkgrid.configure(src,sticky="ns")
tkgrid.configure(src1,sticky="ew")
for (i in 1:NROW(temp)){tkinsert(TXT,"end",paste(temp[i],"\n", sep=""))}
tkconfigure(TXT, state="disabled")
tkpack(frame.preview, fill="x")
    
    frame2 <- tkframe(tf, relief = "groove")
    varnames.cbut <- tkcheckbutton(frame2, text = "Variables names on the first row of data file", variable = varnames)
    tkgrid(varnames.cbut, columnspan = 3, sticky = "w")
    
    rownames.cbut <- tkcheckbutton(frame2, text = "Years on the first column of data file", variable = rowNames)
    tkgrid(rownames.cbut, columnspan = 3, sticky = "w")
    
    show.data.cbut <- tkcheckbutton(frame2, text = "Show data frame", variable = show.data)
    tkgrid(show.data.cbut, sticky="w")
    sepFrame <- tkframe(frame2, relief = "groove")
    sep.entry <- tkentry(sepFrame, textvariable = otherSepVar, width = 5)
    tkgrid(tklabel(sepFrame, text = "Field separator:", foreground = "blue"))
    tkgrid(tkradiobutton(sepFrame, text = "Default", value = 1, variable = sepVar), sticky = "w")
    tkgrid(tkradiobutton(sepFrame, text = "Commas", value = 2,  variable = sepVar), sticky = "w")
    tkgrid(tkradiobutton(sepFrame, text = "Semicolon", value = 3, variable = sepVar), sticky = "w")
    tkgrid(tkradiobutton(sepFrame, text = "Tab", value = 4, variable = sepVar), sticky = "w")
    tkgrid(tkradiobutton(sepFrame, text = "Other", value = 5, variable = sepVar), sep.entry, sticky = "w")
    decSepFrame <- tkframe(frame2, relief = "groove" )
    tkgrid(tklabel(decSepFrame, text = "Decimal separator:", foreground = "blue"))
    tkgrid(tkradiobutton(decSepFrame, text = "Period [.]", value = 1,  variable = decSepVar), sticky = "w")
    tkgrid(tkradiobutton(decSepFrame, text = "Comma [,]", value = 2,   variable = decSepVar), sticky = "w")
    tkgrid(sepFrame,  sticky = "w")
    
    frame.exit<-  tkframe(tf, relief = "groove")
    fr.exit.space         <-   tklabel(frame.exit, text = " ")
    ok.but <- tkbutton(frame.exit, text =     "      Ok      ", command = function() choosefic())
    cancel.but <- tkbutton(frame.exit, text = "    Cancel    ", command = function() tkdestroy(tf))
    tkgrid( cancel.but,fr.exit.space, ok.but)
    tkgrid(decSepFrame,  frame.exit,  sticky = "w")
    tkpack(frame2, fill="x")
    tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tf, "<KeyPress-Return>", function() choosefic())
    tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
    tkwait.variable(done)
    if (tclvalue(done) == "2")
        return(0)
    tkdestroy(tf)
}
#  readTable()->a
