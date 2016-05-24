############################################
###########################################
sucmatrixDbox<- function()
{
#    Library("abind")
    env <- environment()
    initializeDialog(title=gettextRcmdr("Successional Model"))
    dsname <- tclVar("Do_Not_Save")
    ## incluido
	entryDsname <- tkentry(top, width="20", textvariable=dsname)
	tmaxVar <- tclVar("20")
	tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
	clVar <- tclVar("20")
	lnVar <- tclVar("20")
	clEntry <- tkentry(top, width = "4", textvariable = clVar)
	lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
	######
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    ##
    s.outerTableFrame <- tkframe(top)
    assign("s.tableFrame", tkframe(s.outerTableFrame), envir=env)
    ##
##################################
   setUpTable <- function(...)
   {
        tkdestroy(get(".tableFrame", envir=env))
        assign(".tableFrame", tkframe(outerTableFrame), envir=env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(rowsValue))
        make.col.names <- "labelRcmdr(.tableFrame, text='')"
        for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep="")
            assign(col.varname, tclVar(paste("st",j,"(t)", sep="")), envir=env) #name show at coluns
            make.col.names <- paste(make.col.names, ", ", "ttkentry(.tableFrame, width='5', textvariable=",  col.varname, ")", sep="")
            }
        eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
        for (i in 1:nrows)
        {
            varname <- paste(".tab.", i, ".1", sep="")
            assign(varname, tclVar("") , envir=env)
            row.varname <- paste(".rowname.", i, sep="")
            assign(row.varname, tclVar(paste("st",i,"(t+1)", sep="")), envir=env) ## row names show (first table) 
            make.row <- paste("ttkentry(.tableFrame, width='7', textvariable=",
                row.varname, ")", sep="")
            make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                varname, ")", sep="")
            for (j in 2:ncols){
                varname <- paste(".tab.", i, ".", j, sep="")
                assign(varname, tclVar(""), envir=env)
                make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=", varname, ")", sep="")
                }
            eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
          }
        tkgrid(get(".tableFrame", envir=env), sticky="w")
############## entry vector: time zero
        tkdestroy(get("s.tableFrame", envir=env))
        assign("s.tableFrame", tkframe(s.outerTableFrame), envir=env)
        s.nrows <- as.numeric(1)
        s.ncols <- as.numeric(tclvalue(rowsValue))
        s.make.col.names <- "labelRcmdr(s.tableFrame, text='')"
        for (j in 1:s.ncols) 
        {
            s.col.varname <- paste(".scolname.", j, sep="")
            assign(s.col.varname, tclVar(paste("st",j, sep="")), envir=env)
            s.make.col.names <- paste(s.make.col.names, ", ", "ttkentry(s.tableFrame, width='5', textvariable=", s.col.varname, ")", sep="")
        }
        eval(parse(text=paste("tkgrid(", s.make.col.names, ")", sep="")), envir=env)
        for (i in 1:s.nrows)
        {
	    s.varname <- paste(".stab.", i, ".1", sep="") 
            assign(s.varname, tclVar("") , envir=env)
	    s.row.varname <- paste(".srowname.", i, sep="")
            assign(s.row.varname, tclVar("prop"), envir=env)
######################################################## names at row sec.table
            s.make.row <- paste("ttkentry(s.tableFrame, width='7', textvariable=",s.row.varname, ")", sep="")
	    s.make.row <- paste(s.make.row, ", ", "ttkentry(s.tableFrame, width='5', textvariable=", s.varname, ")", sep="")
            for (j in 2:ncols)
            {
                s.varname <- paste(".stab.1.", j, sep="")
                assign(s.varname, tclVar(""), envir=env)
                s.make.row <- paste(s.make.row, ", ", "ttkentry(s.tableFrame, width='5', textvariable=",
                s.varname, ")", sep="")
            }
        eval(parse(text=paste("tkgrid(", s.make.row, ")", sep="")), envir=env)
        }
        tkgrid(get("s.tableFrame", envir=env), sticky="w")
   }
#####################
    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("2")
    rowsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=rowsValue,resolution=1, orient="horizontal", command=setUpTable)
    rowsShow <- labelRcmdr(rowColFrame, textvariable=rowsValue, width=2, justify="right")
    #colsValue <- tclVar("2")
    #colsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=colsValue,
 #       resolution=1, orient="horizontal", command=setUpTable)
    #colsShow <- labelRcmdr(rowColFrame, textvariable=colsValue, width=2, justify="right")
    onOK <- function()
    {
        tmax <- round(as.numeric(tclvalue(tmaxVar)))
        if (is.na(tmax) || tmax <= 0) 
        {
        errorCondition(message = "Number of simulations must be a positive integer")
        return()
        }
		cl <- round(as.numeric(tclvalue(clVar)))
        if (is.na(cl) || cl <= 0) 
        {
        errorCondition(message = "Number of coluns on the simulated arena must be a positive integer.")
        return()
        }
		ln <- round(as.numeric(tclvalue(lnVar)))
        if (is.na(ln) || ln <= 0) 
        {
        errorCondition("Number of lines on the simulated arena must be a positive integer.")
        return()
        }
		nrows <- as.numeric(tclvalue(rowsValue))
      ncols <- as.numeric(tclvalue(rowsValue))
      cell <- 0
      s.cell<-0
      counts <- rep(0, nrows^2)
      s.counts<-rep(0, ncols)
      row.names <- rep("", nrows)
      col.names <- rep("", nrows)
      s.col.names <- rep("", nrows)
#### transition matrix
        for (i in 1:nrows) row.names[i] <-eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
        for (j in 1:ncols) col.names[j] <-eval(parse(text=paste("tclvalue(", paste(".colname.", j, sep=""),")", sep="")))
        for (i in 1:nrows)
        {
            for (j in 1:ncols)
            {
                cell <- cell+1
                varname <- paste(".tab.", i, ".", j, sep="")
                counts[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep="")))) ## aqui ele guarda os valores das celulas
             }
         }
        counts <- na.omit(counts)
        t.counts<-matrix(counts, nrows, nrows, byrow=TRUE)
        if ((sum(counts>1)+ sum(counts<0))>0 )
        {
            errorCondition(recall=sucmatrixDbox, message=sprintf(gettextRcmdr("Transitions probabilities must be between 0 to 1 ")))
            return()
        }
#		t.counts<-matrix(counts, nrows, nrows, byrow=TRUE)
	sum.col<-apply(t.counts, 2, sum)
		if (sum(sum.col==1) != nrows)
		{
		errorCondition(recall=sucmatrixDbox, message=sprintf(gettextRcmdr("Transitions for each stage at time t (columns) must sum 1. ADJUSTED BY TOTALS")))
		t.counts=t.counts/sum.col     
		}        
      	if (length(unique(row.names)) != nrows)
      	{
            errorCondition(recall=sucmatrixDbox, message=gettextRcmdr("Row names are not unique."))
            return()
          }
        if (length(unique(col.names)) != ncols){
            errorCondition(recall=sucmatrixDbox, message=gettextRcmdr("Column names are not unique."))
            return()
           }
########################    
######## Time 0 vector
########################
        #s.row.names[1] <-eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
        for (j in 1:ncols)
        { 
        s.col.names[j] <-eval(parse(text=paste("tclvalue(", paste(".scolname.", j, sep=""),")", sep="")))
        s.cell <- s.cell+1
        s.varname <- paste(".stab.1.", j, sep="")
        s.counts[s.cell] <- as.numeric(eval(parse(text=paste("tclvalue(", s.varname,")", sep="")))) ## aqui ele guarda os valores das celulas
        }
        s.counts <- na.omit(s.counts)
		if (sum(s.counts)!=1)
		{
		errorCondition(message=sprintf(gettextRcmdr("Proportion sum for all stage at initial time (columns) must sum 1\n VALUES ADJUSTED BY TOTAL")))
		s.counts=s.counts/sum(s.counts)       
		}        
      if (length(unique(s.col.names)) != ncols)
      {
         errorCondition(recall=sucmatrixDbox, message=gettextRcmdr("Initial proportions row names are not unique."))
       return()
       }
###########################
		#sn.counts=s.counts*cl*ln      
      closeDialog()
### transition matrix
	.Table <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", nrows,", byrow=TRUE)", sep="")
#         assign(".Table", justDoIt(command), envir=.GlobalEnv)
#         logger(paste(".Table <- ", command, sep=""))
#	rownames(.Table)<- row.names
#         command <- paste("c(",paste(paste("'", row.names, "'", sep=""), collapse=", "), ")", sep="")
#         justDoIt(paste("rownames(.Table) <- ", command, sep=""))
#         logger(paste("rownames(.Table) <- ", command, sep=""))
#	 colnames(.Table)<-  col.names
#         command <- paste("c(",paste(paste("'", col.names, "'", sep=""), collapse=", "), ")", sep="")
#         justDoIt(paste("colnames(.Table) <- ", command, sep=""))
#         logger(paste("colnames(.Table) <- ", command, sep=""))
#         doItAndPrint(".Table  # Transition Probalities")
#################################
        .Nt <- paste("c(", paste(s.counts, collapse=","), ")", sep="")
#         assign(".Nt", justDoIt(command), envir=.GlobalEnv)
#         logger(paste(".Nt <- ", command, sep=""))
#        names(.Nt)<-  s.col.names
#         command <- paste("c(",paste(paste("'", s.col.names, "'", sep=""), collapse=", "), ")", sep="")
#         justDoIt(paste("names(.Nt) <- ", command, sep=""))
#         logger(paste("names(.Nt) <- ", command, sep=""))
#         doItAndPrint(".Nt  # Initial proportion of patchs per stage")
################################# 
##sucMatrix(mat.trans, init.prop, ln, cl, tmax)
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("sucMatrix(mat.trans=", .Table ,", init.prop  =", .Nt,", tmax = ",tmax, ", rw = ", ln,", cl = ", cl,")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<- sucMatrix(mat.trans=", .Table ,", init.prop  =", .Nt,", tmax = ",tmax, ", rw = ", ln,", cl = ", cl,")", sep = "")
		  }
########
		doItAndPrint(command)
	#logger("remove(.Table)")
        #remove(.Table, envir=.GlobalEnv)
        #logger("remove(.Nt)")
        #remove(.Nt, envir=.GlobalEnv)
      tkfocus(CommanderWindow())
      }
    OKCancelHelp(helpSubject="sucMatrix")
##############
tkgrid(tklabel(top, text="Enter name for simulation data set :", fg="blue"), sticky="w")
#    tkgrid(tklabel(top, text="\tEnter name for data set:"), entryDsname, sticky="e")
tkgrid(entryDsname,sticky="e" )
## incluido
tkgrid(tklabel(top, text="Simulation Arena Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns"), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows"), lnEntry, sticky = "e")
#tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
####    
tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter Transitions Probabilities: "), fg="blue"), sticky="w")
    tkgrid(labelRcmdr(rowColFrame, text=gettextRcmdr("Number of stages:")), rowsSlider, rowsShow, sticky="w")
   tkgrid(rowColFrame, sticky="w")

    tkgrid(labelRcmdr(top, text=gettextRcmdr("Columns: stages at time t"), fg="red"), sticky="e")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Rows: stages at \n\ttime t+1"), fg="red"), sticky="w")
    tkgrid(outerTableFrame, sticky="e")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Initial stages proportions: "), fg="blue"), sticky="w")
    tkgrid(s.outerTableFrame, sticky="e")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=10, columns=2)
}

##############################################
##########################################
regnichoDbox<-function () 
{
initializeDialog(title = gettextRcmdr("Niche Regeneration"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
####
tmaxVar <- tclVar("100")
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
######
c1Var <- tclVar("0.7")
c1Entry <- tkentry(top, width = "4", textvariable = c1Var)
c2Var <- tclVar("1.1")
c2Entry <- tkentry(top, width = "4", textvariable = c2Var)
ecVar <- tclVar("0.05")
ecEntry=tkscale(top, from=0, to=10, showvalue=TRUE, variable=ecVar, resolution=0.01, orient="horizontal")
dstVar <- tclVar("0.05")
dstEntry=tkscale(top, from=0, to=1, showvalue=TRUE, variable=dstVar, resolution=0.01, orient="horizontal")
pErVar <- tclVar("0.20")
pErEntry <-tkscale(top, from=0, to=1, showvalue=TRUE, variable=pErVar, resolution=0.01, orient="horizontal")
pScVar <- tclVar("0.10")
pScEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=pScVar, resolution=0.01, orient="horizontal")
pMxVar <- tclVar("0.10")
pMxEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=pMxVar, resolution=0.01, orient="horizontal")
pRsVar <- tclVar("0.10")
pRsEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=pScVar, resolution=0.01, orient="horizontal")
	onOK <- function() 
	{
        closeDialog()
        tmax=as.numeric(tclvalue(tmaxVar))
        cl=as.numeric(tclvalue(clVar))
        ln=as.numeric(tclvalue(lnVar))
        npatch=cl*ln
			if (sum(is.na(c(tmax,npatch)))>0 || tmax <= 0 || npatch <= 0) 
          {
            errorCondition("Number of simulations, columns and rows must be positive integers")
            return()
          }
        c1 <- as.numeric(tclvalue(c1Var))
        c2 <- as.numeric(tclvalue(c2Var))
        ec=as.numeric(tclvalue(ecVar))
        dst=as.numeric(tclvalue(dstVar))
			if (sum(is.na(c(c1,c2)))>0 || c1 <= 0 || c2 <= 0) 
          {
            errorCondition(message = "Colonization rate for both species must be positive ")
            return()
          }
        pEr <- as.numeric(tclvalue(pErVar))
        pSc <- as.numeric(tclvalue(pScVar))
        pMx <- as.numeric(tclvalue(pMxVar))
        pRs <- as.numeric(tclvalue(pRsVar))
        ptot<- pEr+pSc+pMx+pRs
        if (ptot > 1) 
        {
            errorCondition(message = "Proportion of patches occupied should sum less than one\n VALUES ADJUSTED BY TOTAL LESS 10% LEFT EMPTY")
        pEr=(pEr/ptot)*0.9 
        pSc=(pSc/ptot)*0.9 
        pMx=(pMx/ptot)*0.9 
        pRs=(pRs/ptot)*0.9
        }
 ############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("regNicho(tmax= ",tmax, ", rw= ",ln, ", cl = ", cl,", c1 = ", c1,", c2 = ", c2,", ec = ", ec,", dst = ", dst,", er = ", pEr,", sc =", pSc,", mx =", pMx,", rs =", pRs,")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue, " <- regNicho(tmax= ",tmax, ", rw= ",ln, ", cl = ", cl,", c1 = ", c1,", c2 = ", c2,", ec = ", ec,", dst = ", dst,", er = ", pEr,", sc =", pSc,", mx =", pMx,", rs =", pRs,")", sep = "")
		  }
#test1=regNicho(tmax=50, ln=100, cl=100, c1=0.2, c2=0.8, ec=0.5, m=0.04,  Er=0.08, Sc=0.02, Mx=0, Rs=0)
########
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "regNicho")
tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
##
tkgrid(tklabel(top, text="Simulation Arena Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns"), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows"), lnEntry, sticky = "e")
#
tkgrid(tklabel(top, text="Initial Stages Proportions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Early Stage (only sp2) "), pErEntry, sticky = "se")
tkgrid(tklabel(top, text = "Susceptible (only sp1)  "), pScEntry, sticky = "se")
tkgrid(tklabel(top, text = "Mixed (sp1 and sp2)  "), pMxEntry, sticky = "se")
tkgrid(tklabel(top, text = "Resistant (sp1)"), pRsEntry, sticky = "se")
#
tkgrid(tklabel(top, text="Colonization rates :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Better competitor (sp1) "), c1Entry, sticky = "e")
tkgrid(tklabel(top, text = "Inferior competitor (sp2)  "), c2Entry, sticky = "e")
tkgrid(tklabel(top, text="General parameters:", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Competitive exclusion:   "), ecEntry, sticky = "se")
tkgrid(tklabel(top, text = "Disturbance (mortality):  "), dstEntry, sticky = "se")
#
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
#
tkgrid.configure(pErEntry, sticky = "w")
tkgrid.configure(pScEntry, sticky = "w")
tkgrid.configure(pMxEntry, sticky = "w")
tkgrid.configure(pRsEntry, sticky = "w")
#
tkgrid.configure(c1Entry, sticky = "w")
tkgrid.configure(c2Entry, sticky = "w")
tkgrid.configure(ecEntry, sticky = "w")
tkgrid.configure(dstEntry, sticky = "w")
dialogSuffix(rows = 13, columns = 2, focus = tmaxEntry)
}
######################
#############
#teste1=comCompete(tmax=200,ln=100,cl=100, rq=10, fi=0.2, fsp1=0.2, pe=0.04, fr=0, int=0)
comcompDbox<-function () 
{
initializeDialog(title = gettextRcmdr("Trade-off"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
tmaxVar <- tclVar("100")
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
####
rqVar <- tclVar("10")
fiVar <- tclVar("0.10")
fiEntry=tkscale(top, from=0, to=1, showvalue=TRUE, variable=fiVar, resolution=0.01, orient="horizontal")
fsp1Var <- tclVar("0.20")
peVar <- tclVar("0.10")
frVar <- tclVar("0.10")
frEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=frVar, resolution=0.01, orient="horizontal")
intVar <- tclVar("0.10")
intEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=intVar, resolution=0.01, orient="horizontal")
command=paste(paste("gr.toff(rq = ", as.numeric(tclvalue(rqVar)), ", fsp1 = ", as.numeric(tclvalue(fsp1Var)), ", pe = ", as.numeric(tclvalue(peVar)), ",add=FALSE)"))
doItAndPrint(command)
	set.gtoff=function(...)
	{
	#command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", nrows,", byrow=TRUE)", sep="") 
	#gr.toff=function(riq, fsp1,pe,add=FALSE,...)
	command=paste("gr.toff(rq = ", as.numeric(tclvalue(rqVar)), ", fsp1 = ", as.numeric(tclvalue(fsp1Var)), ", pe = ", as.numeric(tclvalue(peVar)), ",add=TRUE)")
	doItAndPrint(command)
	}
rqEntry <- tkscale(top, from=2, to=30, showvalue=TRUE, variable=rqVar, resolution=1, orient="horizontal", command=set.gtoff)
fsp1Entry <-tkscale(top, from=0, to=1, showvalue=TRUE, variable=fsp1Var, resolution=0.01, orient="horizontal",command=set.gtoff)
peEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=peVar, resolution=0.01, orient="horizontal", command=set.gtoff)
#cantoVar <- tclVar("1")
#cantoBox <- tkcheckbutton(top, variable = cantoVar)
onOK <- function() 
	{
	command="dev.off(dev.cur()); dev.new()"
	doItAndPrint(command)
	closeDialog()
	tmax=as.numeric(tclvalue(tmaxVar))
	cl=as.numeric(tclvalue(clVar))
	ln=as.numeric(tclvalue(lnVar))
	npatch=cl*ln
		if (sum(is.na(c(tmax,npatch)))>0 || tmax <= 0 || npatch <= 0)
		{
		errorCondition("Number of simulations, columns and rows must be positive integers")
		return()
		}
	rq <- as.numeric(tclvalue(rqVar))
	fi <- as.numeric(tclvalue(fiVar))
	fsp1 <- as.numeric(tclvalue(fsp1Var))
	pe<- as.numeric(tclvalue(peVar))
	fr <- as.numeric(tclvalue(frVar))
	int <- as.numeric(tclvalue(intVar))
#	cantoVF <- as.logical(as.numeric(tclvalue(cantoVar)))
############ Comando
##
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("comCompete(tmax= ",tmax, ", rw= ",ln, ", cl = ", cl,", S = ", rq,", fi = ", fi,", fsp1 = ", fsp1,", pe = ", pe,", fr = ", fr,", int =", int,")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue, "<- comCompete(tmax= ",tmax, ",rw= ",ln, ", cl = ", cl,", S = ", rq,", fi = ", fi,", fsp1 = ", fsp1,", pe = ", pe,", fr = ", fr,", int =", int,")", sep = "")
		  }
########comCompete(tmax=200,ln=100,cl=100, rq=10, fi=0.2, fsp1=0.2, pe=0.04, fr=0, int=0)
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "comCompete")
tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
##
tkgrid(tklabel(top, text="Simulation Arena Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns"), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows"), lnEntry, sticky = "e")
#
tkgrid(tklabel(top, text="Initial Parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Occupied patches "), fiEntry, sticky = "se")
tkgrid(tklabel(top, text = "Number of species "), rqEntry, sticky = "se")
tkgrid(tklabel(top, text = "Best competitor abundance (sp1) "), fsp1Entry, sticky = "se")
tkgrid(tklabel(top, text = "Mortality rate  "), peEntry, sticky = "se")
#
tkgrid(tklabel(top, text="Disturbance :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Frequency "), frEntry, sticky = "e")
tkgrid(tklabel(top, text = "Intensity "), intEntry, sticky = "e")
#tkgrid(tklabel(top, text = "Show simulation frames"), cantoBox, sticky = "e")
#
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
#
tkgrid.configure(rqEntry, sticky = "w")
tkgrid.configure(fiEntry, sticky = "w")
tkgrid.configure(fsp1Entry, sticky = "w")
tkgrid.configure(peEntry, sticky = "w")
#
tkgrid.configure(frEntry, sticky = "w")
tkgrid.configure(intEntry, sticky = "w")
#tkgrid.configure(cantoBox, sticky = "w")
dialogSuffix(rows = 12, columns = 2, focus = tmaxEntry)
}
#########################
