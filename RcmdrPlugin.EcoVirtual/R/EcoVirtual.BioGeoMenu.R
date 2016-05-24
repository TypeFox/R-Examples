### biogeographical functionsv ###

############################
############################
archipDbox<-function () 
{
initializeDialog(title = gettextRcmdr("The Archipelago"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
tmaxVar <- tclVar("100")
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
nIslVar <- tclVar("10")
nIslEntry <- tkentry(top, width = "3", textvariable = nIslVar)
ctVar <- tclVar("100")
ctEntry <- tkentry(top, width = "4", textvariable = ctVar)
ar.minVar <- tclVar("10")
ar.maxVar <- tclVar("100")
ar.minEntry <- tkscale(top, from=1, to=100, showvalue=TRUE, variable=ar.minVar, resolution=10, orient="horizontal")
ar.maxEntry <- tkscale(top, from=100, to=1000, showvalue=TRUE, variable=ar.maxVar, resolution=100, orient="horizontal")
rqVar <- tclVar("100")
fsp1Var <- tclVar("0.00")
####
command=paste(paste("gr.abund(rq = ", as.numeric(tclvalue(rqVar)), ", fsp1 = ", as.numeric(tclvalue(fsp1Var)), ",add=FALSE)"))
doItAndPrint(command)
########################
	set.gr.abund=function(...)
	{
	#command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", nrows,", byrow=TRUE)", sep="") 
	#gr.toff=function(riq, fsp1,pe,add=FALSE,...)
	command=paste("gr.abund(rq = ", as.numeric(tclvalue(rqVar)), ", fsp1 = ", as.numeric(tclvalue(fsp1Var)),  ",add=TRUE)")
	doItAndPrint(command)
	}
rqEntry <- tkscale(top, from=10, to=10000, showvalue=TRUE, variable=rqVar, resolution=10, orient="horizontal", command=set.gr.abund)
fsp1Entry <-tkscale(top, from=0, to=1, showvalue=TRUE, variable=fsp1Var, resolution=0.01, orient="horizontal",command=set.gr.abund)
cantoVar <- tclVar("1")
cantoBox <- tkcheckbutton(top, variable = cantoVar)
################ passando as variáveis
onOK <- function() 
	{
	command="dev.off(dev.cur()); dev.new()"
	doItAndPrint(command)
	closeDialog()
	tmax=as.numeric(tclvalue(tmaxVar))
	nIsl=as.numeric(tclvalue(nIslVar))
	ct=as.numeric(tclvalue(ctVar))
	ar.min=as.numeric(tclvalue(ar.minVar))
	ar.max=as.numeric(tclvalue(ar.maxVar))
		if (ar.min>ar.max)
		{
		errorCondition("Maximum area < minimum area ")
		return()
		}
	rq <- as.numeric(tclvalue(rqVar))
	fsp1 <- as.numeric(tclvalue(fsp1Var))
	rank <- 1:rq
	abund <- fsp1*(1-fsp1)^(rank-1)
	cantoVF <- as.logical(as.numeric(tclvalue(cantoVar)))
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("archip(n.isl= ", nIsl, ",ar.min= ", ar.min,",ar.max= ",ar.max, ",S= ",rq, ", seed.rain = ", ct,", abund = ", fsp1,", tmax = ", tmax,", anima = ", cantoVF,")" , sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue, "<- archip(n.isl= ", nIsl, ",ar.min= ", ar.min,",ar.max= ",ar.max, ",S= ",rq, ", seed.rain= ", ct,", abund = ", fsp1,", tmax = ", tmax,", anima = ", cantoVF,")" , sep = "")
		  }
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
############ 
OKCancelHelp(helpSubject = "archip")
tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
###########
tkgrid(tklabel(top, text="Archipelago Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Number of Islands"), nIslEntry, sticky = "e")
tkgrid(tklabel(top, text = "Smaller Island (square side)"), ar.minEntry, sticky = "e")
tkgrid(tklabel(top, text = "Bigger Island (square side)"), ar.maxEntry, sticky = "e")
#########
tkgrid(tklabel(top, text="Mainland Parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Number of species "), rqEntry, sticky = "se")
tkgrid(tklabel(top, text = "Abundance dominance"), fsp1Entry, sticky = "se")
tkgrid(tklabel(top, text = "Seed rain size "), ctEntry, sticky = "se")
tkgrid(tklabel(top, text = "Show simulation frames"), cantoBox, sticky = "e")
#
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(nIslEntry, sticky = "w")
tkgrid.configure(ar.minEntry, sticky = "w")
tkgrid.configure(ar.maxEntry, sticky = "w")

#
tkgrid.configure(rqEntry, sticky = "w")
tkgrid.configure(fsp1Entry, sticky = "w")
tkgrid.configure(ctEntry, sticky = "w")
#
tkgrid.configure(cantoBox, sticky = "w")
dialogSuffix(rows = 10, columns = 2, focus = tmaxEntry)
}


######################################
######################################
ColExtDbox<-function () 
{
#animaColExt(minimo=0.01, maximo=1, interv=0.01, Ext="crs", Col="dcr")
initializeDialog(title = gettextRcmdr("Colonization/Extinction Balance"))
#### 
#radioButtons(top, "col", buttons=c("crs", "dcr","fix" ), 
#		labels=gettextRcmdr(c("Increase", "Decrease", "Fixed")), title=gettextRcmdr("Colonization"))
#radioButtons(top, "ext", buttons=c("crs", "dcr","fix" ), 
#		labels=gettextRcmdr(c("Increase", "Decrease", "Fixed")), title=gettextRcmdr("Colonization"))
ColVar <- tclVar("crs")
    ColCrsButton <- ttkradiobutton(top, variable=ColVar, value="crs")
    ColDcrButton <- ttkradiobutton(top, variable=ColVar, value="dcr")
    ColFixButton <- ttkradiobutton(top, variable=ColVar, value="fix")
ExtVar <- tclVar("fix")
    ExtCrsButton <- ttkradiobutton(top, variable=ExtVar, value="crs")
    ExtDcrButton <- ttkradiobutton(top, variable=ExtVar, value="dcr")
    ExtFixButton <- ttkradiobutton(top, variable=ExtVar, value="fix")
onOK <- function(){
        closeDialog()
         Col <- tclvalue(ColVar)
#			Col <- if (col == "crs") "crs"
#			else if (col == "dcr") "dcr"
#			else  "fix"
#			
			Ext <- tclvalue(ExtVar)
#			Ext <- if (ext == "crs") "crs"
#			else if (ext == "dcr") "dcr"
#			else  "fix"
#        
        doItAndPrint(paste("animaColExt(Ext= '", Ext,  "', Col ='", Col,"')", sep=""))
        }
OKCancelHelp(helpSubject = "animaColExt")
###########
tkgrid(tklabel(top, text="Colonization  :", fg="blue"), sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Increase")), ColCrsButton, sticky="e")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Decrease")), ColDcrButton, sticky="e")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Stable")), ColFixButton, sticky="e")
tkgrid(tklabel(top, text="Extinction  :", fg="blue"), sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Increase")), ExtCrsButton, sticky="e")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Decrease")), ExtDcrButton, sticky="e")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Stable")), ExtFixButton, sticky="e")    
################################################################
tkgrid(buttonsFrame, sticky="w", columnspan=2)
    tkgrid.configure(ColCrsButton, sticky="w")
    tkgrid.configure(ColDcrButton, sticky="w")
    tkgrid.configure(ColFixButton, sticky="w")
    tkgrid.configure(ExtCrsButton, sticky="w")
    tkgrid.configure(ExtDcrButton, sticky="w")
    tkgrid.configure(ExtFixButton, sticky="w")
    dialogSuffix(rows=9, columns=2, focus=ColCrsButton)
}
#########################################
#########################################
#bioGeoIsl=function(area , dist , P , peso.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01)
bioGeoIslDbox= function () 
{
    env <- environment()
    initializeDialog(title = gettextRcmdr("Island Biogeographical Model"))
    dsname <- tclVar("Do_Not_Save")
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    PVar <- tclVar("100")
    PEntry <- tkentry(top, width = "4", textvariable = PVar)
    #ncolsVar <- tclVar("2")
    #ncolsValue <- tkentry(top, width = "4", textvariable = PVar)
    bVar <- tclVar("-0.01")
    bEntry <- tkscale(top, from = -1, to = 0, showvalue = TRUE, 
        variable = bVar, resolution = 0.01, orient = "horizontal")
    dVar <- tclVar("-0.01")
    dEntry <- tkscale(top, from = -1, to = 0, showvalue = TRUE, 
        variable = dVar, resolution = 0.01, orient = "horizontal")
    fVar <- tclVar("0.01")
    fEntry <- tkscale(top, from = 0, to = 1, showvalue = TRUE, 
        variable = fVar, resolution = 0.01, orient = "horizontal")
    hVar <- tclVar("0.01")
    hEntry <- tkscale(top, from = 0, to = 1, showvalue = TRUE, 
        variable = hVar, resolution = 0.01, orient = "horizontal")
    peso.AVar <- tclVar("0.5")
    peso.AEntry <- tkscale(top, from = 0, to = 1, showvalue = TRUE, 
        variable = peso.AVar, resolution = 0.01, orient = "horizontal")
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir = env)
    setUpTable <- function(...) {
        tkdestroy(get(".tableFrame", envir = env))
        assign(".tableFrame", tkframe(outerTableFrame), envir = env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(2)
        make.col.names <- "labelRcmdr(.tableFrame, text='')"
        assign(".colname.1", tclVar("Distance"), envir = env)
        assign(".colname.2", tclVar("Size"), envir = env)
        	for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep = "")
            make.col.names <- paste(make.col.names, ", ", "ttkentry(.tableFrame, width='8', textvariable= ",  col.varname, ")", sep = "")
        	}
        eval(parse(text = paste("tkgrid(", make.col.names, ")", sep = "")), envir = env)
        for (i in 1:nrows) {
            varname <- paste(".tab.", i, ".1", sep = "")
            assign(varname, tclVar(""), envir = env)
            row.varname <- paste(".rowname.", i, sep = "")
            assign(row.varname, tclVar(paste("Island", i, sep = "_")), 
                envir = env)
            make.row <- paste("ttkentry(.tableFrame, width='8', textvariable=", 
                row.varname, ")", sep = "")
            make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='8', textvariable=", 
                varname, ")", sep = "")
            for (j in 2:ncols) {
                varname <- paste(".tab.", i, ".", j, sep = "")
                assign(varname, tclVar(""), envir = env)
                make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='8', textvariable=", 
                  varname, ")", sep = "")
            }
            eval(parse(text = paste("tkgrid(", make.row, ")", 
                sep = "")), envir = env)
        }
        tkgrid(get(".tableFrame", envir = env), sticky = "w")
    }
    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("2")
    rowsSlider <- tkscale(rowColFrame, from = 2, to = 10, showvalue = FALSE, variable = rowsValue, resolution = 1, orient = "horizontal", command = setUpTable)
    rowsShow <- labelRcmdr(rowColFrame, textvariable = rowsValue, width = 2, justify = "right")
    onOK <- function() {
        PVar <- as.numeric(tclvalue(PVar))
        bVar <- as.numeric(tclvalue(bVar))
        dVar <- as.numeric(tclvalue(dVar))
        fVar <- as.numeric(tclvalue(fVar))
        hVar <- as.numeric(tclvalue(hVar))
        peso.AVar <- as.numeric(tclvalue(peso.AVar))
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(2)
        cell <- 0
        dist <- rep(0, nrows)
        size <- rep(0, nrows)
        row.names <- rep("", nrows)
        col.names <- rep("", nrows)
        for (i in 1:nrows) row.names[i] <- eval(parse(text = paste("tclvalue(", 
            paste(".rowname.", i, sep = ""), ")", sep = "")))
        for (j in 1:ncols) col.names[j] <- eval(parse(text = paste("tclvalue(", 
            paste(".colname.", j, sep = ""), ")", sep = "")))
        for (i in 1:nrows) {
            cell <- cell + 1
            varname1 <- paste(".tab.", i, ".1", sep = "")
            varname2 <- paste(".tab.", i, ".2", sep = "")
            dist[cell] <- as.numeric(eval(parse(text = paste("tclvalue(", 
                varname1, ")", sep = ""))))
            size[cell] <- as.numeric(eval(parse(text = paste("tclvalue(", 
                varname2, ")", sep = ""))))
        }
        .size<-paste("c(", paste(size, collapse = ","), ")", sep = "")
#          command <- paste("c(", paste(size, collapse = ","), ")", 
#              sep = "")
#          assign(".size", justDoIt(command), envir = env)
#         logger(paste(".size <- ", command, sep = ""))
#         doItAndPrint(".size  # sizes")
#	command <- paste("c(", paste(size, collapse = ","), ")", sep = "")
#         command <- paste("c(", paste(dist, collapse = ","), ")", 
#             sep = "")
#         assign(".dist", justDoIt(command), envir = .GlobalEnv)
#         logger(paste(".dist <- ", command, sep = ""))
#         doItAndPrint(".dist  # distances")
	.dist<-paste("c(", paste(dist, collapse = ","), ")", sep = "")
	dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue == "") {
            command <- paste("bioGeoIsl(area=",.size, ",dist  =", .dist,", P = ", 
                PVar, ", b.e = ", bVar, ", d.i = ", dVar, ",f.i =", 
                fVar, ",h.e= ", hVar, ")", sep = "")
        }
        else {
            command <- paste(dsnameValue, "<- bioGeoIsl(area=", 
                .size, ", dist  = ", .dist, ", P = ", PVar, ", b.e = ", 
                bVar, ", d.i = ", dVar, ",f.i =", fVar, ",h.e= ", hVar, 
                ")", sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "bioGeoIsl")
    tkgrid(tklabel(top, text = "Enter name for simulation data set :", 
        fg = "blue"), sticky = "w")
    tkgrid(entryDsname, sticky = "e")
    tkgrid(tklabel(top, text = "Models Parameters :", fg = "blue"), 
        sticky = "w")
    tkgrid(tklabel(top, text = "Mainland Number of Species"), PEntry, sticky = "w")
    tkgrid(tklabel(top, text = "Extinction/Area coefficient "), bEntry, sticky = "w")
    tkgrid(tklabel(top, text = "Extinction/Distance coefficient "),  hEntry, sticky = "w")
    tkgrid(tklabel(top, text = "Colonization/Area coefficient "),  fEntry, sticky = "w")
    tkgrid(tklabel(top, text = "Colonization/Distance coefficient "), dEntry, sticky = "w")
    tkgrid(tklabel(top, text = "Ratio Area/Distance effect"),  peso.AEntry, sticky = "w")
    tkgrid(labelRcmdr(top, text = gettextRcmdr("Island Size and Distance: "), fg = "blue"), sticky = "w")
    tkgrid(labelRcmdr(rowColFrame, text = gettextRcmdr("Number of Islands:")),  rowsSlider, rowsShow, sticky = "w")
    tkgrid(rowColFrame, sticky = "w")
    tkgrid(outerTableFrame, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 10, columns = 2)
}
#bioGeoIslDbox()

##########################################
##########################################
gr.abund=function(rq, fsp1, add=FALSE,...)
{
	rank=1:rq
	px= fsp1*(1-fsp1)^(rank-1)
		if(add==FALSE)
		{
		toff<-dev.new( width=5, height=5)
		}
	old<-par(mar=c(3,3,3,3))
	plot(px~rank, ylim=c(0,fsp1),type="b", bty="n",  ann=FALSE, cex.axis=0.8)
	par(new=TRUE)
	mtext("Specie abundance rank", 1, 2, cex=0.9)
	mtext("Abundance", 2, 2, cex=0.9)
	mtext("Relative Species Abundance ", 3, 0, cex=1.2)
	par(old)
}


#############################
#############################
randWalkDbox=function () 
{
    initializeDialog(title = gettextRcmdr("Random Walk Simulation"))
    dsname <- tclVar("Do_Not_Save")
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    nsVar <- tclVar("10")
    nsEntry <- tkentry(top, width = "4", textvariable = nsVar)
    stVar <- tclVar("1")
    stEntry <- tkentry(top, width = "2", textvariable = stVar)
    midVar <- tclVar("200")
    midEntry <- tkentry(top, width = "3", textvariable = midVar)
    cantoVar <- tclVar("0")
    cantoBox <- tkcheckbutton(top, variable = cantoVar)
    tfVar <- tclVar("1e5")
    tfVarSlider <- tkentry(top, width = "9", textvariable = tfVar)
    #tkscale(top, from = 1e3, to = 1e6, showvalue = TRUE, 
    #    variable = tfVar, resolution = 1e3, orient = "horizontal")
    onOK <- function() {
#        closeDialog()
        tf <- round(as.numeric(tclvalue(tfVar)))
        cantoVF <- as.logical(as.numeric(tclvalue(cantoVar)))
        ns <- round(as.numeric(tclvalue(nsVar)))
        st <- round(as.numeric(tclvalue(stVar)))
        mid <- as.numeric(tclvalue(midVar))
        dsnameValue <- trim.blanks(tclvalue(dsname))
##        randWalk <- function(n=1,step=1,ciclo=1e5,x1max=NULL, alleq=FALSE)
        if (dsnameValue == "Do_Not_Save" | dsnameValue == "") {
            command <- paste("randWalk(S = ", ns, ", step = ", 
                st, ", tmax = ", tf, ", x1max =", mid, ", alleq = ", cantoVF, 
                ")", sep = "")
        }
        else {
            command <- paste(dsnameValue, "<-randWalk(S = ", ns, ", step = ", 
                st, ", tmax = ", tf, ", x1max =", mid, ", alleq = ", cantoVF, 
                ")", sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "randWalk")
    tkgrid(tklabel(top, text = "Enter name for data set: "), 
        entryDsname, sticky = "e")
    tkgrid(tklabel(top, text = "Random Walk Parameters", 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(top, text = "Number of Species  "), nsEntry, sticky = "e")
    tkgrid(tklabel(top, text = "Step size  "), stEntry, sticky = "e")
    tkgrid(tklabel(top, text = "Maximum initial distance  "), midEntry, sticky = "e")
    tkgrid(tklabel(top, text = "Initial distance equal"), cantoBox, sticky = "e")
    tkgrid(tklabel(top, text = "Final time  "), tfVarSlider, sticky = "e")


    tkgrid.configure(entryDsname, sticky = "sw")
    tkgrid.configure(nsEntry, sticky = "sw")
    tkgrid.configure(stEntry, sticky = "sw")
    tkgrid.configure(midEntry, sticky = "sw")
    tkgrid.configure(tfVarSlider, sticky = "sw")
    tkgrid.configure(cantoBox, sticky = "sw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 8, columns = 2, focus = nsEntry)
}


############################
############################
extGameDbox=function() 
{
    initializeDialog(title = gettextRcmdr("Zero Sum Game"))
    dsname <- tclVar("Do_Not_Save")
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    mVar <- tclVar(20)
    mVarSlider <- tkscale(top, from = 10, to = 100, showvalue = TRUE, 
        variable = mVar, resolution = 1, orient = "horizontal")
    apVar <- tclVar(1)
    apVarSlider <- tkscale(top, from = 1, to = 10, showvalue = TRUE, 
        variable = apVar, resolution = 1, orient = "horizontal")
    mtVar<-tclVar("1")
    mtEntry <- tkentry(top, width = "1", textvariable = mtVar)
    onOK <- function() {
        closeDialog()
        m <- as.numeric(tclvalue(mVar))
        ap <- as.numeric(tclvalue(apVar))
        mt <- as.numeric(tclvalue(mtVar))
        dsnameValue <- trim.blanks(tclvalue(dsname))
##extGame <- function(aposta=1,total=100, tmax=5
        if (dsnameValue == "Do_Not_Save" | dsnameValue == "") {
            command <- paste("extGame(bet = ", ap, ", total = ", 
                m, ", tmax = ", mt,")", sep = "")
        }
        else {
            command <- paste(dsnameValue, "<-extGame(bet = ", ap, ", total = ", m, ", tmax = ", mt,")", sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "extGame")
    tkgrid(tklabel(top, text = "Enter name for data set: "), 
        entryDsname, sticky = "e")
    tkgrid(tklabel(top, text = "Game Parameters", 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(top, text = "Total amount  "), mVarSlider , sticky = "e")
    tkgrid(tklabel(top, text = "bet size  "), apVarSlider , sticky = "e")
    tkgrid(tklabel(top, text = "Maximum game time  "), mtEntry, sticky = "e")
 
    tkgrid.configure(entryDsname, sticky = "sw")
    tkgrid.configure(mVarSlider, sticky = "sw")
    tkgrid.configure(apVarSlider, sticky = "sw")
    tkgrid.configure(mtEntry, sticky = "sw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = mVarSlider)
}

#############################################
#############################################
hubDbox=function () 
{
    initializeDialog(title = gettextRcmdr("Neutral Model Simulation"))
    dsname <- tclVar("Do_Not_Save")
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    nsVar <- tclVar(10)
    nsVarSlider <- tkscale(top, from = 10, to = 1000, showvalue = TRUE, 
        variable = nsVar, resolution = 10, orient = "horizontal") 
    jiVar <- tclVar(10)
    jiVarSlider <- tkscale(top, from = 1, to = 100, showvalue = TRUE, 
        variable = jiVar, resolution = 1, orient = "horizontal")
    dVar <- tclVar("1")    
    dEntry <- tkentry(top, width = "2", textvariable = dVar)    
    cicloVar <- tclVar(1000)
    cicloSlider <- tkscale(top, from = 1e3, to = 1e5, showvalue = TRUE, 
        variable = cicloVar, resolution = 1000, orient = "horizontal")
    migVar <- tclVar("0")
    migBox <- tkcheckbutton(top, variable = migVar)
    spVar <- tclVar("0")
    spBox <- tkcheckbutton(top, variable = spVar)
    animaVar <- tclVar("1")
    animaBox <- tkcheckbutton(top, variable = animaVar)
    onOK <- function() {
#        closeDialog()
        S<- as.numeric(tclvalue(nsVar))
        ji<- as.numeric(tclvalue(jiVar))
        D<- as.numeric(tclvalue(dVar))
        ciclos<-as.numeric(tclvalue(cicloVar))
        migVF <- as.logical(as.numeric(tclvalue(migVar)))
        spVF <- as.logical(as.numeric(tclvalue(spVar)))
        animaVF <- as.logical(as.numeric(tclvalue(animaVar)))
        if(spVF)
        {
        closeDialog()
        hubDbox3()
        stop()
        }
			if(!spVF & migVF)
        {
        closeDialog()
        hubDbox2()
        stop()
        }
        dsnameValue <- trim.blanks(tclvalue(dsname))
##        simHub1=function(S= 100, j=10, D=1, ciclo=1e4){
        if (dsnameValue == "Do_Not_Save" | dsnameValue == "") {
            command <- paste("simHub1(S = ", S, ", j = ", 
                ji, ", D = ", D,", cycles = ", ciclos, ", anima = ",animaVF ,")", sep = "")
        }
        else {
            command <- paste(dsnameValue, "<-simHub1(S = ", S, ", j = ", 
                ji,", D = ", D, ", cycles = ", ciclos, ", anima = ",animaVF , ")", sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "simHub")
    tkgrid(tklabel(top, text = "Enter name for data set: "), 
        entryDsname, sticky = "e")
    tkgrid(tklabel(top, text = "Neutral Model Parameters", 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(top, text = "Number of Species  "), nsVarSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Individuals per species  "), jiVarSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Cycles per Simulation "), cicloSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Number of dead per cycle  "), dEntry, sticky = "e")
    tkgrid(tklabel(top, text = "Immigration  "), migBox, sticky = "e")
    tkgrid(tklabel(top, text = "Immigration and Speciation   "), spBox, sticky = "s")
    tkgrid(tklabel(top, text = "Show simulation frames "), animaBox, sticky = "s")

    tkgrid.configure(entryDsname, sticky = "sw")
    tkgrid.configure(nsVarSlider, sticky = "sw")
    tkgrid.configure(jiVarSlider, sticky = "sw")
    tkgrid.configure(cicloSlider, sticky = "sw")
    tkgrid.configure(dEntry, sticky = "sw")
    tkgrid.configure(migBox, sticky = "sw")
    tkgrid.configure(spBox, sticky = "sw")
    tkgrid.configure(animaBox, sticky = "sw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "e")
    dialogSuffix(rows = 9, columns = 2, focus = nsVarSlider)
}

##############################################
#############################################
hubDbox2=function () 
{
    initializeDialog(title = gettextRcmdr("Neutral Model Simulation"))
    dsname <- tclVar("Do_Not_Save")
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    nsVar <- tclVar(10)
    nsVarSlider <- tkscale(top, from = 10, to = 1000, showvalue = TRUE, 
        variable = nsVar, resolution = 10, orient = "horizontal") 
    jiVar <- tclVar(10)
    jiVarSlider <- tkscale(top, from = 1, to = 100, showvalue = TRUE, 
        variable = jiVar, resolution = 1, orient = "horizontal")
    dVar <- tclVar("1")    
    dEntry <- tkentry(top, width = "2", textvariable = dVar)    
    cicloVar <- tclVar(1000)
    cicloSlider <- tkscale(top, from = 1e3, to = 1e5, showvalue = TRUE, 
        variable = cicloVar, resolution = 1000, orient = "horizontal")
    mig1Var <- tclVar("1e-4")
    mig1Slider <- tkentry(top, width = "10", textvariable = mig1Var)
#    mig1Var <- tclVar(0)
#    mig1Slider <- tkscale(top, from = 0, to = 1, showvalue = TRUE, 
#        variable = mig1Var, resolution = 1e-9, orient = "horizontal")
    animaVar <- tclVar("1")
    animaBox <- tkcheckbutton(top, variable = animaVar)
    onOK <- function() {
#        closeDialog()
        S<- as.numeric(tclvalue(nsVar))
        ji<- as.numeric(tclvalue(jiVar))
        D<- as.numeric(tclvalue(dVar))
        ciclos<-as.numeric(tclvalue(cicloVar))
        animaVF <- as.logical(as.numeric(tclvalue(animaVar)))
        migra<-as.numeric(tclvalue(mig1Var))
        dsnameValue <- trim.blanks(tclvalue(dsname))
##        simHub1=function(S= 100, j=10, D=1, ciclo=1e4){
        if (dsnameValue == "Do_Not_Save" | dsnameValue == "") {
            command <- paste("simHub2(S = ", S, ", j = ", 
                ji,  ", D = ", D, ", cycles = ", ciclos, ", m = ",migra ,", anima = ",animaVF ,")", sep = "")
        }
        else {
            command <- paste(dsnameValue, "<-simHub2(S = ", S, ", j = ", 
                ji,  ", D = ", D,", cycles = ", ciclos, ", m = ",migra, ", anima = ",animaVF , ")", sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "simHub")
    tkgrid(tklabel(top, text = "Enter name for data set: "), 
        entryDsname, sticky = "e")
    tkgrid(tklabel(top, text = "Neutral Model Parameters", 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(top, text = "Number of Species  "), nsVarSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Individuals per species  "), jiVarSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Number of dead per cycle  "), dEntry, sticky = "e")
    tkgrid(tklabel(top, text = "Cycles per Simulation "), cicloSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Immigration rate  "), mig1Slider, sticky = "e")
    tkgrid(tklabel(top, text = "Show simulation frames "), animaBox, sticky = "s")
    tkgrid.configure(entryDsname, sticky = "sw")
    tkgrid.configure(nsVarSlider, sticky = "sw")
    tkgrid.configure(jiVarSlider, sticky = "sw")
    tkgrid.configure(dEntry, sticky = "sw")
    tkgrid.configure(mig1Slider, sticky = "sw")
    tkgrid.configure(cicloSlider, sticky = "sw")
    tkgrid.configure(animaBox, sticky = "sw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "e")
    dialogSuffix(rows = 9, columns = 2, focus = mig1Slider)
}


##############################################
#############################################
#simHub3=function(Sm=200, jm=20, S= 100, j=10, D=1, ciclo=1e4, m=0.01, nu=0.001, anima=TRUE)
hubDbox3=function() 
{
    initializeDialog(title = gettextRcmdr("Neutral Model Simulation"))
    dsname <- tclVar("Do_Not_Save")
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    nsVar <- tclVar(10)
    nsVarSlider <- tkscale(top, from = 10, to = 1000, showvalue = TRUE, 
        variable = nsVar, resolution = 10, orient = "horizontal") 
    jiVar <- tclVar(10)
    jiVarSlider <- tkscale(top, from = 1, to = 100, showvalue = TRUE, 
        variable = jiVar, resolution = 1, orient = "horizontal")
    ## Metacomunity    
    SmVar <- tclVar(10)
    SmSlider <- tkscale(top, from = 10, to = 1000, showvalue = TRUE, 
        variable = SmVar, resolution = 10, orient = "horizontal") 
    jmVar <- tclVar(10)
    jmSlider <- tkscale(top, from = 1, to = 100, showvalue = TRUE, 
        variable = jmVar, resolution = 1, orient = "horizontal")

    dVar <- tclVar("1")    
    dEntry <- tkentry(top, width = "2", textvariable = dVar)    
    cicloVar <- tclVar(1000)
    cicloSlider <- tkscale(top, from = 1e3, to = 1e5, showvalue = TRUE, 
        variable = cicloVar, resolution = 1000, orient = "horizontal")
    mig1Var <- tclVar("1e-4")
    mig1Slider <- tkentry(top, width = "10", textvariable = mig1Var)
    #tkscale(top, from = 0, to = 1, showvalue = TRUE, 
        #variable = mig1Var, resolution = 1e-9, orient = "horizontal")
    nuVar <- tclVar("1e-4")
    nuSlider <- tkentry(top, width = "10", textvariable = nuVar)
    #tkscale(top, from = 0, to = 1, showvalue = TRUE, 
        #variable = nuVar, resolution = 1e-9, orient = "horizontal")
    animaVar <- tclVar("1")
    animaBox <- tkcheckbutton(top, variable = animaVar)
    onOK <- function() {
#        closeDialog()
        S<- as.numeric(tclvalue(nsVar))
        Sm<-as.numeric(tclvalue(SmVar))
        ji<- as.numeric(tclvalue(jiVar))
        jm<- as.numeric(tclvalue(jmVar))
        D<- as.numeric(tclvalue(dVar))
        ciclos<-as.numeric(tclvalue(cicloVar))
        animaVF <- as.logical(as.numeric(tclvalue(animaVar)))
        migra<-as.numeric(tclvalue(mig1Var))
        nu<-as.numeric(tclvalue(nuVar))
        dsnameValue <- trim.blanks(tclvalue(dsname))
##        simHub1=function(S= 100, j=10, D=1, ciclo=1e4){
        if (dsnameValue == "Do_Not_Save" | dsnameValue == "") {
            command <- paste("simHub3(Sm = ", Sm,  ", jm = ", jm,",S = ", S, ", j = ", 
                ji,  ", D = ", D, ", cycles = ", ciclos, ", m = ", migra ,", nu = ", nu ,", anima = ",animaVF ,")", sep = "")
        }
        else {
            command <- paste(dsnameValue, "<-simHub3(Sm = ", Sm,  ", jm = ", jm,",S = ", S, ", j = ", 
                ji,  ", D = ", D,", cycles = ", ciclos, ", m = ",migra,", nu = ", nu , ", anima = ",animaVF , ")", sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "simHub")
    tkgrid(tklabel(top, text = "Enter name for data set: "), 
        entryDsname, sticky = "e")
    tkgrid(tklabel(top, text = "Regional Communnity Parameters", 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(top, text = "Number of Species (Metacommunity)  "), SmSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Individuals per species (Metacommunity)  "), jmSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Immigration rate  "), mig1Slider, sticky = "e")
    tkgrid(tklabel(top, text = "Speciation rate  "), nuSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Local Community Parameters", 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(top, text = "Number of local Species  "), nsVarSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Individuals per species (local)  "), jiVarSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Number of dead per cycle  "), dEntry, sticky = "e")
    tkgrid(tklabel(top, text = "Cycles per Simulation "), cicloSlider, sticky = "e")
    tkgrid(tklabel(top, text = "Show simulation frames "), animaBox, sticky = "s")
    tkgrid.configure(entryDsname, sticky = "sw")
    tkgrid.configure(SmSlider, sticky = "sw")
    tkgrid.configure(jmSlider, sticky = "sw")
    tkgrid.configure(mig1Slider, sticky = "sw")
    tkgrid.configure(nuSlider, sticky = "sw")
    tkgrid.configure(nsVarSlider, sticky = "sw")
    tkgrid.configure(jiVarSlider, sticky = "sw")
    tkgrid.configure(dEntry, sticky = "sw")
    tkgrid.configure(cicloSlider, sticky = "sw")
    tkgrid.configure(animaBox, sticky = "sw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "e")
    dialogSuffix(rows = 10, columns = 2, focus = nuSlider)
}

