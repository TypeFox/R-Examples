#************************************#
# SECOND ATTEMPT TO CREATE POWER.GUI #
#************************************#


require(tcltk)      			# Load the TclTk package
tclRequire("Iwidgets")			# needed for tabbed notebook 
tclRequire("BWidget")			# needed for drop down combo boxes 
# tclRequire("Tktable")			# needed for create tables

power.gui <- tktoplevel()		# Create a new toplevel window
tktitle(power.gui) <- "PoweR"  	# Give the window a title

tkpack(tab.nb <- tkwidget(power.gui, "iwidgets::tabnotebook"))		# tab.nb = tabnotebook


# fonts
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14,weight="bold")


#-------------------------------------------------------------------#
# a menu consists of a horizontal menubar, with associated submenus #
#-------------------------------------------------------------------#


### first create the horizontal menubar 
mb <- tkmenu(power.gui)
tkconfigure(power.gui, menu = mb)


### now put in the menubar some menubuttons 
m.file <- tkmenu(mb, tearoff = FALSE)
m.edit <- tkmenu(mb, tearoff = FALSE)
m.opts <- tkmenu(mb, tearoff = FALSE)
m.help <- tkmenu(mb, tearoff = FALSE)



### now create the submenus 

## the FILE menus 
tkadd(m.file, "command", label = "Load", command = function() {  })
tkadd(m.file, "command", label = "Save", command = function() {  })
tkadd(m.file, "command", label = "Quit", command = function() { tkdestroy(power.gui) })


## the EDIT menus



## the OPTIONS menus



## the HELP menus


# function HELP.PACKAGE   
help.package <- function() {

  # 1st solution by using tktext
  
  # tt <- tktoplevel()
  # tktitle(tt) <- "Help window"
  # toto <- tktext(tt, height=10)
  # tkpack(toto)
  
  # tkinsert(toto,"end",paste("help(package=\"PoweR\")\n"))
 
  # code.package <- tclvalue(tkget(toto,"0.0","end"))
  # e.package <- try(parse(text=code.package))
  ## cat(code.package)
  ## print(eval(e.package))
  # eval(e.package)
  # tkdestroy(tt)
  
  # 2nd solution 
  
  eval(parse(text=paste("help(package=\"PoweR\")\n")))
  
}
# end of HELP.PACKAGE function			

tkadd(m.help, "command", label = "Help", command = help.package)
tkadd(m.help, "command", label = "About", command = function() { tkmessageBox(title = "About", message = "Authors :\n
Pierre Lafaye de Micheaux \n
Viet Anh Tran \n", icon = "info") })
																				   

### continue the creation of the horizontal menubar
tkadd(mb, "cascade", label = "File", menu = m.file)
tkadd(mb, "cascade", label = "Edit", menu = m.edit)
tkadd(mb, "cascade", label = "Options", menu = m.opts)
tkadd(mb, "cascade", label = "Help", menu = m.help)




#------------------------------ END OF MENU ------------------------------#



#-------------#
# Tabnotebook #
#-------------#



### configure size of tab.nb
tkconfigure(tab.nb,                         
            tabpos="n",
            # width=665,
			# height=670,
            width=800,
			height=600,
			angle=0,
            bevelamount=2,
            gap=2,
            margin=2,
            tabborders=0,
            tabbackground="gray",
            background="lightgray",
            backdrop="white")
			
			
### list name of tabs			
nm <- c("Generate sample","Compute statistic","Critical values","Compute power","Examples")



#---------------------------------#
# tab1 = Generate sample from law #
#---------------------------------#


tbn1 <- tclvalue(tkadd(tab.nb, label=nm[1]))
tkpack(tbw1 <- .Tk.newwin(tbn1))
tkpack(fr1 <- tkframe(tbw1))
tkpack(lb1 <- tklabel(fr1, text=paste("Generate random samples from a law added in the package \n"), font=fontHeading))
ID <- paste(tab.nb$ID, evalq(num.subwin <- num.subwin+1, tab.nb$env), sep=".")
win <- .Tk.newwin(ID)
assign(ID, tbw1, envir = tab.nb$env)
assign("parent", tab.nb, envir = tbw1$env)


## 
# The previous example was a very minimalistic approach, with no concern with the actual 
# appearance of the final product. One simple way to improve the appearance is to add a 
# little space around the widgets when we pack them. The tkpack function accepts two arguments, 
# padx= and pady= to control the amount of space around the widgets when they are packed. 
# Each of these arguments takes a scalar or a vector of length 2. When a scalar is given, 
# that much space (in pixels) is added on either side of the widget in the specified dimension. 
# When a list is given the first number is the amount of space to the left (padx) or the top (pady), 
# and the second number is the amount of space to the right (padx) or the bottom (pady). 
# Like most aspects of GUI design, experimentation is often necessary to find the best solution. 
# The following code produces a version of the coin power GUI with additional padding:


## FRAME LEFT.FR ##

left.fr <- tkframe(fr1)

# to use with function LOAD and SAVE
wfile <- ""

law.index <- tclVar("")
n <- tclVar("")
parlaw <- tclVar("NULL")


# function RESET and RESET BUTTON
reset <- function() {

         tclvalue(law.index) <- ""
         tclvalue(n) <- ""
         tclvalue(parlaw) <- "NULL"

}
 
reset.but <- tkbutton(left.fr, text="Reset", command=reset)


# function SUBMIT and SUBMIT BUTTON	   
submit <- function() {

  # delete/clear txt before inserting new lines
  tkdelete(txt,"0.0","end")   
  
  # read data from input
  tkinsert(txt,"end",paste("law.index <-",tclvalue(law.index),"\n",sep=" "))
  
  tkinsert(txt,"end",paste("n <-",tclvalue(n),"\n",sep=" "))
  
  if (tclvalue(parlaw) == "NULL") {
    tkinsert(txt,"end",paste("law.pars <- NULL","\n",sep=" "))
  } else {
    tkinsert(txt,"end",paste("law.pars <- c(",tclvalue(parlaw),") \n",sep=""))
  } 
  
  tkinsert(txt,"end",paste("gensample(law.index,n,law.pars) \n\n"))

}

submit.but <- tkbutton(left.fr, text="Submit", command=submit)
  
  
# Law index, n, parlaw1, parlaw2, parlaw3, parlaw4 entries
f11 <- tkframe(left.fr)
tkpack(tklabel(f11,text='law index',width=10),side='left',pady=c(60,10))
tkpack(tkentry(f11,width=20,textvariable=law.index),side='left',padx=c(0,0),pady=c(60,10))
tkpack(info11.but <- tkbutton(f11, text="?", command=function() {
	tkmessageBox(message="Law index as given by function getindex(). length(law)=1 \n
Example : law.index <- 2 \n
Should write in box : 2") })
	,side='left',padx=c(5,0), pady=c(60,10))

	
f12 <- tkframe(left.fr)
tkpack(tklabel(f12,text='sample size',width=10),side='left',pady=c(20,10))
tkpack(tkentry(f12,width=20,textvariable=n),side='left',padx=c(0,0),pady=c(20,10))
tkpack(info12.but <- tkbutton(f12, text="?", command=function() {
	tkmessageBox(message="Number of observations to generate \n
Example : n <- 1000 \n
Should write in box : 1000") })
	,side='left',padx=c(5,0), pady=c(20,10))
	

f13 <- tkframe(left.fr)
tkpack(tklabel(f13,text='law.pars',width=10),side='left',pady=c(20,10))
tkpack(tkentry(f13,width=20,textvariable=parlaw),side='left',padx=c(0,0),pady=c(20,10))
tkpack(info13.but <- tkbutton(f13, text="?", command=function() {
	tkmessageBox(message="NULL or a vector of length at most 4 containing 4 possible parameters 
to generate random values from distribution law(parlaw[j],j<=4). 
If NULL, the default parameter values for this law wil be used. \n
Example : law.pars <- c(0,1) \n
Should write in box : 0,1") })
	,side='left',padx=c(5,0), pady=c(20,10))


tkpack(f11,side='top')
tkpack(f12,side='top')
tkpack(f13,side='top')


# function HELP.GENSAMPLE  
help.gensample <- function() {

  # 1st solution by using tktext
  
  # tt.gensample <- tktoplevel()
  # tktitle(tt.gensample) <- "Help window"
  # toto.gensample <- tktext(tt.gensample, height=10)
  # tkpack(toto.gensample)
  
  # tkinsert(toto.gensample,"end",paste("help(gensample)\n"))
 
  # code.gensample <- tclvalue(tkget(toto.gensample,"0.0","end"))
  # e.gensample <- try(parse(text=code.gensample))
  # print(eval(e.gensample))
  # tkdestroy(tt.gensample)
  
  # 2nd solution
  print(eval(parse(text=paste("help(gensample)\n"))))
  
}
# end of HELP.GENSAMPLE function	

tkpack(info.but <- tkbutton(left.fr, text="?", command=help.gensample),side='left',padx=c(30,0), pady=c(30,20))

tkpack(reset.but,side='left',padx=c(35,0), pady=c(30,20))
tkpack(submit.but,side='left',padx=c(35,0), pady=c(30,20))

tkpack(left.fr,side='left',padx=c(30,10))


## END OF LEFT.FR ##



## FRAME RIGHT.FR ##


right.fr <- tkframe(fr1)


# Command editor 
cmd.edit <- tklabel(right.fr,text='Command editor',width=20,font=fontTextLabel)
txt <- tktext(right.fr, width=60, height=10)

tkpack(cmd.edit,side='top',padx=c(0,30),pady=c(5,0))
tkpack(txt,side='top',padx=c(10,30),pady=c(10,10))


# function LOAD and LOAD BUTTON
load <- function() {
  file <- tclvalue(tkgetOpenFile())
  if (!length(file)) return()
  chn <- tclopen(file, "r")
  tkinsert(txt, "0.0", tclvalue(tclread(chn)))
  tclclose(chn)
  wfile <<- file
}	     

load.but <- tkbutton(right.fr, text="Load", command=load) 


# function SAVE and SAVE BUTTON
save <- function() {
  file <- tclvalue(tkgetSaveFile(
    initialfile=tclvalue(tclfile.tail(wfile)),
    initialdir=tclvalue(tclfile.dir(wfile))))
  if (!length(file)) return()
  chn <- tclopen(file, "w")
  tclputs(chn, tclvalue(tkget(txt,"0.0","end")))
  tclclose(chn)
  wfile <<- file
}
 
save.but <- tkbutton(right.fr, text="Save", command=save) 


# function RUN and RUN BUTTON
run <- function() {
  code <- tclvalue(tkget(txt,"0.0","end"))
  e <- try(parse(text=code))
  if (inherits(e, "try-error")) {
    tkmessageBox(message="Syntax error",icon="error")
  return()
  }
  cat("Executing from script window:",
      "-----", code, "Result:", sep="\n")
  print(eval(e))
}
 
run.but <- tkbutton(right.fr, text="Run", command=run) 
  
tkpack(load.but,side='left',padx=c(130,0),pady=c(20,20))
tkpack(save.but,side='left',padx=c(50,0),pady=c(20,20))
tkpack(run.but,side='left',padx=c(50,0),pady=c(20,20))


tkpack(right.fr,side='right',padx=c(20,20))


## END OF RIGHT.FR ##




#------------------------------ END OF TAB1 ------------------------------#




#--------------------------#
# tab2 = Compute statistic #
#--------------------------#


tbn2 <- tclvalue(tkadd(tab.nb, label=nm[2]))
tkpack(tbw2 <- .Tk.newwin(tbn2))
tkpack(fr2 <- tkframe(tbw2))
tkpack(lb2 <- tklabel(fr2, text=paste("Perform the test statistic for the given value of stat \n"), font=fontHeading))
ID <- paste(tab.nb$ID, evalq(num.subwin <- num.subwin+1, tab.nb$env), sep=".")
win <- .Tk.newwin(ID)
assign(ID, tbw2, envir = tab.nb$env)
assign("parent", tab.nb, envir = tbw2$env)



## FRAME LEFT2.FR ##

left2.fr <- tkframe(fr2)

# to use with function LOAD and SAVE
wfile2 <- ""

law.index2 <- tclVar("")
n2 <- tclVar("")
parlaw2 <- tclVar("NULL")
stat.index2 <- tclVar("")
level2 <- tclVar("0.05")
crit2L <- tclVar("NULL")
crit2R <- tclVar("NULL")
alter2 <- tclVar("")
parstat2 <- tclVar("NULL")


# function RESET2 and RESET2 BUTTON
reset2 <- function() {
         tclvalue(law.index2) <- ""
         tclvalue(n2) <- ""
         tclvalue(parlaw2) <- "NULL"
         tclvalue(stat.index2) <- ""
         tclvalue(level2) <- "0.05"
         tclvalue(crit2L) <- "NULL"
         tclvalue(crit2R) <- "NULL"
         tclvalue(alter2) <- ""
         tclvalue(parstat2) <- "NULL"
       }

reset2.but <- tkbutton(left2.fr, text="Reset", command=reset2)


# function SUBMIT2 and SUBMIT2 BUTTON	   
submit2 <- function() {	  

  # delete/clear txt2 before inserting new lines
  tkdelete(txt2,"0.0","end")    

  # read data from input
  
  tkinsert(txt2,"end",paste("# Here we begin by generating a sample from law \n"))
  
  tkinsert(txt2,"end",paste("law.index <-",tclvalue(law.index2),"\n",sep=" "))
  
  tkinsert(txt2,"end",paste("n <-",tclvalue(n2),"\n",sep=" "))
  
  # test if parlaw2 is NULL
  if (tclvalue(parlaw2) == "NULL") {
    tkinsert(txt2,"end",paste("law.pars <- NULL","\n",sep=" "))
  } else {
    tkinsert(txt2,"end",paste("law.pars <- c(",tclvalue(parlaw2),") \n",sep=""))
  } 

  tkinsert(txt2,"end",paste("data <- gensample(law.index,n,law.pars)$sample \n"))
  
  tkinsert(txt2,"end",paste("# Now we compute test statistic \n"))
  
  tkinsert(txt2,"end",paste("stat.index <-",tclvalue(stat.index2),"\n",sep=" "))
  
  tkinsert(txt2,"end",paste("levels <- c(",tclvalue(level2),")\n",sep=""))
  
  tkinsert(txt2,"end",paste("critvalL <-",tclvalue(crit2L),"\n",sep=" "))
  
  tkinsert(txt2,"end",paste("critvalR <-",tclvalue(crit2R),"\n",sep=" "))
  
  tkinsert(txt2,"end",paste("alter <-",tclvalue(alter2),"\n",sep=" "))
  
  tkinsert(txt2,"end",paste("stat.pars <- c(",tclvalue(parstat2),")\n",sep=""))
  
  tkinsert(txt2,"end",paste("statcompute(stat.index,data,levels,critvalL,critvalR,alter,stat.pars) \n\n"))
  
}

submit2.but <- tkbutton(left2.fr, text="Submit", command=submit2)


# Law index, n, parlaw1, parlaw2, parlaw3, parlaw4 entries
# Stat index, level, critvalL, critvalR, alter, stat.pars entries
f21 <- tkframe(left2.fr)
tkpack(tklabel(f21,text='law index',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f21,width=20,textvariable=law.index2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info21.but <- tkbutton(f21, text="?", command=function() {
	tkmessageBox(message="Law index as given by function getindex(). length(law)=1 \n
Example : law.index <- 2 \n
Should write in box : 2") })
	,side='left',padx=c(5,0), pady=c(10,5))
	
	
f22 <- tkframe(left2.fr)
tkpack(tklabel(f22,text='sample size',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f22,width=20,textvariable=n2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info22.but <- tkbutton(f22, text="?", command=function() {
	tkmessageBox(message="Number of observations to generate \n
Example : n <- 1000 \n
Should write in box : 1000") })
	,side='left',padx=c(5,0), pady=c(10,5))


f23 <- tkframe(left2.fr)
tkpack(tklabel(f23,text='law.pars',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f23,width=20,textvariable=parlaw2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info23.but <- tkbutton(f23, text="?", command=function() {
	tkmessageBox(message="NULL or a vector of length at most 4 containing 4 possible parameters 
to generate random values from distribution law(parlaw[j],j<=4). 
If NULL, the default parameter values for this law wil be used. \n
Example : law.pars <- c(0,1) \n
Should write in box : 0,1") })
	,side='left',padx=c(5,0), pady=c(10,5))


f24 <- tkframe(left2.fr)
tkpack(tklabel(f24,text='stat index',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f24,width=20,textvariable=stat.index2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info24.but <- tkbutton(f24, text="?", command=function() {
	tkmessageBox(message="One stat index as given by function getindex() \n
Example : stat.index <- 10 \n
Should write in box : 10") })
	,side='left',padx=c(5,0), pady=c(10,5))


f25 <- tkframe(left2.fr)
tkpack(tklabel(f25,text='levels',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f25,width=20,textvariable=level2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info25.but <- tkbutton(f25, text="?", command=function() {
	tkmessageBox(message="Vector of desired significance levels for the test \n
Example : levels <- c(0.05,0.01) \n
Should write in box : 0.05,0.01") })
	,side='left',padx=c(5,0), pady=c(10,5))


f26 <- tkframe(left2.fr)
tkpack(tklabel(f26,text='critvalL',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f26,width=20,textvariable=crit2L),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info26.but <- tkbutton(f26, text="?", command=function() {
	tkmessageBox(message="NULL or vector of left critival values \n
Example : critvalL <- 8.692 \n
Should write in box : 8.692") })
	,side='left',padx=c(5,0), pady=c(10,5))


f27 <- tkframe(left2.fr)
tkpack(tklabel(f27,text='critvalR',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f27,width=20,textvariable=crit2R),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info27.but <- tkbutton(f27, text="?", command=function() {
	tkmessageBox(message="NULL or vector of right critival values \n
Example : critvalR <- 9.286 \n
Should write in box : 9.286") })
	,side='left',padx=c(5,0), pady=c(10,5))


f28 <- tkframe(left2.fr)
tkpack(tklabel(f28,text='alter',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f28,width=20,textvariable=alter2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info28.but <- tkbutton(f28, text="?", command=function() {
	tkmessageBox(message="Type of test for each statistical test : 
0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 
3: bilateral test that rejects H0 only for large values of the test statistic, 
4: bilateral test that rejects H0 only for small values of the test statistic \n
Example : alter <- 3 \n
Should write in box : 3") })
	,side='left',padx=c(5,0), pady=c(10,5))


f29 <- tkframe(left2.fr)
tkpack(tklabel(f29,text='stat.pars',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f29,width=20,textvariable=parstat2),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info29.but <- tkbutton(f29, text="?", command=function() {
	tkmessageBox(message="A vector of parameters. 
If NULL, the default parameter values for this statistic wil be used. \n
Example : stat.pars <- c(2) \n
Should write in box : 2") })
	,side='left',padx=c(5,0), pady=c(10,5))



tkpack(f21,side='top')
tkpack(f22,side='top')
tkpack(f23,side='top')
tkpack(f24,side='top')
tkpack(f25,side='top')
tkpack(f26,side='top')
tkpack(f27,side='top')
tkpack(f28,side='top')
tkpack(f29,side='top')


# function HELP.STATCOMPUTE  
help.statcompute <- function() {

  print(eval(parse(text=paste("help(statcompute)\n")))) 
 
}
# end of HELP.STATCOMPUTE function	

tkpack(info2.but <- tkbutton(left2.fr, text="?", command=help.statcompute),side='left',padx=c(30,0), pady=c(20,20))

tkpack(reset2.but,side='left',padx=c(30,0), pady=c(20,20))
tkpack(submit2.but,side='left',padx=c(30,0), pady=c(20,20))

tkpack(left2.fr,side='left',padx=c(30,10))


## END OF LEFT2.FR ##



## FRAME RIGHT2.FR ##


right2.fr <- tkframe(fr2)


# Command editor 
cmd2.edit <- tklabel(right2.fr,text='Command editor',width=20,font=fontTextLabel)
txt2 <- tktext(right2.fr, width=60, height=19)

tkpack(cmd2.edit,side='top',padx=c(0,25),pady=c(5,0))
tkpack(txt2,side='top',padx=c(10,20),pady=c(10,10))
  

# function LOAD2 and LOAD2 BUTTON
load2 <- function() {
  file2 <- tclvalue(tkgetOpenFile())
  if (!length(file2)) return()
  chn2 <- tclopen(file2, "r")
  tkinsert(txt2, "0.0", tclvalue(tclread(chn2)))
  tclclose(chn2)
  wfile2 <<- file2
}

load2.but <- tkbutton(right2.fr, text="Load", command=load2) 


# function SAVE2 and SAVE2 BUTTON
save2 <- function() {
  file2 <- tclvalue(tkgetSaveFile(
    initialfile=tclvalue(tclfile.tail(wfile2)),
    initialdir=tclvalue(tclfile.dir(wfile2))))
  if (!length(file2)) return()
  chn2 <- tclopen(file2, "w")
  tclputs(chn2, tclvalue(tkget(txt2,"0.0","end")))
  tclclose(chn2)
  wfile2 <<- file2
}
 
save2.but <- tkbutton(right2.fr, text="Save", command=save2) 


# function RUN2 and RUN2 BUTTON
run2 <- function() {
  code2 <- tclvalue(tkget(txt2,"0.0","end"))
  e2 <- try(parse(text=code2))
  if (inherits(e2, "try-error")) {
    tkmessageBox(message="Syntax error",
                 icon="error")
  return()
  }
  cat("Executing from script window:",
      "-----", code2, "Result:", sep="\n")
  print(eval(e2))
}
  
run2.but <- tkbutton(right2.fr, text="Run", command=run2) 
  
tkpack(load2.but,side='left',padx=c(130,0),pady=c(20,20))
tkpack(save2.but,side='left',padx=c(60,0),pady=c(20,20))
tkpack(run2.but,side='left',padx=c(60,0),pady=c(20,20))


tkpack(right2.fr,side='right',padx=c(20,30))


## END OF RIGHT2.FR ##




#------------------------------ END OF TAB2 ------------------------------#



#------------------------#
# tab3 = Critical values #
#------------------------#


tbn3 <- tclvalue(tkadd(tab.nb, label=nm[3]))
tkpack(tbw3 <- .Tk.newwin(tbn3))
tkpack(fr3 <- tkframe(tbw3))
tkpack(lb3 <- tklabel(fr3, text=paste("Computation of critical values for several test statistics \n"), font=fontHeading))
ID <- paste(tab.nb$ID, evalq(num.subwin<-num.subwin+1, tab.nb$env), sep=".")
win <- .Tk.newwin(ID)
assign(ID, tbw3, envir = tab.nb$env)
assign("parent", tab.nb, envir = tbw3$env)



## FRAME LEFT3.FR ##

left3.fr <- tkframe(fr3)

# to use with function LOAD and SAVE
wfile3 <- ""

law.index3 <- tclVar("")
M3 <- tclVar("")
vectn3 <- tclVar("")
stat.index3 <- tclVar("")
level3 <- tclVar("")
alter3 <- tclVar("")
parlaw3 <- tclVar("NULL")
parstat3 <- tclVar("NULL")
model3 <- tclVar("NULL")

# function RESET3 and RESET3 BUTTON
reset3 <- function() {
         tclvalue(law.index3) <- ""
         tclvalue(M3) <- ""
         tclvalue(vectn3) <- ""
         tclvalue(stat.index3) <- ""
         tclvalue(level3) <- ""
         tclvalue(alter3) <- ""
         tclvalue(parlaw3) <- "NULL"
         tclvalue(parstat3) <- "NULL"
         tclvalue(model3) <- "NULL"
       }
 
reset3.but <- tkbutton(left3.fr, text="Reset", command=reset3)


# function SUBMIT3 and SUBMIT3 BUTTON	   
submit3 <- function() {	   
  
  # delete/clear txt3 before inserting new lines
  tkdelete(txt3,"0.0","end")   
  
  # read data from input
  tkinsert(txt3,"end",paste("law.index <-",tclvalue(law.index3),"\n",sep=" "))
  tkinsert(txt3,"end",paste("M <-",tclvalue(M3),"\n",sep=" "))
  tkinsert(txt3,"end",paste("vectn <- c(",tclvalue(vectn3),")\n",sep=""))
  tkinsert(txt3,"end",paste("stat.indices <- c(",tclvalue(stat.index3),")\n",sep=""))
  tkinsert(txt3,"end",paste("levels <- c(",tclvalue(level3),")\n",sep=""))
  
  # tkinsert(txt3,"end",paste("alter <- list(",tclvalue(alter3),")\n",sep=""))
  
  tkinsert(txt3,"end",paste("alter <- list(",paste("stat",as.numeric(strsplit(tclvalue(stat.index3),",")[[1]]),"=",as.numeric(strsplit(tclvalue(alter3),",")[[1]]),sep="",collapse=","),")\n",sep=""))
  
  if (tclvalue(parlaw3) == "NULL") {
    tkinsert(txt3,"end",paste("law.pars <- ",tclvalue(parlaw3),"\n",sep=""))
  } else {
    tkinsert(txt3,"end",paste("law.pars <- c(",tclvalue(parlaw3),")\n",sep=""))
  }
  
  if (tclvalue(parstat3) == "NULL") {
    tkinsert(txt3,"end",paste("parstats <- ",tclvalue(parstat3),"\n",sep=""))
  } else {
    # tkinsert(txt3,"end",paste("parstats <- list(",tclvalue(parstat3),")\n",sep=""))
	tkinsert(txt3,"end",paste("parstats <- list(",paste("stat",as.numeric(strsplit(tclvalue(stat.index3),",")[[1]]),"=c",strsplit(tclvalue(parstat3)," ")[[1]],sep="",collapse=","),")\n",sep=""))
  }
  
  tkinsert(txt3,"end",paste("model <-",tclvalue(model3),"\n",sep=" "))
  tkinsert(txt3,"end",paste("critvalues <- many.crit(law.index,stat.indices,M,vectn,levels,alter,law.pars,parstats,model) \n"))
  tkinsert(txt3,"end",paste("critvalues \n\n"))
  
}

submit3.but <- tkbutton(left3.fr, text="Submit", command=submit3)


# Law index, n, parlaw1, parlaw2, parlaw3, parlaw4 entries
# Stat index, level, M, model, alter, parstat entries
f31 <- tkframe(left3.fr)
tkpack(tklabel(f31,text='law index',width=10),side='left',padx=c(0,0),pady=c(5,5))
tkpack(tkentry(f31,width=20,textvariable=law.index3),side='left',padx=c(0,0),pady=c(5,5))
tkpack(info31.but <- tkbutton(f31, text="?", command=function() {
	tkmessageBox(message="Law index as given by function getindex(). length(law)=1 \n
Example : law.index <- 2 \n
Should write in box : 2") })
	,side='left',padx=c(5,0), pady=c(5,5))

f32 <- tkframe(left3.fr)
tkpack(tklabel(f32,text='M',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f32,width=20,textvariable=M3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info32.but <- tkbutton(f32, text="?", command=function() {
	tkmessageBox(message="Number of Monte Carlo repetitions to use \n 
Example : M <- 10000 \n
Should write in box : 10000") })
	,side='left',padx=c(5,0), pady=c(10,5))

f33 <- tkframe(left3.fr)
tkpack(tklabel(f33,text='sample size',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f33,width=20,textvariable=vectn3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info33.but <- tkbutton(f33, text="?", command=function() {
	tkmessageBox(message="Vector of number of observations for the samples to be generated. length(n)>=1 \n
Example : vectn <- c(10,20,50,100) \n
Should write in box : 10,20,50,100") })
	,side='left',padx=c(5,0), pady=c(10,5))

f34 <- tkframe(left3.fr)
tkpack(tklabel(f34,text='stat indices',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f34,width=20,textvariable=stat.index3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info34.but <- tkbutton(f34, text="?", command=function() {
	tkmessageBox(message="Vector of stat indices as given by function getindex(). length(stat)>=1 \n
Example : stat.indices <- c(10,15) \n
Should write in box : 10,15") })
	,side='left',padx=c(5,0), pady=c(10,5))

f35 <- tkframe(left3.fr)
tkpack(tklabel(f35,text='levels',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f35,width=20,textvariable=level3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info35.but <- tkbutton(f35, text="?", command=function() {
	tkmessageBox(message="Vector of required significance level values. length(level)>=1 \n
Example : levels <- c(0.1,0.05,0.01) \n
Should write in box : 0.1,0.05,0.01") })
	,side='left',padx=c(5,0), pady=c(10,5))

f36 <- tkframe(left3.fr)
tkpack(tklabel(f36,text='alter',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f36,width=20,textvariable=alter3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info36.but <- tkbutton(f36, text="?", command=function() {
	tkmessageBox(message="Named-list with type of test for each statistical test : alter[[\"statj\"]]=0,1,2,3 or 4; 
for each j in stat.indices 
(0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 
3: bilateral test that rejects H0 only for large values of the test statistic, 
4: bilateral test that rejects H0 only for small values of the test statistic) \n
Example : alter <- list(stat10=3,stat15=3) \n
Should write in box : 3,3") })
	,side='left',padx=c(5,0), pady=c(10,5))

f37 <- tkframe(left3.fr)
tkpack(tklabel(f37,text='law.pars',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f37,width=20,textvariable=parlaw3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info37.but <- tkbutton(f37, text="?", command=function() {
	tkmessageBox(message="NULL or a vector of length at most 4 containing 4 possible parameters to 
generate random values from distribution law(parlaw[j],j<=4) \n
Example : law.pars <- c(0,1) \n
Should write in box : 0,1") })
	,side='left',padx=c(5,0), pady=c(10,5))

f38 <- tkframe(left3.fr)
tkpack(tklabel(f38,text='parstats',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f38,width=20,textvariable=parstat3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info38.but <- tkbutton(f38, text="?", command=function() {
	tkmessageBox(message="Named-list of parameter values for each statistic to simulate. 
The names of the list should be statj, j taken in stat.indices. 
If NULL, the default parameter values for these statistics wil be used. \n
Example : parstats <- list(stat42=c(NA),stat47=c(2),stat49=c(0.5)) \n
Should write in box : (NA) (2) (0.5)") })
	,side='left',padx=c(5,0), pady=c(10,5))

f39 <- tkframe(left3.fr)
tkpack(tklabel(f39,text='model',width=10),side='left',padx=c(0,0),pady=c(10,5))
tkpack(tkentry(f39,width=20,textvariable=model3),side='left',padx=c(0,0),pady=c(10,5))
tkpack(info39.but <- tkbutton(f39, text="?", command=function() {
	tkmessageBox(message="If NULL, no model is used. If an integer i>0, the model coded in the C function modelei is used. \n
Example : model <- NULL") })
	,side='left',padx=c(5,0), pady=c(10,5))


tkpack(f31,side='top')
tkpack(f32,side='top')
tkpack(f33,side='top')
tkpack(f34,side='top')
tkpack(f35,side='top')
tkpack(f36,side='top')
tkpack(f37,side='top')
tkpack(f38,side='top')
tkpack(f39,side='top')


# function HELP.MANYCRIT  
help.manycrit <- function() {

  print(eval(parse(text=paste("help(many.crit)\n"))))

}
# end of HELP.MANYCRIT function	

tkpack(info3.but <- tkbutton(left3.fr, text="?", command=help.manycrit),side='left',padx=c(30,0), pady=c(10,20))

tkpack(reset3.but,side='left',padx=c(35,0), pady=c(10,20))
tkpack(submit3.but,side='left',padx=c(35,0), pady=c(10,20))

tkpack(left3.fr,side='left',padx=c(20,20))


## END OF LEFT3.FR ##



## FRAME RIGHT3.FR ##


right3.fr <- tkframe(fr3)


# Command editor 
cmd3.edit <- tklabel(right3.fr,text='Command editor',width=20,font=fontTextLabel)
txt3 <- tktext(right3.fr, width=58, height=18)

tkpack(cmd3.edit,side='top',padx=c(0,30),pady=c(5,0))
tkpack(txt3,side='top',padx=c(15,30),pady=c(10,20))
  

# function LOAD3 and LOAD3 BUTTON
load3 <- function() {
  file3 <- tclvalue(tkgetOpenFile())
  if (!length(file3)) return()
  chn3 <- tclopen(file3, "r")
  tkinsert(txt3, "0.0", tclvalue(tclread(chn3)))
  tclclose(chn3)
  wfile3 <<- file3
}

load3.but <- tkbutton(right3.fr, text="Load", command=load3) 


# function SAVE3 and SAVE3 BUTTON
save3 <- function() {
  file3 <- tclvalue(tkgetSaveFile(
    initialfile=tclvalue(tclfile.tail(wfile3)),
    initialdir=tclvalue(tclfile.dir(wfile3))))
  if (!length(file3)) return()
  chn3 <- tclopen(file3, "w")
  tclputs(chn3, tclvalue(tkget(txt3,"0.0","end")))
  tclclose(chn3)
  wfile3 <<- file3
}
 
save3.but <- tkbutton(right3.fr, text="Save", command=save3) 


# function RUN3 and RUN3 BUTTON
run3 <- function() {
  code3 <- tclvalue(tkget(txt3,"0.0","end"))
  e3 <- try(parse(text=code3))
  if (inherits(e3, "try-error")) {
    tkmessageBox(message="Syntax error",
                 icon="error")
  return()
  }
  cat("Executing from script window:",
      "-----", code3, "Result:", sep="\n")
  print(eval(e3))
}
  
run3.but <- tkbutton(right3.fr, text="Run", command=run3)


# function LATEX.CRIT  
latex.crit <- function() {
 
  tkinsert(txt3,"end",paste("print.critvalues(critvalues,digits=3,latex.output=TRUE) \n"))
  code.crit <- tclvalue(tkget(txt3,"0.0","end"))
  e.crit <- try(parse(text=code.crit))
  # cat(code.crit)
  # print(eval(e.crit))
  eval(e.crit)

}
# end of LATEX.CRIT function		 
  
tkpack(load3.but,side='left',padx=c(100,0),pady=c(20,20))
tkpack(save3.but,side='left',padx=c(50,0),pady=c(20,20))
tkpack(run3.but,side='left',padx=c(50,0),pady=c(20,20))

tkpack(latex3.but <- tkbutton(right3.fr, text="LaTeX", command=latex.crit),side='left',padx=c(50,0), pady=c(20,20))

tkpack(right3.fr,side='right',padx=c(10,15))


## END OF RIGHT3.FR ##




#------------------------------ END OF TAB3 ------------------------------#



#----------------------#
# tab4 = Compute power #
#----------------------#


tbn4 <- tclvalue(tkadd(tab.nb, label=nm[4]))
tkpack(tbw4 <- .Tk.newwin(tbn4))
tkpack(fr4 <- tkframe(tbw4))
tkpack(lb4 <- tklabel(fr4, text=paste("Computation of power and level tables for hypothesis tests \n"), font=fontHeading))
ID <- paste(tab.nb$ID, evalq(num.subwin<-num.subwin+1, tab.nb$env), sep=".")
win <- .Tk.newwin(ID)
assign(ID, tbw4, envir = tab.nb$env)
assign("parent", tab.nb, envir = tbw4$env)



## FRAME LEFT4.FR ##

left4.fr <- tkframe(fr4)

# to use with function LOAD and SAVE
wfile4 <- ""


# definition of variables
law4 <- tclVar("")
parlaw4 <- tclVar("NULL")
vectn4 <- tclVar("")
M4 <- tclVar("")
law.index4 <- tclVar("")
stat.index4 <- tclVar("")
level4 <- tclVar("")
critval4 <- tclVar("NULL")
alter4 <- tclVar("")
parlaws4 <- tclVar("NULL")
parstats4 <- tclVar("NULL")
nbclus4 <- tclVar("1")
model4 <- tclVar("NULL")	


# function RESET4 and RESET4 BUTTON
reset4 <- function() {
		 
tclvalue(law4) <- ""
		 tclvalue(parlaw4) <- "NULL"
		 tclvalue(vectn4) <- ""
         tclvalue(M4) <- ""
		 tclvalue(law.index4) <- ""
         tclvalue(stat.index4) <- ""
		 tclvalue(level4) <- ""
		 tclvalue(critval4) <- "NULL"
		 tclvalue(alter4) <- ""
		 tclvalue(parlaws4) <- "NULL"
		 tclvalue(parstats4) <- "NULL"
		 tclvalue(nbclus4) <- "1"
         tclvalue(model4) <- "NULL"
		
}
 
reset4.but <- tkbutton(left4.fr, text="Reset", command=reset4)


# function SUBMIT4 and SUBMIT4 BUTTON	   
submit4 <- function() {	

  # delete/clear txt4 before inserting new lines
  tkdelete(txt4,"0.0","end")   
  
  # read data from input
  
  # law for critical values
  tkinsert(txt4,"end",paste("# Here we begin by computing the critical values \n",sep=""))
  
  tkinsert(txt4,"end",paste("law.index <-",tclvalue(law4),"\n",sep=" "))
  
  if (tclvalue(parlaw4) == "NULL") {
    tkinsert(txt4,"end",paste("law.pars <- ",tclvalue(parlaw4),"\n",sep=""))
  } else {
    tkinsert(txt4,"end",paste("law.pars <- c(",tclvalue(parlaw4),")\n",sep=""))
  }  
  
  tkinsert(txt4,"end",paste("M <-",tclvalue(M4),"\n",sep=" "))
  
  tkinsert(txt4,"end",paste("vectn <- c(",tclvalue(vectn4),")\n",sep=""))
  
  tkinsert(txt4,"end",paste("stat.indices <- c(",tclvalue(stat.index4),")\n",sep=""))
  
  tkinsert(txt4,"end",paste("levels <- c(",tclvalue(level4),")\n",sep=""))
  
  tkinsert(txt4,"end",paste("alter <- list(",paste("stat",as.numeric(strsplit(tclvalue(stat.index4),",")[[1]]),"=",as.numeric(strsplit(tclvalue(alter4),",")[[1]]),sep="",collapse=","),")\n",sep=""))

  if (tclvalue(parstats4) == "NULL") {
    tkinsert(txt4,"end",paste("parstats <- ",tclvalue(parstats4),"\n",sep=""))
  } else {
    tkinsert(txt4,"end",paste("parstats <- list(",paste("stat",as.numeric(strsplit(tclvalue(stat.index4),",")[[1]]),"=c",strsplit(tclvalue(parstats4)," ")[[1]],sep="",collapse=","),")\n",sep=""))
  }
  
  # compute critical values for future use
  tkinsert(txt4,"end",paste("critval <- many.crit(law.index,stat.indices,M,vectn,levels,alter,law.pars,parstats) \n",sep=""))
  
  tkinsert(txt4,"end",paste("# Now we compute power for hypothesis tests vs alternative laws by using calculated critical values \n",sep=""))
  
  # list of alternative laws
  tkinsert(txt4,"end",paste("law.indices <- c(",tclvalue(law.index4),")\n",sep=""))
  
  if (tclvalue(parlaws4) == "NULL") {
    tkinsert(txt4,"end",paste("parlaws <- ",tclvalue(parlaws4),"\n",sep=""))
  } else {
    tkinsert(txt4,"end",paste("parlaws <- list(",paste("law",as.numeric(strsplit(tclvalue(law.index4),",")[[1]]),"=c",strsplit(tclvalue(parlaws4)," ")[[1]],sep="",collapse=","),")\n",sep=""))
  }    
  
  tkinsert(txt4,"end",paste("nbclus <-",tclvalue(nbclus4),"\n",sep=" "))
    
  tkinsert(txt4,"end",paste("model <-",tclvalue(model4),"\n",sep=" "))
  
  tkinsert(txt4,"end",paste("pow <- powcomp.fast(law.indices,stat.indices,vectn,M,levels,critval,alter,parlaws,parstats,nbclus,model,law.index,law.pars) \n"))
  
  tkinsert(txt4,"end","pow \n\n")
  
}

submit4.but <- tkbutton(left4.fr, text="Submit", command=submit4)


# variables' entries

f41 <- tkframe(left4.fr)
tkpack(tklabel(f41,text='law',width=10),side='left',padx=c(0,0),pady=c(10,0))
tkpack(tkentry(f41,width=20,textvariable=law4),side='left',padx=c(0,0),pady=c(10,0))
tkpack(info41.but <- tkbutton(f41, text="?", command=function() {
	tkmessageBox(message="Law index as given by function getindex(). length(law)=1 \n
Example : law.index <- 2 \n
Should write in box : 2") })
	,side='left',padx=c(5,0), pady=c(10,0))

	
f42 <- tkframe(left4.fr)
tkpack(tklabel(f42,text='law.pars',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f42,width=20,textvariable=parlaw4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info42.but <- tkbutton(f42, text="?", command=function() {
	tkmessageBox(message="NULL or a vector of length at most 4 containing 4 possible parameters to 
generate random values from distribution law(parlaw[j],j<=4) \n
Example : law.pars <- c(0,1) \n
Should write in box : 0,1") })
	,side='left',padx=c(5,0), pady=c(7,0))

	
f43 <- tkframe(left4.fr)
tkpack(tklabel(f43,text='M',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f43,width=20,textvariable=M4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info43.but <- tkbutton(f43, text="?", command=function() {
	tkmessageBox(message="Number of Monte Carlo repetitions to use \n 
Example : M <- 10000 \n
Should write in box : 10000") })
	,side='left',padx=c(5,0), pady=c(7,0))

	
f44 <- tkframe(left4.fr)
tkpack(tklabel(f44,text='sample size',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f44,width=20,textvariable=vectn4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info44.but <- tkbutton(f44, text="?", command=function() {
	tkmessageBox(message="Vector of number of observations for the samples to be generated. length(n)>=1 \n
Example : vectn <- c(10,20,50,100) \n
Should write in box : 10,20,50,100") })
	,side='left',padx=c(5,0), pady=c(7,0))

	
f45 <- tkframe(left4.fr)
tkpack(tklabel(f45,text='stat indices',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f45,width=20,textvariable=stat.index4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info45.but <- tkbutton(f45, text="?", command=function() {
	tkmessageBox(message="Vector of stat indices as given by the function getindex(). length(stat)>=1 \n
Example : stat.indices <- c(10,15) \n
Should write in box : 10,15") })
	,side='left',padx=c(5,0), pady=c(7,0))

	
f46 <- tkframe(left4.fr)
tkpack(tklabel(f46,text='levels',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f46,width=20,textvariable=level4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info46.but <- tkbutton(f46, text="?", command=function() {
	tkmessageBox(message="Vector of required level values. length(level)>=1 \n
Example : levels <- c(0.1,0.05,0.01) \n
Should write in box : 0.1,0.05,0.01") })
	,side='left',padx=c(5,0), pady=c(7,0))

	
f47 <- tkframe(left4.fr)
tkpack(tklabel(f47,text='alter',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f47,width=20,textvariable=alter4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info47.but <- tkbutton(f47, text="?", command=function() {
	tkmessageBox(message="Named-list with type of test for each statistical test : alter[[\"statj\"]]=0,1,2,3 or 4; 
for each j in stat.indices 
(0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 
3: bilateral test that rejects H0 only for large values of the test statistic, 
4: bilateral test that rejects H0 only for small values of the test statistic) \n
Example : alter <- list(stat10=3,stat15=3) \n
Should write in box : 3,3") })
	,side='left',padx=c(5,0), pady=c(7,0))

	
f48 <- tkframe(left4.fr)
tkpack(tklabel(f48,text='parstats',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f48,width=20,textvariable=parstats4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info48.but <- tkbutton(f48, text="?", command=function() {
	tkmessageBox(message="Named-list of parameter values for each statistic to simulate. 
The names of the list should be statj, j taken in stat.indices. 
If NULL, the default parameter values for these statistics wil be used. \n
Example : parstat <- list(stat42=c(NA),stat47=c(2),stat49=c(0.5)) \n
Should write in box : (NA) (2) (0.5)") })
	,side='left',padx=c(5,0), pady=c(7,0))
	

f49 <- tkframe(left4.fr)
tkpack(tklabel(f49,text='law indices',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f49,width=20,textvariable=law.index4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info49.but <- tkbutton(f49, text="?", command=function() {
	tkmessageBox(message="Vector of law indices as given by the function getindex() for alternative. 
length(laws)>=1 \n
Example : law.indices <- c(2,3,5,35) \n
Should write in box : 2,3,5,35") })
	,side='left',padx=c(5,0), pady=c(7,0))
	
	
f410 <- tkframe(left4.fr)
tkpack(tklabel(f410,text='parlaws',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f410,width=20,textvariable=parlaws4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info410.but <- tkbutton(f410, text="?", command=function() {
	tkmessageBox(message="Named-list of parameter values for each law to simulate. 
The names of the list should be lawj, j taken in law.indices. 
If NULL, the default parameter values for these laws wil be used. \n
Example : parslaws <- list(law2=c(0,1),law3=c(9,2),law5=c(8,6)) \n
Should write in box : (0,1) (9,2) (8,6)") })
	,side='left',padx=c(5,0), pady=c(7,0))
	
	
f411 <- tkframe(left4.fr)
tkpack(tklabel(f411,text='nbclus',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f411,width=20,textvariable=nbclus4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info411.but <- tkbutton(f411, text="?", command=function() {
	tkmessageBox(message="Number of slaves to use for the computation on a cluster. 
This needs parallel or Rmpi package to be installed and functionnal on the system. 
Also the daemon mpd sould be started. \n
Example : nbclus <- 1 \n
Should write in box : 1") })
	,side='left',padx=c(5,0), pady=c(7,0))


f412 <- tkframe(left4.fr)
tkpack(tklabel(f412,text='model',width=10),side='left',padx=c(0,0),pady=c(7,0))
tkpack(tkentry(f412,width=20,textvariable=model4),side='left',padx=c(0,0),pady=c(7,0))
tkpack(info412.but <- tkbutton(f412, text="?", command=function() {
	tkmessageBox(message="If NULL, no model is used. If an integer i>0, the model coded in the C function modelei is used. \n
Example : model <- NULL") })
	,side='left',padx=c(5,0), pady=c(7,0))


tkpack(f41,side='top')
tkpack(f42,side='top')
tkpack(f43,side='top')
tkpack(f44,side='top')
tkpack(f45,side='top')
tkpack(f46,side='top')
tkpack(f47,side='top')
tkpack(f48,side='top')
tkpack(f49,side='top')
tkpack(f410,side='top')
tkpack(f411,side='top')
tkpack(f412,side='top')


# function HELP.POWCOMPFAST  
help.powcompfast <- function() {
  
  print(eval(parse(text=paste("help(powcomp.fast)\n"))))

}
# end of HELP.POWCOMPFAST function	

tkpack(info4.but <- tkbutton(left4.fr, text="?", command=help.powcompfast),side='left',padx=c(30,0), pady=c(25,20))

tkpack(reset4.but,side='left',padx=c(35,0), pady=c(25,20))
tkpack(submit4.but,side='left',padx=c(35,0), pady=c(25,20))

tkpack(left4.fr,side='left',padx=c(0,10))


## END OF LEFT4.FR ##



## FRAME RIGHT4.FR ##


right4.fr <- tkframe(fr4)


# Command editor 
cmd4.edit <- tklabel(right4.fr,text='Command editor',width=20,font=fontTextLabel)
txt4 <- tktext(right4.fr, width=65, height=22)

tkpack(cmd4.edit,side='top',padx=c(0,30),pady=c(5,0))
tkpack(txt4,side='top',padx=c(15,30),pady=c(5,10))
  

# function LOAD4 and LOAD4 BUTTON
load4 <- function() {
  file4 <- tclvalue(tkgetOpenFile())
  if (!length(file4)) return()
  chn4 <- tclopen(file4, "r")
  tkinsert(txt4, "0.0", tclvalue(tclread(chn4)))
  tclclose(chn4)
  wfile4 <<- file4
}

load4.but <- tkbutton(right4.fr, text="Load", command=load4) 


# function SAVE4 and SAVE4 BUTTON
save4 <- function() {
  file4 <- tclvalue(tkgetSaveFile(
    initialfile=tclvalue(tclfile.tail(wfile4)),
    initialdir=tclvalue(tclfile.dir(wfile4))))
  if (!length(file4)) return()
  chn4 <- tclopen(file4, "w")
  tclputs(chn4, tclvalue(tkget(txt4,"0.0","end")))
  tclclose(chn4)
  wfile4 <<- file4
}
 
save4.but <- tkbutton(right4.fr, text="Save", command=save4) 


# function RUN4 and RUN4 BUTTON
run4 <- function() {
  code4 <- tclvalue(tkget(txt4,"0.0","end"))
  e4 <- try(parse(text=code4))
  if (inherits(e4, "try-error")) {
    tkmessageBox(message="Syntax error",
                 icon="error")
  return()
  }
  cat("Executing from script window:",
      "-----", code4, "Result:", sep="\n")
  print(eval(e4))
}
  
run4.but <- tkbutton(right4.fr, text="Run", command=run4)


# function LATEX.POW  
latex.pow <- function() {
 
  tkinsert(txt4,"end",paste("print.power(pow,digits=2,latex.output=TRUE) \n"))
  code.pow <- tclvalue(tkget(txt4,"0.0","end"))
  e.pow <- try(parse(text=code.pow))

  eval(e.pow)

}
# end of LATEX.POW function		 
  
tkpack(load4.but,side='left',padx=c(110,0),pady=c(15,20))
tkpack(save4.but,side='left',padx=c(60,0),pady=c(15,20))
tkpack(run4.but,side='left',padx=c(60,0),pady=c(15,20))

tkpack(latex4.but <- tkbutton(right4.fr, text="LaTeX", command=latex.pow),side='left',padx=c(60,0), pady=c(15,20))

tkpack(right4.fr,side='right')


## END OF RIGHT4.FR ##




#------------------------------ END OF TAB4 ------------------------------#




#-----------------#
# tab5 = Examples #
#-----------------#


tbn5 <- tclvalue(tkadd(tab.nb, label=nm[5]))
tkpack(tbw5 <- .Tk.newwin(tbn5))
tkpack(fr5 <- tkframe(tbw5))
tkpack(lb5 <- tklabel(fr5, text=paste("Reproduce examples from articles \n"), font=fontHeading))
ID <- paste(tab.nb$ID, evalq(num.subwin<-num.subwin+1, tab.nb$env), sep=".")
win <- .Tk.newwin(ID)
assign(ID, tbw5, envir = tab.nb$env)
assign("parent", tab.nb, envir = tbw5$env)


# to use with function LOAD and SAVE
wfile5 <- ""


## FRAME LEFT5.FR ##

left5.fr <- tkframe(fr5)


# list box
scr <- tkscrollbar(left5.fr, repeatinterval=5, command=function(...)tkyview(tl,...))
tl<-tklistbox(left5.fr,height=15,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
tkgrid(tklabel(left5.fr, text="List of examples", font=fontTextLabel), padx=c(10,10),pady=c(10,10))
tkgrid(tl,scr)
tkgrid.configure(scr,rowspan=4,sticky="nsw")


# read the directory "examples" for all possible file names 
examples <- sub(".txt","",list.files(paste(path.package("PoweR"),"/examples/",sep=""),pattern = "^.*\\.txt"))


for (i in (1:length(examples))) {

  tkinsert(tl,"end",examples[i])

}
# tkselection.set(tl,2)  # Default is the third from list.  Indexing starts at zero.


# function OnOK et OK BUTTON
OnOK <- function() {

  exampleChoice <- examples[as.numeric(tkcurselection(tl))+1]
 
  
  # delete/clear txt5 before inserting new lines
  tkdelete(txt5,"0.0","end")

  # read the example file corresponding from the directory "PoweR/inst/examples/"
  tmp <- scan(file=paste(paste(path.package("PoweR"),"/examples/",sep=""),exampleChoice,".txt",sep=""),what=c("character"),sep="\n",multi.line=TRUE)
	
  # write lines from example file to text editor
  for (i in 1:length(tmp)) {
	tkinsert(txt5,"end",paste(tmp[i]))
	tkinsert(txt5,"end",paste("\n"))
  }
  
  
}

# end of function OnOK

OK.but <- tkbutton(left5.fr,text="   Submit   ",command=OnOK)
tkgrid(OK.but, pady=c(20,20))



tkpack(left5.fr,side='left',padx=c(20,10))


## END OF LEFT5.FR ##



## FRAME RIGHT5.FR ##


right5.fr <- tkframe(fr5)


# Command editor 
cmd5.edit <- tklabel(right5.fr,text='Command editor',width=20,font=fontTextLabel)
txt5 <- tktext(right5.fr, width=65, height=22)

tkpack(cmd5.edit,side='top',padx=c(0,30),pady=c(5,0))
tkpack(txt5,side='top',padx=c(15,30),pady=c(5,10))
  

# function LOAD5 and LOAD5 BUTTON
load5 <- function() {
  file5 <- tclvalue(tkgetOpenFile())
  if (!length(file5)) return()
  chn5 <- tclopen(file5, "r")
  tkinsert(txt5, "0.0", tclvalue(tclread(chn5)))
  tclclose(chn5)
  wfile5 <<- file5
}

load5.but <- tkbutton(right5.fr, text="Load", command=load5) 


# function SAVE5 and SAVE5 BUTTON
save5 <- function() {
  file5 <- tclvalue(tkgetSaveFile(
    initialfile=tclvalue(tclfile.tail(wfile5)),
    initialdir=tclvalue(tclfile.dir(wfile5))))
  if (!length(file5)) return()
  chn5 <- tclopen(file5, "w")
  tclputs(chn5, tclvalue(tkget(txt5,"0.0","end")))
  tclclose(chn5)
  wfile5 <<- file5
}
 
save5.but <- tkbutton(right5.fr, text="Save", command=save5) 


# function RUN5 and RUN5 BUTTON
run5 <- function() {
  code5 <- tclvalue(tkget(txt5,"0.0","end"))
  e5 <- try(parse(text=code5))
  if (inherits(e5, "try-error")) {
    tkmessageBox(message="Syntax error",
                 icon="error")
  return()
  }
  cat("Executing from script window:",
      "-----", code5, "Result:", sep="\n")
  print(eval(e5))
}
  
run5.but <- tkbutton(right5.fr, text="Run", command=run5)


tkpack(load5.but,side='left',padx=c(160,0),pady=c(15,20))
tkpack(save5.but,side='left',padx=c(60,0),pady=c(15,20))
tkpack(run5.but,side='left',padx=c(60,0),pady=c(15,20))


tkpack(right5.fr,side='right',padx=c(10,20))


## END OF RIGHT5.FR ##




#------------------------------ END OF TAB5 ------------------------------#


tkbind(power.gui, "<Destroy>", function() tkdestroy(tab.nb))
tkselect(tab.nb, 0)



### change theme ? NOT WORKING

# list all themes available, default is 'vista'
# tcl('ttk::style', 'theme', 'use')
# change to the theme chosen, here e.g. 'winnative'
# tcl('ttk::style', 'theme', 'use', 'winnative')




### popup table for info about arguments ? WORKING PARTIALLY

# displayInTable <- function(tclarray,title="",height=-1,width=-1,nrow=-1,ncol=-1)
# {
  # require(tcltk)
  # tt <- tktoplevel()
  # title <- "Info"
  # tclRequire("Tktable")
  # tkwm.title(tt,title)
  # table1 <- tkwidget(tt,"table",rows=nrow,cols=ncol,titlerows=1,titlecols=0,
                     # height=height+1,width=width+1,colwidth=80,
                     # xscrollcommand=function(...) tkset(xscr,...),yscrollcommand=function(...) tkset(yscr,...))
  # xscr <-tkscrollbar(tt,orient="horizontal", command=function(...)tkxview(table1,...))
  # yscr <- tkscrollbar(tt,command=function(...)tkyview(table1,...))

  # tkgrid(table1,yscr)
  # tkgrid.configure(yscr,sticky="nsw")
  # tkgrid(xscr,sticky="new")
  # tkconfigure(table1,variable=tclarray,background="white",selectmode="extended")
  # return (table1)
# }


## Define a matrix :

# c1 <- c("Variable","\"Law index\"","Sample size","Parameter")
# c2 <- c('Detail','law index as given by the function getindex(), should take value from 0 to 35','number of observations to generate','parameter(s) for the law')
# matrix.info <- cbind(c1,c2)

## Define a Tcl array and initialize it to that matrix :
# require(tcltk)
# tclArray1 <- tclArray()
# for (i in (1:dim(matrix.info)[1]))
  # for (j in (1:dim(matrix.info)[2]))
    # tclArray1[[i-1,j-1]] <- matrix.info[i,j]
 
	
# table1 <- displayInTable(tclArray1,nrow=dim(matrix.info)[1],ncol=dim(matrix.info)[2])


