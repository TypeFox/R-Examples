
#Function to assist in reading in data and calling PRIM

sd.start <- function(...){

#... : arguments to be passed onto sdprim, if it is called

cat("\n","\n")
cat("Welcome to the Scenario Discovery toolkit","\n","\n")

cat("Copyright (C) 2009  Evolving Logic","\n")
nicecat("This program comes with ABSOLUTELY NO WARRANTY; for details type 'sdwarranty()' when you are back at the R command prompt.  This is free software, and you are welcome to redistribute it under certain conditions; when you are not in this dialogue, type 'RShowDoc('COPYING')' for the complete GNU General Public License, which states these conditions.")

cat("\n")
cat("Note, unless you are asked specifically, you do not need to include quotes","\n", 
"around your responses to dialogues, including filenames and directories.", "\n", 
"(Nor do you need them when the dialogue uses single quotes to tell you possible
 responses - as in 'y' or 'n')","\n","\n",
"However, directory specifications should use forward slashes or double", "\n", 
"backslashes, and all dialogues require you to press ENTER to end the dialogue.","\n")
cat("\n","\n")


#Current dialogue:
cat("====== SETTING THE WORKING DIRECTORY ======","\n")
cat("The current working directory is:", "\n")
cat(getwd())
cat("\n","\n")

dirnotset <- TRUE

while(dirnotset){
  
  wkdir <- readline(cat("Please enter the directory for you analysis, or enter 'a' to accept","\n","the current directory.","\n")) 
  
  cat("\n")
  
  if (wkdir != "a"){
  
    tryout <-  try(setwd(wkdir))
    if(class(tryout)=="try-error"){
      cat("Sorry, can't seem to change to the directory you entered.  Please try again.")
      cat("\n")
      cat("\n")
    } else {
      dirnotset <- FALSE 
      cat("Working directory changed to","\n",getwd(),"\n","\n")
    }
  } else {
    dirnotset <- FALSE
    cat("Directory remains the same as above.","\n","\n")
  }

}

cat("If you ever need to change the directory later, the R command is:","\n",
"setwd(\"the new path\") (using quotes around the path)","\n","\n")

cat("====== READING IN YOUR DATA ======","\n")

nicecat("Future versions will give a choice between csv, rda files and mySQL databases.")
cat("\n")  
nicecat("mySQL is not yet implemented, so for now, please enter the name of the csv or rda file 
containing your data (no quotes, but include the .csv or .rda extension).")
cat("\n")
filename <- readline(nicecat("Note, it is assumed that your csv file begins with
a row of column names.  If that is not the case, enter 'doh' instead of the filename."))

if(filename!="doh"){
  
  nc <- nchar(filename)
  lastthree <- substr(filename,(nc-2),nc)
  
  if(lastthree=="rda"){
  
    gotloaded <- load(filename)
    mydata <- get(gotloaded)
    rm(list=gotloaded) #this takes out the original dataset, not the character 'gotloaded'
                       #would be better served by a "rename" function, though can't seem
                       #to find one.
      
  } else if(lastthree=="csv"){
  
    mydata <- read.csv(filename)
  
  }
  
} else if (filename=="doh"){

  cat("\n")
  cat("If your csv file does not start with a row of column names, you can accept","\n",
  "'bland' names from R (such as 'V1', 'V2', 'V3', etc, or you can exit this","\n",
  "program and modify your input file to include column names, or you can","\n", 
  "specify them by hand inside this program.","\n")
  
  howtofix <- readline(cat("Which would you like to do?  Enter 'bland', 'exit', or 'specify'.","\n"))
  
  if (howtofix=="exit"){
    
    stop("You chose to exit and fix your column names.  Rerun this program when done.")
  
  } else if (howtofix=="bland"){
  
    filename <- readline(cat("You chose to accept bland names given by R.  Please enter your csv filename:","\n"))
    mydata <- read.csv(filename,header=FALSE)
  
  } else if(howtofix=="specify"){
  
    filename <- readline(cat("You chose to specify your column names.  Please enter your csv filename first:","\n"))
    mydata <- read.csv(filename,header=FALSE)

    satisfied <- FALSE

    while(!satisfied){
  
      vars <- readline(cat("Your data appears to contain",ncol(mydata),"variables (columns).","\n",
      "Please enter that many variable names separated by commas.  For example:","\n",  "Myvar1,Coolvar2,spiffy3,regret ","\n"))
      
      tcnames <- strsplit(vars,",")[[1]]
  
      ok <- readline(cat("You entered the following column names:","\n",tcnames,"\n","Is this ok? ('y' or 'n')","\n"))
      
      if(ok=='n'){
        cat("You're not satisfied, try re-entering again.","\n")
        satisfied <- FALSE
      }
        
      if(ok=='y'){
        cat("Great.  Let's move on.","\n")  
        satisfied <- TRUE
      }
    }
  }
}



cat("\n")
cat("You chose the following file:","\n")
cat(getwd(),"/",filename,sep="","\n","\n")

cat("As read in by R, your data has",nrow(mydata),"rows (cases) and",ncol(mydata),"columns (variables).")
cat("\n")
cat("\n")

    nas <- is.na(mydata)
    nans <- apply(mydata,2,is.nan)
    infs <- apply(mydata,2,is.infinite)
    chars <- any(apply(mydata,2,is.character))

    if(chars){
      stop(cat("ERROR!
You have character values in your input file, which at present are
too difficult for R and this program to fix or handle appropriately.  This 
program will now exit so you can inspect your data and translate or remove 
character values.  If they are character representations of binary or ORDINAL 
variables (eg \"win\" and \"lose\" or \"low\", \"med\",\"high\"), this program
should handle them appropriately provided they are converted to numbers prior 
to reading them in.  If they represent unordered categorical variables, this 
program is not yet capable of handling that.  If you need to rely heavily on 
true categorical variables, please consider using the SuperGEM implementation
of PRIM (available free online).","\n"))
    }

if(any(c(nas,nans,infs))){
  cat("WARNING: The data contains non-finite values, specifically some mix of","\n",
  " 'NA', 'NaN', or +/- 'Inf'.","\n","These could quite likely break the algorithms.","\n")
  inspect <- readline(cat("To find out more, enter 'inspect', otherwise enter 'fine' to take your chances","\n"))
  
  cat("\n")
  
  if(inspect=="inspect"){
                  
    if(any(nas)){
      cat("There are",sum(nas),"'NA' values (including 'NaN'), distributed throughout the 
      following columns:",c(1:ncol(mydata))[apply(nas,2,any)])
      cat("\n")
    }
    
    if(any(nans)){
     cat("There are",sum(nans),"'NaN' values, distributed throughout the following columns:",c(1:ncol(mydata))[apply(nans,2,any)])
      cat("\n")
    }
    
    if(any(infs)){
      cat("There are",sum(infs),"infinite values, distributed throughout the following columns:",c(1:ncol(mydata))[apply(infs,2,any)])
      cat("\n")
    }
    

    cat("\n")    
    cat("Chances are these values will mess up the scenario discovery tools, ","\n",
    "but you have three options:  
   
    Ignore them and see what happens, 
    drop rows containing weird values from the dataset and see what happens, 
    drop variables (columns) containing weird values, or
    exit this program and fix your data however you choose.","\n","\n")
    
    whattodo <- readline(cat("Which would you like? Enter 'ignore', 'droprows', 'dropvars', or 'fix'.","\n"))
    
    anybad <- !((!nas)*(!infs))
    
    if(whattodo=="fix"){
    
      stop("You chose to fix it.  Restart this program when you are done.")
    
    } else if(whattodo=="droprows"){
    
      cat("You chose to drop records that have weird data.","\n")
      badrows <- apply(anybad,1,any)
      dropok <- readline(cat("This will drop",sum(badrows),"out of",nrow(mydata),"rows.  Is that ok?  (Enter 'y' or 'n')","\n"))
      
      if(dropok=='y'){
        
        mydata <- mydata[!badrows,]
      
      } else if(dropok=='n'){
      
        cat("\n")
        stop(cat("You have bad data and you said dropping that data was not ok.  
        Please rerun this program after you fixed your data."))
    
      } 
    
    } else if(whattodo=="dropvars"){
    
      badcols <- apply(anybad,2,any)
    
      cat("You chose to drop variables containing weird entries.","\n")
      dropok <- readline(cat("This will drop the following columns out of",ncol(mydata),"columns total:",
      c(1:ncol(mydata))[badcols],"\n", 
      "Is that ok?  (Enter 'y' or 'n')","\n"))
      
      if(dropok=='y'){
        
        mydata <- mydata[,!badcols]
      
      } else if(dropok=='n'){
      
        stop(cat("You have bad data and you said dropping that data was not ok. 
         Please rerun this program after you fixed your data."))
    
      } 
      
    } else if(whattodo=="ignore"){
    
      cat("You chose to keep the data as is and see what happens.")
      
    }
        
    cat("\n")
    cat("After making any modifications describe above, your data now 
    contains",nrow(mydata),"rows and",ncol(mydata),"columns.")
    cat("\n")
    cat("\n")
   
  }
}



cat("The variable names (and column number) are as follows:","\n")
print(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata))," name")))

cat("\n")

datagood <- readline(cat("Are the size and variables as expected?  If so enter 'y', if not enter 'n'. 
(Note that you will be able to specify a subset of variables you'd like to use
in the analysis later.)","\n"))

if(datagood=="n"){cat("\n","If the data does not look as you expected it to, there are a few possible 
causes, though 'NA' related causes should have been taken care of above.
  
If you have one fewer row than expected, you maybe lied when you implied that 
your data had column names in the first row.  Inspect your csv file to see if 
that is the case.  If so, either name your columns and start this program over,
or keep them unnamed and start the program over and enter 'doh' at the
appropriate place in dialogue, and R will give your data bland names 
(like \"V1\",\"V2\",\"V10\").")

keepgoing <- readline(cat("To continue anyway, enter 'y', to stop enter 'n'.","\n"))

if(keepgoing=="n"){
  stop("Alright, see you next time.")
}


}   


#******************************************************************
## CODE TO AUTOMATICALLY REMOVE COLUMNS CONTAINING ONLY 1 VALUE:

x <- mydata
xfull <- x   #preserve original x
repcols <- logical(length=ncol(x)) #make vector indicating which have duplicates

#identify duplicates
for (i in 1:ncol(x)){
	if (isTRUE(all.equal(xfull[,i],rep(x[1,i],nrow(x))))){repcols[i] <- TRUE}
}

if(any(repcols)){
  
  cat("\n")
  cat("====  WARNING  ====","\n")
  cat("The following columns were identified as containing entirely single values:","\n")
  
  for(i in which(repcols)){
  
    cat("Column", i, "corresponding to the variable", colnames(x)[i],"\n")

  }
  cat("\n")
  cat("Because they show no variation, it is impossible for these variables to
play a role in the box definition, and they will be removed and/or ignored to
avoid fouling up the algorithms.")

   mydata <- x[,!repcols]

} else{
  cat("\n","Checked for columns with all identical entries - none found.")
}

#************************************************************
#Check for duplicate columns:

cat("\n")
cat("\n")
dupout <- dupcolchecker(mydata)

#This structure for checking the output of dupcolchecker is dumb...

if(length(dupout)==1){
  if(dupout=="nodups"){
  #everything is fine
  } 
} else if (is.vector(dupout)){
  nicecat("Here's the vector of columns that should be deleted.  These remain in 
  the dataset for now and may cause problems.")
  print(dupout)
} else {
  mydata <- dupout
}
  
  mydata <- as.matrix(mydata)

 
##*****************************

#Display variables and select input vars


varsatisfied <- FALSE

while(!varsatisfied){

cat("======  SELECTING VARIABLES FOR ANALYSIS  ======","\n","\n")
                 
satisfied <- FALSE
                  
while(!satisfied){                                                         
  cat("To remind you, your data has the following variables:","\n")
  cat("\n")
  print(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata))," name")))
  
  cat("\n")
  
  sensicalinput <- FALSE
  
  while(!sensicalinput){
  
  sensicalinput <- TRUE
  
  nicecat("Please specify which variables you would like to have the algorithms
  consider as input variables.  If you would like to specify all variables from
  column number X to column Y, you can enter 'X:Y'.  If you only want specific
  variables, you can enter them separated by commas, for example:  '1,2,3,6,7'
  (You can also mix, like '1,3,5:7'.)")
  cat("\n")
  inputvars <- readline(nicecat("Also, note that if you modified columns above, you 
  will need to specify the new column numbers since they may have changed from the original dataset."))
  
  
  cat("\n")
  
  ##RECAST INPUT VARS!!!
  
  tinputvars <- as.vector(strsplit(inputvars,",")[[1]])
  
  hascolon <- grep(":",tinputvars)
  
  inputvars <- c()
    
  for (i in 1:length(tinputvars)){
  
    if(i %in% hascolon){
      cols <- as.numeric(as.vector(strsplit(tinputvars[i],":")[[1]]))
      inputvars <- c(inputvars,c(cols[1]:cols[2])) 
    } else {
      inputvars <- c(inputvars,as.numeric(tinputvars[i])) 
    }
  }
      
  if(length(inputvars)==1){
    cat("\n","You selected only one input variable, which doesn't make sense for using PRIM.","\n",
    "Try again.","\n")
    sensicalinput <- FALSE 
  }
  
  }
  
  cat("You selected the following input variables:,","\n")
  cat("\n")
   
  print(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata)),"name"))[inputvars,])
  #CHECK THAT THIS DOES THE COL NUMBER CORRECTLY
  cat("\n")
  good <- readline(cat("Is this correct?  Enter 'y' to accept or 'n' to try again.","\n"))
  
  if (good=='y'){
    satisfied <- TRUE
  }

}

#Display remaining vars and specify output vars to consider:

satisfied <- FALSE

while(!satisfied){       

  cat("Here are the remaining variables (those not already designated as 
input variables:","\n")

  cat("\n")
  coln <- 1:ncol(mydata)
  
  print(matrix(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata)),"name"))[-inputvars,],ncol=1,dimnames=list(coln[-inputvars],"name")))

  cat("\n")

  outputvars <- readline(cat("Please specify which variables you would like to have the algorithms
consider as OUTPUT variables, using the same format as for input variables.","\n"))
  
  cat("\n")
  
  ##Translate their character input into actual vector: 
  
  toutputvars <- as.vector(strsplit(outputvars,",")[[1]])
  
  hascolon <- grep(":",toutputvars)
  
  outputvars <- c()
    
  for (i in 1:length(toutputvars)){
  
    if(i %in% hascolon){
      cols <- as.numeric(as.vector(strsplit(toutputvars[i],":")[[1]]))
      outputvars <- c(outputvars,c(cols[1]:cols[2])) 
    } else {
      outputvars <- c(outputvars,as.numeric(toutputvars[i])) 
    }
  }
    
  
  
  cat("You selected the following output variables:","\n")
   
  print(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata))," name"))[outputvars,])
  #CHECK THAT THIS DOES THE COL NUMBER CORRECTLY
  cat("\n")
  good <- readline(cat("Is this correct?  Enter 'y' to accept or 'n' to try again.","\n"))
  
  if (good=='y'){
    satisfied <- TRUE
  }

}

####### PLOT HISTS AND CDFS
#
#if(length(outputvars)>1){
#
#  cat("Would you like the histograms and cdfs for the different output variables
#  to include all output variables on each plot, or have individual plots?","\n")
#  complot <- readline(cat("(Enter 'common' or 'ind')","\n"))
#  
#}


#PLOT HISTOGRAMS/CDF's
#For all new windows:
cat("\n")
cat("====== PLOTTING HISTOGRAMS AND CDFs ======","\n")
cat("There should now be one or more graphics windows open displaying useful","\n")
cat("histograms and cdfs to help you choose a threshold for your data.","\n",
"reactivate this window to continue the dialogue.")

oldpar <- par(no.readonly = TRUE)

graphics.off()
for (i in 1:length(outputvars)){
  xname = colnames(mydata)[outputvars[i]]
  options("device")$device(width=14)
  par(mfcol=c(1,2))
  hist(mydata[,outputvars[i]],xlab=xname,main = paste("Histogram of",xname),breaks=10)
  plot.ecdf(mydata[,outputvars[i]],xlab=xname, main = paste("Empirical CDF of",xname), ylab=paste("Fn(",xname,")",sep=""), verticals=TRUE)
}

#RESTORE GRAPHICS PARAMETERS:
par(oldpar)




cat("\n")
cat("\n")
cat("====== CHOOSING THRESHOLDS ======")

threshvals    <- vector(length=length(outputvars))
threshhighlow <- vector(length=length(outputvars))

for (i in 1:length(outputvars)){
  
  threshsatisfied <- FALSE
  
  while(!threshsatisfied){ 
    cat("\n")
    uplownone <- readline(cat("For variable",colnames(mydata)[outputvars[i]],"would you like to define a:
    1. Upper Threshold (values above are 1, values below are 0)
    2. Lower Threshold (values below are 1, values above are 0)
    3. No Threshold
    (Enter '1', '2' or '3')","\n"))
    
    if(!(uplownone=="1" | uplownone=="2" | uplownone=="3")){
      cat("Invalid entry, please try again.","\n")
    } else {  
    
    if(uplownone=="3"){
      cat("You chose not to define a threshold.  Moving on...","\n")
      threshhighlow[i] <- "3"
      threshsatisfied <- TRUE
    } else {
    if(uplownone=="1"){
      ulnchar <- "upper"
      upthresh <- readline(cat("You chose to define an upper threshold.  Please enter it now:","\n"))
      
      upthresh <- as.numeric(upthresh)
      at <- sum(mydata[,outputvars[i]]>upthresh)
      tot <- nrow(mydata)
      cat("Using",upthresh,"as a threshold, you will have:","\n",
      at,"out of",tot,"above the threshold (",100*at/tot,"percent)","\n",
      (tot-at),"out of",tot,"below the threshold (",100*(1-at/tot),"percent)","\n")
      cat("\n")
  
    }
    
    if(uplownone=="2"){
      ulnchar <- "lower"
      lthresh <- readline(cat("You chose to define a lower threshold.  Please enter it now:","\n"))
      
      lthresh <- as.numeric(lthresh)
      at <- sum(mydata[,outputvars[i]]<lthresh)
      tot <- nrow(mydata)
      cat("Using",lthresh,"as a lower threshold, you will have:","\n",
      at,"out of",tot,"below the threshold (",100*at/tot,"percent becoming '1's)","\n",
      (tot-at),"out of",tot,"above the threshold (",100*(1-at/tot),"percent becoming '0's)","\n")
  
      cat("\n")
  
    }
  
    cat("Are you satisfied with this threshold for this variable?","\n")
    keeporno <- readline(cat("Enter 'y' or 'n'","\n"))
    
    if(keeporno=="y"){
      if(uplownone=="1"){
        threshvals[i]   <- upthresh
        threshhighlow[i] <- uplownone
      } else{
        threshvals[i]   <- lthresh
        threshhighlow[i] <- uplownone
      }
      
      threshsatisfied <- TRUE
    } else{
    cat("\n")
    nicecat("You said you were unsatisfied - try again, 
    or move on by selecting option 3 below.")
    }
    }
    }
  }
  
}

infinalq <- TRUE

while(infinalq){

  cat("\n")
  cat("Would you like to:
  
  1) save your cleaned data so you can call sdprim directly in the future
  2) begin the PRIM analysis using 'sdprim', calling it from this program? 
  3) redefine variables?
  4) start over from the very beginning, possibly with a new dataset?
  5) exit this dialogue, and possibly begin the PRIM analysis on your own using
  the 'sdprim' command?  (Recommended if you wish to change some parameters 
  from their defaults.)
  ","\n")
  whattodo <- readline(cat("(Enter '1', '2','3','4' or '5')","\n"))
  cat("\n")
  
  infinalq <- FALSE
  varsatisfied <- TRUE
  
  if (whattodo =="1"){
  
    inputneeded <- TRUE
    infinalq <- TRUE
  
    while(inputneeded){
    
      nicecat("If your dataset was modified above,
      you might like to save a copy of the modified version so that in the future you can
      run sdprim directly, without needing to 'fix' your data. 
      If you would like to save a copy, either a csv or rda file, then enter the filename you would like to save it as, including
      the extension, and it will automatically save in appropriate format.  Otherwise enter 'n'.")
      cat("\n")
      
          
      savetype <- readline(nicecat("(Note that this will not save the thresholded data as an addition to your dataset - rather, it will simply save the original
      dataset, but without the rows and columns and entries that caused problems.)"))
      
      nc <- nchar(savetype)
      
      if(nc > 4){ #They in theory specified a filename
        
        ext <- substr(savetype, nc-3,nc)
        
        if(ext==".csv" | ext==".rda"){
          if(ext==".csv"){
            write.csv(mydata, file=savetype, row.names=FALSE)
          } else { 
            save(mydata, file=savetype)
          }
          inputneeded <- FALSE
          cat("\n")
          cat("Modified data saved as file", savetype,"\n")
          cat("\n")
         } else {
          nicecat("Sorry, you input a filename that is either missing an extension, or 
          does not have an extension that matches the currently allowed forms of either rda or csv.
          You will now be returned to the start of the data-saving dialogue.")
          cat("\n")
          cat("\n")
        } 
      
      } else {  #They did not want to save
      
        inputneeded <- FALSE
        nicecat("You decided not to save.")
        
      }
    
    }
  
  }
  
  else if (whattodo =="2"){
  
    cat("\n")
  
  
  
    cat("\n")
    cat(" ====== Running PRIM Analysis ====== ","\n")
    
    if(length(outputvars)!=1){
        
      cat("\n")
      cat("Which output variable would you like?","\n")
      cat("(If you defined a threshold it will be used, if not the output variable
      should already be binary.)","\n")
      
      coln <- 1:ncol(mydata)
  
      print(matrix(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata)),"name"))[outputvars,],ncol=1,dimnames=list(coln[outputvars],"name")))
      cat("\n")
      outputvar <- as.numeric(readline(cat("Enter the number next to the variable you would like to use","\n")))
    }else{
      cat("You chose to run PRIM and you only had only identified one output variable:","\n")
      print(colnames(mydata)[outputvars])
      cat("\n")
      outputvar <- outputvars
    }                               
    
    vec <- mydata[,outputvar]
    
    #check for being binary:
    looksbinary <- all(vec==0 | vec==1)
    
    theind <- which(outputvars==outputvar)
    
    if(!(threshhighlow[theind]=="1" | threshhighlow[theind]=="2") & !looksbinary){
      cat("You did not previously select a threshold for this variable.","\n")
      
      cat("Also it does not appear to be binary. Thus you will need to start over,","\n",
        "at the variable selection phase, which you will be returned to now.")
        infinalq <- FALSE
        varsatisfied <- FALSE
       
    } else {
    
    
      if(threshhighlow[theind]=="3"){
      
        cat("You did not previously select a threshold for this variable.","\n")
    
        cat("Fortunately, it appears to already be binary, so it will get used as is.","\n")
        yvec <- mydata[,outputvar]
      
      } else if(threshhighlow[theind]=="1"){
        cat("Earlier, you defined an upper threshold of",threshvals[theind],"\n",
        "This will be used for the PRIM analysis.","\n")
        
        yvec <- 1*mydata[,outputvar]>threshvals[theind]
        
      } else if(threshhighlow[theind]=="2"){
      
         cat("Earlier, you defined a lower threshold of",threshvals[theind],"\n",
        "This will be used for the PRIM analysis.","\n")
        
        yvec <- 1*mydata[,outputvar]<threshvals[theind]
        
      } 
    
      cat("\n")
      cat("\n")
      cat("PRIM will be using the following input variables:","\n")
      cat("\n")
      #PRINT INPUT VARIABLES
      print(matrix(colnames(mydata),ncol=1, dimnames=list(c(1:ncol(mydata)),"name"))[inputvars,])
      
      cat("\n")
      cat("And the following output variable:  ") 
      print(colnames(mydata)[outputvar])
      cat("\n")
      
      #CHECK FOR NEW ARGUMENTS:
      cat("\n")
      
      cat("The default arguments for sdprim are listed below:","\n")
      c("\n")
      sdpformals <- formals(sdprim)
      
      #note, using 3 skips the problematic "x" and "x[,ncol(x)]" arguments
      
      for (i in 3:length(sdpformals)){
        if(!identical(sdpformals[[i]],"")){
          if(is.character(sdpformals[[i]])){
            cat("  ",names(sdpformals)[i],"=");cat(c(),sdpformals[[i]],c(),sep="\"","\n")
          } else if(is.null(sdpformals[[i]])){
            cat("  ",names(sdpformals)[i],"=","NULL","\n")
          } else if(sdpformals[[i]]==""){
            #Don't print
          } else{
            cat("  ",names(sdpformals)[i],"=",sdpformals[[i]],"\n")
          }
        }
      }
      
      cat("\n")
      
      sdscall <- match.call()
      
      if(length(match.call())>1){
        cat("Based on the arguments passed to sd.start, these will be modified as follows:")
        cat("\n")
#        sdscall <- match.call()
      
        for (i in 2:length(sdscall)){
          cat("  ",names(sdscall)[i],"=",sdscall[[i]],"\n")
        }
      }
      cat("\n")#
      
  
      
      mod <- readline(cat("Would you like to modify them further? (enter 'y' or 'n')","\n"))
      
      if(mod=="y"){
          nicecat("Please enter the arguments you would like to change, specifying the argument name
          followed by an = sign, followed by the value of the argument, with each argument separated by a comma.")
          cat("For example, you might enter:","\n")
          cat("peel.alpha=.2, outfile=\"somethingnew.txt\"","\n")
          newargs <- readline()                                    


          eacharg <- strsplit(newargs,",")
          argvals <- strsplit(eacharg[[1]],"=")

          #remove spaces and escape characters          
          backtogether <- argvals
          
          for (i in 1:length(argvals)){
            
            for (j in 1:length(argvals[[i]])){
            
              allchars <- strsplit(argvals[[i]][j],"")[[1]]
              notspaces <- which(allchars!=" " & allchars!="\"")
              backtogether[[i]][j] <- paste(allchars[notspaces],collapse="")
            
            }
            
            
            currwarn <- getOption("warn")
            options(warn=-1) #suppress warnings from as.numeric(character)
            
            assig <- backtogether[[i]][2]
            
            #test whether character, logical or numeric
            if(!is.na(as.logical(assig))){#then is logical
              sdscall[[backtogether[[i]][1]]] <- as.logical(assig)
            } else if (!is.na(as.numeric(assig))){ #is numeric, not character
              sdscall[[backtogether[[i]][1]]] <- as.numeric(assig)
            } #default is to remain as character
            
            
            options(warn=currwarn)
            
          }
      
      }
          
      cat("\n")
      cat("\n")
      cat("\n")
      cat("====== NOW ENTERING SDPRIM FUNCTION ======","\n")
      cat("\n")
      
      flush.console()
      
      sdscall$x <- mydata[,inputvars]
      sdscall$y <- yvec
      sdscall[[1]] <- sdprim
#      sdprimresults <- sdprim(x=mydata[,inputvars],y=yvec, ...)
#      print(sdscall)
  
  
      sdprimresults <- eval(sdscall)
      
      return(sdprimresults)
  
    }
    
  } else if (whattodo=="3"){
    
    infinalq <- FALSE
    varsatisfied <- FALSE
    cat("\n")
    cat("You chose to redefine variables.","\n")
  
  } else if (whattodo=="4"){
  
    cat("\n")
    cat("You chose to start over from the beginning.","\n")
    cat("For now, please implement that yourself after this program exits by typing
    'sd.start()' at the command line you are about to see.  Thanks.")
    cat("\n")
    infinalq <- FALSE
    
  } else if (whattodo=="5"){
  
    cat("\n")
    cat("Have fun.","\n")
    cat("\n")
    infinalq <- FALSE
    
  } else {cat("\n")
      
      cat("Didn't recognize what you entered - please try again.","\n")
      infinalq <- TRUE
  }

}

}

}  
   
   
  #

#
#what <- readline(cat("Enter '1', '2' or '3'","\n"))
#
#if (what=="1"){
#  
#  cat("Which output variable would you like to use?")
#  
#
#
#
#  outvar <- readline(cat("Enter the number rather than the variable name","\n"))
#
#
##   
##  } else if(uplownone=="2"){
##    ulnchar <- "lower"       
##  } else if(uplownone=="3"){
##    ulnchar <- "No threshold 
#

  


####CHANGE BELOW TO INCLUDE INPUT/OUTPUT VARS
#cat("If you did any significant cleaning to your data, you might want to save it
#as '.rda' file so you can just load it up no problem next time.")
#cat("\n")
#saverda <- readline("To do this, enter the filename (including the rda extension)
#or enter 'n' to continue without saving.")
#
#if(saverda!='n'){
#  save(mydata, file=saverda)
#}
#
##when get to prim, printout default args