# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"read.npc.vpc.results" <-
  function(vpc.results=NULL,
           npc.results=NULL,
           verbose=FALSE,
           ...) {


    ## Make sure we have the necessary variables defined
    if(is.null(vpc.results) & is.null(npc.results)){
      cat(paste("Both the arguments vpc.results and npc.results are NULL\n"))
      cat(paste("One of these must be defined\n"))
      return(NULL)
    }

    ## Make sure we have the necessary variables defined
    if(!is.null(vpc.results) & !is.null(npc.results)){
      cat(paste("Both the arguments vpc.results and npc.results are defined\n"))
      cat(paste("ONLY one of these may be defined\n"))
      return(NULL)
    }

    vpc <- FALSE
    npc <- FALSE
    if(!is.null(vpc.results))  vpc <- TRUE
    if(!is.null(npc.results))  npc <- TRUE

    if(vpc)  filename <- vpc.results
    if(npc)  filename <- npc.results
    if(is.readable.file(filename)){
      if(verbose) cat(paste("    Reading",filename,"\n"))
      file.scan <- scan(filename,sep="\n",
                        what = character(),
                        quiet=TRUE,
                        blank.lines.skip = FALSE)
    } else {
      cat(paste(filename,"was not found for reading\n"))
      return(NULL)
    }

    blank.line.pat <- "^$"
    if(vpc){
      table.start.pat <- "VPC results"
      table.head.pat <- "<="
    }
    if(npc){
      table.start.pat <- "NPC results"
      table.head.pat <- "points below PI"
    }
    table.start <- grep(table.start.pat,file.scan)
    num.tables <- length(table.start)
    table.head <- grep(table.head.pat,file.scan)

    ## get end of table
    blank.line <- grep(blank.line.pat,file.scan)
    table.stop <- c()
    for(i in 1:num.tables){
      for(j in 1:length(blank.line)){
        if (table.start[i]>blank.line[j]) next
        if (table.start[i]<blank.line[j]) {
          table.stop <- c(table.stop,blank.line[j]-1)
          break
        }
      }
    }
    table.rows.to.read <- table.stop-table.head

    ## get the DV and IDV for the plot
    dv.pat <- "Dependent variable"
    idv.pat <- "Independent variable"
    mod.pat <- "Modelfile"
    dv.idv.table.start <- grep(dv.pat,file.scan)
    dv.idv.table.stop <- NULL
    for(j in 1:length(blank.line)){
      if (dv.idv.table.start>blank.line[j]) next
      if (dv.idv.table.start<blank.line[j]) {
        dv.idv.table.stop <- blank.line[j]-1
        break
      }
    }
    dv.idv.table <- read.table(filename,
                               skip=dv.idv.table.start-1,
                               nrows=dv.idv.table.stop - dv.idv.table.start,
                               sep=",",comment.char="",
                               header=T,strip.white=TRUE)

    dv.var <- paste(dv.idv.table[[grep("Dependent.variable",names(dv.idv.table))]])
    ##dv.var <- paste(dv.idv.table$Dependent.variable)
    model.file <- paste(dv.idv.table[[grep("Modelfile",names(dv.idv.table))]])
    ##model.file <- paste(dv.idv.table$Modelfile)
    if (vpc) idv.var <- paste(dv.idv.table[[grep("Independent.variable",names(dv.idv.table))]])
    ##if (vpc) idv.var <- paste(dv.idv.table$Independent.variable)

    ## get categorical boundaries if present
    cat.tables <- F
    cen.tables <- F
    cat.boundaries <- NULL
    lloq <- NA
    uloq <- NA
    pred.corr <- FALSE
    var.corr <- FALSE
    add.feats.row <- grep("Additional.feature",file.scan)
    ##if(!is.null(dv.idv.table$Additional.features)){
    if(length(add.feats.row)!=0){
      for(i in 1:length(add.feats.row)){
        ##if(dv.idv.table$Additional.features=="Categorization"){
        if(length(grep("Categorization",file.scan[add.feats.row[i]+1]))!=0){
          ##if(dv.idv.table[[add.feats.row[i]]]=="Categorization"){
          cat.tables <- T
          boundary.table <- read.table(filename,
                                       skip=add.feats.row[i]+1,
                                       nrows=1,
                                       sep=",",
                                       fill=T,
                                       comment.char="",
                                       strip.white=TRUE,
                                       header=T)

          boundary.rows <- grep("Boundary",names(boundary.table))
          cat.boundaries <- boundary.table[,boundary.rows]
        }
        if(length(grep("Censored.data",file.scan[add.feats.row[i]+1]))!=0){
          ##if(dv.idv.table[[add.feats.row]]=="Censored data"){
          ##if(dv.idv.table$Additional.features=="Censored data"){
          cen.tables <- T
          censored.table <- read.table(filename,
                                       skip=add.feats.row[i]+1,
                                       nrows=1,
                                       sep=",",
                                       fill=T,
                                       comment.char="",
                                       strip.white=TRUE,
                                       header=T)
          lloq <- censored.table$LLOQ
          uloq <- censored.table$ULOQ
          ##lloq.row <- grep("LLOQ",names(dv.idv.table))
          ##uloq.row <- grep("ULOQ",names(dv.idv.table))
          ##lloq <- dv.idv.table[,lloq.row]
          ##uloq <- dv.idv.table[,uloq.row]
        }
        if(length(grep("Prediction.correction",file.scan[add.feats.row[i]+1]))!=0){
          ##if(dv.idv.table[[add.feats.row[i]]]=="Categorization"){
          pred.corr <- T
        }
        if(length(grep("Variability.correction",file.scan[add.feats.row[i]+1]))!=0){
          ##if(dv.idv.table[[add.feats.row[i]]]=="Categorization"){
          var.corr <- T
        }
      }
    }


    ##get the numerous tables
    by.interval <- NULL
    strata.names <- NULL
    if(num.tables>1){
      bin.table <- vector("list",num.tables+1)
      strata.names <- c()
      tmp.interval <- c()
      for(i in 1:num.tables){
        ## get the names of the strata
        strata.pat <- "strata"
        strata.line <- file.scan[table.start[i]]
        strata.start <- regexpr(strata.pat,strata.line)+7
        if(strata.start==6){ ## there is no strata
          tmp.strata <- NULL
          tmp.interval <- NULL
        } else {
          strata.stop <- regexpr(",",strata.line)-1
          if (strata.stop==-2) strata.stop <- regexpr("$",strata.line)
          tmp.strata <- substring(strata.line,strata.start,strata.stop)
          strata.var.stop <- regexpr(" ",tmp.strata)-1
          strata.int.start <- regexpr(" ",tmp.strata)+1
          strata.var <- substring(tmp.strata,1,strata.var.stop)
          strata.int <- substring(tmp.strata,strata.int.start)
          tmp.strata = gsub("\\[",">= ",tmp.strata)
          tmp.strata = gsub("\\]",paste(" >=",strata.var,sep=" "),tmp.strata)
          tmp.strata = gsub("\\(","> ",tmp.strata)
          tmp.strata = gsub("\\)",paste(" >",strata.var,sep=" "),tmp.strata)
          tmp.strata = gsub("\\;"," \\& ",tmp.strata)
          tmp.strata = gsub(" = ", " == ",tmp.strata)
        }
        strata.names <- c(strata.names,tmp.strata)

        if(!is.null(tmp.strata)){
          if(regexpr(">",tmp.strata)!=-1){ #then we have a continuous variable
            ## set the intervals for the continuous variable
            semi.loc <- regexpr("\\;",strata.int)
            lt.GE.loc <- regexpr("\\[",strata.int)
            lt.GT.loc <- regexpr("\\(",strata.int)
            rt.LE.loc <- regexpr("\\]",strata.int)
            rt.LT.loc <- regexpr("\\)",strata.int)

            strata.int.low <- substring(strata.int,1,semi.loc-1)
            strata.int.low <- gsub("\\[","",strata.int.low)
            strata.int.low <- gsub("\\(","",strata.int.low)
            strata.int.low <- gsub(" ","",strata.int.low)
            strata.int.low <- as.numeric(strata.int.low)

            strata.int.high <- substring(strata.int,semi.loc+1)
            strata.int.high <- gsub("\\]","",strata.int.high)
            strata.int.high <- gsub("\\)","",strata.int.high)
            strata.int.high <- gsub(" ","",strata.int.high)
            strata.int.high <- as.numeric(strata.int.high)

            interval.length <- strata.int.high - strata.int.low
            add.to.ends <- interval.length*0.0000001
            if(lt.GT.loc!=-1) strata.int.low <- strata.int.low+add.to.ends
            if(rt.LT.loc!=-1) strata.int.high <- strata.int.high-add.to.ends

            tmp.interval <- c(tmp.interval,strata.int.low,strata.int.high)
          }
        }

        ## get the table
        bin.table[[i]] <- read.table(filename,
                                     skip=table.head[i]-1,
                                     nrows=table.rows.to.read[i],
                                     sep=",",comment.char="",
                                     header=T,strip.white=TRUE,
                                     blank.lines.skip=FALSE)

      }
      if(length(tmp.interval)!=0){
        by.interval <- matrix(tmp.interval,
                              nrow=num.tables,
                              ncol=2,byrow=T)
      }
      bin.table[[num.tables+1]] <- strata.names
    } else {
      bin.table <- read.table(filename,
                              skip=table.head-1,
                              nrows=table.rows.to.read,
                              sep=",",comment.char="",
                              header=T,strip.white=TRUE,
                              blank.lines.skip=FALSE)
    }

    ## rename headers for tables
    for(i in 1:num.tables){
      if (num.tables==1){
        tmp.table <- bin.table
      } else {
        tmp.table <- bin.table[[i]]
      }

      if(vpc){
        tmp.table$X <- NULL
        tmp.table$X.1 <- NULL

        names(tmp.table)[1] <- "lower"
        names(tmp.table)[2] <- "upper"
        names(tmp.table)[3] <- "nobs"
      }

      if(npc){
        names(tmp.table)[1] <- "PI"
        tmp.table$PI <- as.numeric(sub("% PI","",tmp.table$PI))
      }

      tmp.names <- names(tmp.table)
      tmp.names <- sub("X\\.*","",tmp.names)
      tmp.names <- gsub("_",".",tmp.names)
      tmp.names <- gsub("\\.+","\\.",tmp.names)
      tmp.names <- gsub("\\.$","",tmp.names)

      names(tmp.table) <- tmp.names

      ##browser()
      ##names(bin.table)
      ##names(tmp.table)

      if (num.tables==1){
        bin.table <- tmp.table
      } else {
        bin.table[[i]] <- tmp.table
      }
    }

    ## make a categorical, censored and continuous list if needed
    table.multiples = 1
    if(cat.tables) table.multiples=table.multiples+1
    if(cen.tables) table.multiples=table.multiples+1
    if(table.multiples > 1){
      bin.table.cont <- vector("list",num.tables/table.multiples)
      if(cat.tables) bin.table.cat <- vector("list",num.tables/table.multiples)
      if(cen.tables) bin.table.cen <- vector("list",num.tables/table.multiples)

      sub.i <- 0
      for(ii in seq(1,num.tables,by=table.multiples)){
        sub.i <- sub.i+1
        bin.table.cont[[sub.i]] <- bin.table[[ii]]
      }
      if(sub.i==1) bin.table.cont <- bin.table.cont[[sub.i]]

      cen.start = 2
      if(table.multiples==3){
        cat.start = 3
      } else {
        cat.start = 2
      }
      if(cen.tables){
        sub.i <- 0
        for(ii in seq(cen.start,num.tables,by=table.multiples)){
          sub.i <- sub.i+1
          bin.table.cen[[sub.i]] <- bin.table[[ii]]
        }
        if(sub.i==1) bin.table.cen <- bin.table.cen[[sub.i]]
      } else {
        bin.table.cen <- NULL
      }
      if(cat.tables){
        sub.i <- 0
        for(ii in seq(cat.start,num.tables,by=table.multiples)){
          sub.i <- sub.i+1
          bin.table.cat[[sub.i]] <- bin.table[[ii]]
        }
        if(sub.i==1) bin.table.cat <- bin.table.cat[[sub.i]]
      } else {
        bin.table.cat <- NULL
      }

      if(cat.tables){
        num.tables.cat <- num.tables/table.multiples
      } else {
        num.tables.cat <- NULL
      }
      if(cen.tables){
        num.tables.cen <- num.tables/table.multiples
      } else {
        num.tables.cen <- NULL
      }
      num.tables.cont <- num.tables/table.multiples
      strata.names <- strata.names[seq(1,num.tables,by=table.multiples)]

    } else {
      bin.table.cont <- bin.table
      bin.table.cat <- NULL
      bin.table.cen <- NULL
      num.tables.cont <- num.tables
      num.tables.cat <- NULL
      num.tables.cen <- NULL
    }

    if(npc) return(list(model.file=model.file,
                        dv.var=dv.var,
                        idv.var=NULL,
                        num.tables=num.tables,
                        result.tables=bin.table))

    if(vpc) return(list(model.file=model.file,
                        dv.var=dv.var,
                        idv.var=idv.var,
                        num.tables=num.tables.cont,
                        by.interval=by.interval,
                        result.tables=bin.table.cont,
                        strata.names=strata.names,
                        num.tables.cat=num.tables.cat,
                        result.tables.cat=bin.table.cat,
                        cat.boundaries=cat.boundaries,
                        num.tables.cen=num.tables.cen,
                        result.tables.cen=bin.table.cen,
                        lloq=lloq,uloq=uloq,
                        pred.corr=pred.corr,
                        var.corr=var.corr
                        ))
  }
