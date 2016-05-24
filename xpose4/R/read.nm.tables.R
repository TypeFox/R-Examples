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

read.nm.tables <-
  function(table.files=NULL,
           runno=NULL,
           tab.suffix="",
           ##sim.suffix="sim",
           table.names=c("sdtab","mutab","patab","catab",
             "cotab","mytab","extra","xptab"),
           cwres.name=c("cwtab"),
           cwres.suffix="",
           quiet=FALSE,...) {
    
    if (is.null(table.files)){
      if(is.null(runno)) {
        cat(paste("runno must be specified if no table files provided\n"))
        return(NULL)
      }
      match.pos <- match(cwres.name,table.names)
      if (!is.na(match.pos)) table.names <- table.names[-match.pos]
      tab.files <- sapply(table.names,paste,runno,tab.suffix,sep="")
      cwres.files <- sapply(cwres.name,paste,runno,cwres.suffix,tab.suffix,sep="")
      tab.files <- c(tab.files,cwres.files)
    } else {
      tab.files <- table.files
    }
    
    ## Read in the table files
    totab      <- NULL
    totnam     <- NULL
    seen.files <- NULL
    filedim    <- NULL
    
    for(i in 1:length(tab.files)) {
      filename <- tab.files[i]
      if(!is.readable.file(filename)) {
        ##if (!quiet) {cat(filename,"not readable\n")}
        next
        
      } else {
        cat(paste("    Reading",filename,"\n"))

        ## Check which type of separator we have in our tables
        header.line = scan(file=filename,nlines=1,skip=1,what="character",sep="\n",quiet=T)
        sep.char = ""
        if(length(grep(",",header.line))!=0) sep.char = ","
       
        ## Check if we have unequal number of fields in the file
        ## used for multiple simulations
        fields.per.line <- count.fields(filename)
        fields.in.first.line <- fields.per.line[1]
        fields.in.rest <- fields.per.line[-1]
        if((length(unique(fields.in.rest))!=1) ||
           (all(fields.in.first.line==fields.in.rest))){ 
          if(!quiet) {
            cat(paste("Found different number of fields in ",filename,".\n",sep=""))
            cat("This may be due to multiple TABLE and header rows \n")
            cat("caused by running multiple simulations in NONMEM (NSIM > 1).\n")
            cat("Will try to remove these rows. It may take a while...\n")
          }
          tmp   <- readLines(filename, n = -1)
          inds  <- grep("TABLE",tmp)
          if (length(inds)!=1){
            inds  <- inds[c(2:length(inds))]
            inds2 <- inds+1
            tempfile<- paste(filename,".xptmp",sep="")
            write.table(tmp[-c(inds,inds2)],file=tempfile,
                        row.names=FALSE,quote=FALSE)
            assign(paste("n.",filename,sep=""),read.table(tempfile,skip=2,header=T,sep=sep.char))
            unlink(tempfile)
          } else {
            assign(paste("n.",filename,sep=""),read.table(filename,skip=1,header=T,sep=sep.char))
          }
        } else {
          assign(paste("n.",filename,sep=""),read.table(filename,skip=1,header=T,sep=sep.char))
        }
        
        ## Remember the files seen
        ##if(is.null(seen.files)) {
        ##  seen.files <- paste("n.",filename,sep="")
        ##} else {
        seen.files <- c(seen.files,paste("n.",filename,sep=""))
        ##}
      }
    }

    ## Check if we found any table files

    if(any(is.null(seen.files))) {
                                        #if(tab.suffix!=sim.suffix) {
      cat("Couldn't find any table files that match run number",
          runno, "!\n")
      return(NULL)
                                        #} else {
                                        #  cat("Couldn't find any simulation table files that match run number",
                                        #     runno, "!\n")
                                        #}
    }
    
    ## Check if the files have the same length
    for(nfile in seen.files) {
      if(is.null(filedim)) {
        filedim <- nrow(get(nfile))
      } else {
        filedim <- c(filedim,nrow(get(nfile)))
      }
    }
    
    file.df <- data.frame(seen.files=seen.files,filedim=filedim)
    lngths  <- sort(unique(file.df$filedim))

    if(length(lngths) !=1) {
      cat("\nThe table files associated with this run number (",runno,
          ") appear\n")
      cat("to have different lengths.\n")
      cat("You will have to sort this out and try again!\n")
      return(NULL)
    }


    ## Add the tables to totab and replicate the shorter ones to match
    ## the size of the longest one
    maxlngth <- max(file.df$filedim)
    
    ##singdef <-
    ##  c("id","idlab","idv","dv","pred","ipred","iwres","wres","res")
    
    for(ii in 1:nrow(file.df)) {
      filnam <- as.character(file.df[ii,"seen.files"])
      new.df <- get(filnam)
      sz     <- file.df[ii,"filedim"]
      rl     <- maxlngth/sz

      if(any(is.null(totab))) {
        totab <- new.df
      } else {
        totab <- cbind(totab,new.df)
      }
      
      totnam <- c(totnam,names(new.df))
      
      ## store parameters & covariates for Data.R & SData.R
      
##       if(!is.na(pmatch("n.patab", filnam))){ 
##         write(names(new.df), file=".patab.names.tmp")
##       } else {
##         if(!is.na(pmatch("n.catab", filnam))){ 
##           write(names(new.df), file=".catab.names.tmp")
##         } else {
##           if(!is.na(pmatch("n.cotab", filnam))){ 
##             write(names(new.df), file=".cotab.names.tmp")
##           } 
##         }
##       }

      if(!is.na(pmatch("n.patab", filnam))){ 
        write(names(new.df), file=".patab.names.tmp")
      } else {
        if(!is.na(pmatch("n.catab", filnam))){ 
          write(names(new.df), file=".catab.names.tmp")
        } else {
          if(!is.na(pmatch("n.cotab", filnam))){ 
            write(names(new.df), file=".cotab.names.tmp")
          } else {
            if(!is.na(pmatch("n.sdtab", filnam))){ 
              write(names(new.df), file=".sdtab.names.tmp")
            }  
          } 
        }
      }
      
      
    } 

                                        # cat(totnam, "\n")
    
    ## Delete all duplicates

    totab <- totab[, !duplicated(totnam)]
    return(totab)
  }

