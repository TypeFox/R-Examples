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

sqrtm <- function(x) {
  xe <- eigen(x)
  xe1 <- xe$values
  if(all(xe1 >= 0)) {
    xev1 <- diag(sqrt(xe1),nrow=length(xe1))
  } else {
    i=1
    while(i<(length(xe1)+1)){
      if(xe1[i]<0){
        xe1[i]=0
      }
      i=i+1
    }
    xev1 <- diag(sqrt(xe1),nrow=length(xe1))
  }
  xval1 <- cbind(xe$vectors)
  xval1i <- solve(xval1)
  y <- xval1 %*% xev1 %*% xval1i

  ## test with
  ##foo <- matrix(c(5,-4,1,0,0,-4,6,-4,1,0,1,-4,6,-4,1,0,1,-4,6,-4,0,0,1,-4,5),nrow=5,ncol=5)
  ##foo1 <- sqrtm(foo)
  ##foo1
  ##foo1%*%foo1
  ##

}

"is.cwres.readable.file"  <- function(filename )
{
  ## If we are not dealing with R -> Splus
  if(is.null(version$language)) {
    cat("This version of Xpose uses R")
    ## if(platform() == "WIN386") {
    ##   access(filename, 4) == 0
    ## } else {
    ##   filename <- paste("'", filename, "'", sep = "")
    ##   sapply(paste("test -f", filename, "-a -r", filename), unix,
    ##          output = F) == 0
    ## }
  } else {
    return(file.exists(filename)[1])
  }    
}


read.cwres.data <-
  function(filename,
           old.file.convention=FALSE,
           est.tab.suffix=".est",
           deriv.tab.suffix=".deriv",
           ...) {
    
    tables.read <- FALSE
    
    if(old.file.convention){
      if ((is.cwres.readable.file(paste(filename,".50",sep=""))) &&
          ##(is.cwres.readable.file(paste(filename,".51",sep=""))) &&
          (is.cwres.readable.file(paste(filename,".52",sep=""))) &&
          ##(is.cwres.readable.file(paste(filename,".53",sep=""))) &&
          (is.cwres.readable.file(paste(filename,".54",sep=""))) &&
          ##(is.cwres.readable.file(paste(filename,".55",sep=""))) &&
          (is.cwres.readable.file(paste(filename,".56",sep=""))) &&
          ##(is.cwres.readable.file(paste(filename,".57",sep=""))) &&
          (is.cwres.readable.file(filename))) {

        nsim <- 1
        
        num.fields <- count.fields(filename,skip=1)
        if((length(unique(num.fields))!=1) ||
           ## These lines fix problem if num_fields matches fields first header row 
           (unique(num.fields) == 8) ||
           (unique(num.fields) == 3)) { 
          
          tmp   <- readLines(filename, n = -1)
          inds  <- grep("TABLE",tmp)
          if (length(inds)!=1){
            cat("Multiple simulations not supported\n")
            cat("using this old file convention\n")
            return(NULL)
          } else {
            data <- read.table(filename,skip=1,header=T)
          }
        } else {
          data <- read.table(filename,skip=1,header=T)
        }
        size.of.sim <- dim(data)[1]/nsim
        data[,"iteration.number"] <- sort(rep(1:nsim,size.of.sim))

        eta <- vector("list",nsim)
        theta <- vector("list",nsim)
        omega <- vector("list",nsim)
        sigma <- vector("list",nsim)
                   
        ##data <- read.table(filename,skip=1,header=TRUE)
        eta[[1]] <- read.table(paste(filename,".50",sep=""))
        ##ofv <- read.table(paste(filename,".51",sep=""))
        theta[[1]] <- read.table(paste(filename,".52",sep=""))
        ##setheta <- read.table(paste(filename,".53",sep=""))
        omega[[1]] <- read.table(paste(filename,".54",sep=""))
        ##seomega <- read.table(paste(filename,".55",sep=""))
        sigma[[1]] <- read.table(paste(filename,".56",sep=""))
        ##sesigma <- read.table(paste(filename,".57",sep="")) 

        tables.read <- TRUE
      }
    } else { # new file convention
      est.file <- paste(filename,est.tab.suffix,sep="")
      deriv.file <- paste(filename,deriv.tab.suffix,sep="")
      if ((is.cwres.readable.file(est.file)) &&
          (is.cwres.readable.file(deriv.file))) {

        nsim <- 1
        
        #########################
        ## read derivatives table
        #########################
        
        ## Check if we have unequal number of fields in the derivative file
        ## used for multiple simulations
        num.fields <- count.fields(deriv.file,skip=1)
        if((length(unique(num.fields))!=1) ||
           ## These lines fix problem if num_fields matches fields first header row 
           (unique(num.fields) == 8) ||
           (unique(num.fields) == 3)) { 
          
          tmp   <- readLines(deriv.file, n = -1)
          inds  <- grep("TABLE",tmp)
          if (length(inds)!=1){
            inds  <- inds[c(2:length(inds))]
            inds2 <- inds+1
            tempfile<- paste(deriv.file,".xptmp",sep="")
            write.table(tmp[-c(inds,inds2)],file=tempfile,
                        row.names=FALSE,quote=FALSE)
            data <- read.table(tempfile,skip=2,header=T)
            unlink(tempfile)

            ## add iteration label to data
            nsim <- length(inds)+1       
          
          } else {
            data <- read.table(deriv.file,skip=1,header=T)
          }
        } else {
          data <- read.table(deriv.file,skip=1,header=T)
        }
        size.of.sim <- dim(data)[1]/nsim
        data[,"iteration.number"] <- sort(rep(1:nsim,size.of.sim))

        
        #########################
        ## read estimated parameters table
        #########################        
        filename.extra <- est.file
        data.extra <- scan(filename.extra,
                           sep = "\n", what = character(),
                           quiet=TRUE)
        eta.pat <- "^ *ETAS"
        theta.pat <- "^ *THETAS"
        omega.pat <- "^ *OMEGAS"
        sigma.pat <- "^ *SIGMAS"

        eta.pat.line <- grep(eta.pat, data.extra)
        theta.pat.line <- grep(theta.pat, data.extra)
        omega.pat.line <- grep(omega.pat, data.extra)
        sigma.pat.line <- grep(sigma.pat, data.extra)

        pat.lines <- list(eta.pat.line,
                           theta.pat.line,
                           omega.pat.line,
                           sigma.pat.line)
        pattern.lengths <- sapply(pat.lines,length)
        if(!all(pattern.lengths==nsim)){
          cat(paste("The",est.file,"and",deriv.file,
              "files do not match in size\n"))
          return(NULL)
        }

        eta <- vector("list",nsim)
        theta <- vector("list",nsim)
        omega <- vector("list",nsim)
        sigma <- vector("list",nsim)
        for(i in 1:nsim){
          tot.eta.rows <- theta.pat.line[i] - eta.pat.line[i] - 1
          tot.theta.rows <- omega.pat.line[i] - theta.pat.line[i] - 1
          tot.omega.rows <- sigma.pat.line[i] - omega.pat.line[i] - 1
          if(i==nsim){
            tot.sigma.rows <- length(data.extra) - sigma.pat.line[i]
          } else {
            tot.sigma.rows <- eta.pat.line[i+1] - sigma.pat.line[i] - 1
          }
          eta[[i]] <- read.table(filename.extra,skip=eta.pat.line[i],nrows=tot.eta.rows)
          theta[[i]] <- read.table(filename.extra,skip=theta.pat.line[i],nrows=tot.theta.rows)
          omega[[i]] <- read.table(filename.extra,skip=omega.pat.line[i],nrows=tot.omega.rows)
          sigma[[i]] <- read.table(filename.extra,skip=sigma.pat.line[i],nrows=tot.sigma.rows)
        }
        tables.read <- TRUE
      }
    } # end new file convention
    
    if (tables.read){
      
      setClass("nm.data",
               representation(data = "data.frame"
                              ,eta  = "data.frame"
                              ##,ofv = "data.frame"
                              ,theta = "data.frame"
                              ##,setheta = "data.frame"
                              ,omega = "data.frame"
                              ##,seomega = "data.frame"
                              ,sigma = "data.frame"
                              ##,sesigma = "data.frame",
                              ),
               where = .GlobalEnv
               )

      all.data <- vector("list",nsim)
      for(i in 1:nsim){
        for (j in names(eta[[i]])){
          if (!is.numeric(eta[[i]][[j]])) {
            cat(paste("Some values of ETA in the NONMEM produced file",
                      filename.extra,
                      "are non-numeric\n"))
            cat("Attempting to fix the problem\n")
            ## add an 'E' where NONMEM has dropped one
            ## i.e. if the number is  -0.2404305E-105
            ## NONMEM writes  -0.2404305-105
            exp.pat <- "([0-9]+)([+-][0-9]+)"
            repl.exp.pat <- "\\1E\\2"
            #bad.locs <- grep(exp.pat,eta[[i]][[j]])
            new.vals <- gsub(exp.pat,repl.exp.pat,eta[[i]][[j]])
            eta[[i]][j] <- as.numeric(new.vals)
         }
        }
        
        all.data[[i]] <- new("nm.data"
                             ,data=data[data$iteration.number==i,]
                             ,eta=eta[[i]]
                             ##,ofv=ofv
                             ,theta=theta[[i]]
                             ##,setheta=setheta
                             ,omega=omega[[i]]
                             ##,seomega=seomega
                             ,sigma=sigma[[i]]
                             ##,sesigma=sesigma
                             )
      }
      
      if(is.null(data$MDV)){
        cat("Assuming all dataset lines contain a measurement\n")
        cat("No MDV item defined in dataset\n")
      }

      if(is.null(data$DV)){
        dv.found <- FALSE
        cat("No DV item defined in dataset\n")

        ## use this code to det the DV to the value at the position
        ## just beofre PRED RES and WRES
        ##
        ## doesn't work right now!
        ##
        #append.list <- c("PRED","RES","WRES")
        #data.names <- names(data)
        #tmp.length <- length(data.names)
        #if (all(data.names[(tmp.length-3):(tmp.length-1)]==append.list)){
        #  data[,"DV"]=data[[data.names[(tmp.length-4)]]]
        #  dv.found <- TRUE
        #  cat(paste("Using", data.names[(tmp.length-4)],"as the DV value\n"))
        #}
        
        if(!dv.found){
          all.data <- NULL
        }
      }
      
    } else { # no tables read      
      all.data <- NULL
      cat("Some or all of the required CWRES tables are\n", 
          "missing. Please see the online help for details\n", 
          "of what is required (?compute.cwres).\n")
    }
    return(all.data)
  }



ind.cwres <-
  function(ind.data,
           H.names,
           G.names,
           OMEGA,
           SIGMA,
           IND.ETAS,
           ...){
    
    if(is.null(ind.data$MDV)){
      ind.data1 <- ind.data
    } else {
      ind.data1 <- ind.data[ind.data$MDV==0,]
    }
    
    
    if (nrow(ind.data1)!=0) {
      
      ## for testing
      ##cat(paste("ID",ind.data1[1,"ID"],"\n"))
      
      ## create the epsilon linearization matrix H
      H.EPS = as.matrix(subset(ind.data1,select=H.names))
      
      ## create the eta linearization matrix G
      G.ETA = as.matrix(subset(ind.data1,select=G.names))
      
      ## covariance matrix
      ## this is only true for NON-L2 type data
      ## L2 data will have the first term non-diagonal
      TMP <- diag(H.EPS %*% SIGMA %*% t(H.EPS))
      IND.COV = diag(TMP,nrow=length(TMP)) +
        G.ETA %*% OMEGA %*% t(G.ETA)
      
      ## FOCE residuals
      #browser()
      if(any(grepl("^IPRED$",names(ind.data1)))) EXP.F <- as.matrix(ind.data1$IPRED) - G.ETA %*% IND.ETAS
      if(any(grepl("^IPRE$",names(ind.data1)))) EXP.F <- as.matrix(ind.data1$IPRE) - G.ETA %*% IND.ETAS
      #EXP.F <- as.matrix(ind.data1$IPRE) - G.ETA %*% IND.ETAS 
      FOCE.RES <- as.matrix(ind.data1$DV) - EXP.F
      
      
      ## CWRES
      SQRT.IND.COV <- sqrtm(IND.COV)
      IND.CWRES <-  solve(SQRT.IND.COV,FOCE.RES)
      
      ## add zero values back again
      if(is.null(ind.data$MDV)){
      } else {
        CWRES <- rep(0,length(ind.data[,1]))
        ind.data2 <- cbind(ind.data,CWRES)
        ind.data2[ind.data2$MDV==0,"CWRES"] <- IND.CWRES
        IND.CWRES <- as.matrix(ind.data2["CWRES"])
      }      
    }
    else {
      CWRES <- rep(0, length(ind.data[, 1]))
      ind.data2 <- cbind(ind.data, CWRES)
      IND.CWRES <- as.matrix(ind.data2["CWRES"])
    }
    return(IND.CWRES)
  }

      
compute.cwres <-
  function(run.number,
           tab.prefix="cwtab",
           sim.suffix="",
           est.tab.suffix=".est",
           deriv.tab.suffix=".deriv",
           old.file.convention=FALSE,
           id="ALL",
           printToOutfile=TRUE,
           onlyNonZero=TRUE,
           ...){


    out.file = paste(tab.prefix,run.number,sim.suffix,sep="")

    full.dataset <- read.cwres.data(out.file,
                                    old.file.convention=old.file.convention,
                                    est.tab.suffix=est.tab.suffix,
                                    deriv.tab.suffix=deriv.tab.suffix,
                                    ...)
    
    if (is.null(full.dataset)) {
      return()
    }

    num.reps <- length(full.dataset)
    tot.cwres <- c()
    for(rep in 1:num.reps){

      ## get one of the repetitions
      dataset <- full.dataset[[rep]]

      ##    dataset <- read.cwres.data("sdtab1") # only for testing
      first.only.data <- dataset@data[!duplicated(dataset@data$ID),]
      all.etas <- dataset@eta
      all.etas <- cbind(first.only.data["ID"], all.etas)
      

      ## create OMEGA, SIGMA  matricies
      OMEGA <- as.matrix(dataset@omega)
      SIGMA <- as.matrix(dataset@sigma)
      
      ## create the names of the H and G columns in the dataset
      H.names = c()
      i=1
      while(i<(length(dataset@sigma)+1)){
        H.names = c(H.names,paste("H",i,"1",sep=""))
        i=i+1
      }
      
      G.names = c()
      i=1
      while(i<(length(dataset@omega)+1)){
        G.names = c(G.names,paste("G",i,"1",sep=""))
        i=i+1
      }
      

      if(id=="ALL"){
        

        id.vals <- unique(dataset@data$ID)
        CWRES <- c()
        for(i in id.vals){
          #browser()
          #ind.data <- subset(dataset@data,ID==i)
          ind.data <- dataset@data[dataset@data$ID==i,]
          ind.etas <- t(as.matrix(all.etas[all.etas$ID==i,colnames(all.etas)!="ID"]))
          CWRESI <- ind.cwres(ind.data,
                              H.names,
                              G.names,
                              OMEGA,
                              SIGMA,
                              ind.etas,
                              ...)
          CWRES <- c(CWRES,CWRESI)
        }
        CWRES <- as.matrix(CWRES)
        
        
        
        if(printToOutfile==TRUE){
          
          ## set up out file name
          if(old.file.convention){
            filename <- paste(out.file,".cwres",sep="")
          } else {
            filename <- out.file
          }
          ## bind CWRES to data file
                                        #data.cwres <- cbind(dataset@data,CWRES)
          data.cwres <- data.frame("ID"=dataset@data$ID)
          if(!is.null(dataset@data$MDV)) data.cwres$MDV=dataset@data$MDV
          if(!is.null(dataset@data$DV)) data.cwres$DV=dataset@data$DV
          if(any(grepl("^IPRE$",names(dataset@data)))){
            if(!is.null(dataset@data$IPRE)) data.cwres$IPRE=dataset@data$IPRE
          }
          if(any(grepl("^IPRED$",names(dataset@data)))){
            if(!is.null(dataset@data$IPRED)) data.cwres$IPRED=dataset@data$IPRED
          }
          if(!is.null(dataset@data$WRES)) data.cwres$WRES=dataset@data$WRES
          if(!is.null(CWRES)) data.cwres$CWRES=CWRES
          
          ## data.cwres <- data.frame("ID"=dataset@data$ID,
          ##                                  "DV"=dataset@data$DV,
          ##                                  "MDV"=dataset@data$MDV,
          ##                                  "IPRE"=dataset@data$IPRE,
          ##                                  "WRES"=dataset@data$WRES,
          ##                                  "CWRES"=CWRES)

          ## print out dataset
          tmp <- installed.packages(priority="NA")
          if (length(grep("xpose4",tmp))>0){
            xpose.version <- tmp["xpose4","Version"]
            xpose.text <- paste("from Xpose version",xpose.version,sep=" ")
          } else {
            xpose.text <- ""
          }
          if (rep==1) {
            append.table.message <- FALSE
          } else {
            append.table.message <- TRUE
          }
          cat(paste("TABLE for CWRES computed using compute.cwres.R",xpose.text,"on",
                    format(Sys.time(), "%a %b %d, %X, %Y"),"\n"),
              file=filename,
              append=append.table.message)
          newdata <- format(data.cwres,sci=TRUE)
          suppressWarnings(
                           write.table(newdata,filename,row.names=FALSE,sep=" ",
                                       quote=FALSE,append=TRUE)
                           )
        } # end print to out file

        if(onlyNonZero==TRUE){
          if(is.null(dataset@data$MDV)){
          } else {
            ## bind the data to the data file
            data.cwres <- cbind(dataset@data,CWRES)
            tmp <- data.cwres[data.cwres$MDV==0,]
            #tmp <- subset(data.cwres,MDV==0)
            CWRES <- tmp$CWRES
          }
        }
        
      } else { ## end if ID==ALL
        data1 <- dataset@data[dataset@data$ID==id,]
        ind.etas <- t(as.matrix(all.etas[all.etas$ID==id,colnames(all.etas)!="ID"]))
        CWRES <- ind.cwres(data1,
                           H.names,
                           G.names,
                           OMEGA,
                           SIGMA,
                           ind.etas,
                           ...)

        if(onlyNonZero==TRUE){
          if(is.null(data1$MDV)){
          } else {
            ## bind the data to the data file
            data1.cwres <- cbind(data1,CWRES)
            tmp <- data1.cwres[data1.cwres$MDV==0,]
            #tmp <- subset(data1.cwres,MDV==0)
            CWRES <- tmp$CWRES
          }
        }
      } # end if ID=id
     
      tot.cwres <- c(tot.cwres,CWRES)
    }# end num.reps 
    
    return(tot.cwres)

  }



