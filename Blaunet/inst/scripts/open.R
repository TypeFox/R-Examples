loadfile <- function(h,...) {
  f <- gfile(text="Select a file", type="open")
  if (f=="") return;
  if (grepl('[:punct:.]rda',f) | grepl('[:punct:.]Rdata',f) | grepl('[:punct:.]dta',f) | grepl('[:punct:.]csv',f) | grepl('[:punct:.]sav',f) | grepl('[:punct:.]xpt',f) | grepl('[:punct:.]dat',f) | grepl('[:punct:.]DAT',f) | grepl('[:punct:.]txt',f)) {
    if (grepl('[:punct:.]rda',f) | grepl('[:punct:.]Rdata',f)) {
      c1 <- load(f)
      load(f,envir=.GlobalEnv)
      if (class(get(c1))=="data.frame" | class(get(c1))=="matrix") {
        assign("cov",get(c1),envir=.GlobalEnv)
      } else if (class(get(c1))=="list") {
        assign("adj",get(c1)$adj,envir=.GlobalEnv)
        assign("el",get(c1)$el,envir=.GlobalEnv)
        assign("cov",get(c1)$square.data,envir=.GlobalEnv)
      }
      gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
    } else if (grepl('[:punct:.]dta',f)) {
      c1 <- read.dta(f)
      assign("cov",c1,envir=.GlobalEnv)
      gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
    } else if (grepl('[:punct:.]csv',f)) {
      c1 <- read.csv(f,header=T,sep=",")
      assign("cov",c1,envir=.GlobalEnv)
      gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
    } else if (grepl('[:punct:.]sav',f)) {
      c1 <- read.spss(f)
      assign("cov",c1,envir=.GlobalEnv)
      gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
    } else if (grepl('[:punct:.]xpt',f)) {
      c1 <- read.xport(f)
      assign("cov",c1,envir=.GlobalEnv)
      gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
    } else if (grepl('[:punct:.]dat',f) | grepl('[:punct:.]DAT',f) | grepl('[:punct:.]txt',f)) {
      c1 <- read.table(f)
      assign("cov",c1,envir=.GlobalEnv)
      gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
    }
  } else gmessage("Sorry! Unknown file format.", parent = window)
}

loadnet <- function(h,...) {
  f <- gfile(text="Select a file", type="open")
  if (f=="") return;
  if (grepl('[:punct:.]rda',f) | grepl('[:punct:.]Rdata',f) | grepl('[:punct:.]paj',f) | grepl('[:punct:.]dta',f) | grepl('[:punct:.]csv',f) | grepl('[:punct:.]sav',f) | grepl('[:punct:.]xpt',f) | grepl('[:punct:.]dat',f) | grepl('[:punct:.]DAT',f) | grepl('[:punct:.]net',f) | grepl('[:punct:.]txt',f)) {
    if (grepl('[:punct:.]rda',f) | grepl('[:punct:.]Rdata',f)) {
      net <- network(get(load(f)))
    } else if (grepl('[:punct:.]paj',f)) {
      net <- network(read.paj(f))
    } else if (grepl('[:punct:.]dta',f)) {
      net <- network(read.dta(f))
    } else if (grepl('[:punct:.]csv',f)) {
      net <- network(read.csv(f,header=F,sep=","))
    } else if (grepl('[:punct:.]sav',f)) {
      net <- network(read.spss(f))
    } else if (grepl('[:punct:.]xpt',f)) {
      net <- network(read.xport(f))
    } else if (grepl('[:punct:.]dat',f) | grepl('[:punct:.]DAT',f) | grepl('[:punct:.]net',f) | grepl('[:punct:.]txt',f)) {
      net <- network(read.table(f))
    }
    el <- data.frame(as.matrix(net,matrix.type='edgelist'))
    names(el) <- c("i","j") 
    adj <- data.frame(as.matrix(net))
    rownames(adj) <- colnames(adj) <- attr(as.matrix(net,matrix.type='edgelist'),"vnames")
    assign("adj",adj,envir=.GlobalEnv)
    assign("el",el,envir=.GlobalEnv)
    gmessage(paste("Congratulations! Your attribute file ",f," is now loaded",sep=''), parent = window)
  } else gmessage("Sorry! Unknown file format.", parent = window)
}


