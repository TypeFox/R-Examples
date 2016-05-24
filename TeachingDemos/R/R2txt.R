### consider adding option to include errors Can implement by using
### options(error=newfunction) and newfunction would use the
### savehistory command to get the expression and geterrmessage to get
### the error message.  Warnings can be included by checking to see if
### last.warning has changed, use print.warnings to format.

R2txt.vars <- new.env()

R2txt <- function(cmd,res,s,vis) {
  if(R2txt.vars$first) {
      R2txt.vars$first <- FALSE
      if( R2txt.vars$res ) {
          sink()
          close(R2txt.vars$outcon)
          R2txt.vars$outcon <- textConnection(NULL, open='w')
          sink(R2txt.vars$outcon, split=TRUE)
      }
  } else {

      if( R2txt.vars$cmd ){
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', paste("\n",R2txt.vars$continue, sep=''),
                          cmdline)
          cmdline <- gsub('}', paste("\n",R2txt.vars$continue,"}", sep=''),
                          cmdline)
          cat(R2txt.vars$prompt, cmdline, "\n", sep='',
              file=R2txt.vars$con)
      }
      if( R2txt.vars$cmdfile ) {
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', "\n ", cmdline)
          cmdline <- gsub('}', "\n}", cmdline)
          cat(cmdline,"\n", file=R2txt.vars$con2)
      }

      if( R2txt.vars$res ) {
          tmp <- textConnectionValue(R2txt.vars$outcon)
          if(length(tmp)) {
              cat(tmp,sep='\n',file=R2txt.vars$con)
              sink()
              close(R2txt.vars$outcon)
              R2txt.vars$outcon <- textConnection(NULL, open='w')
              sink(R2txt.vars$outcon, split=TRUE)
          }
      }
  }

  TRUE
}

txtStart <- function(file, commands=TRUE, results=TRUE, append=FALSE,
                     cmdfile, visible.only=TRUE) {

  tmp <- TRUE
  if(is.character(file)){
    if(append){
      con <- file(file,open='a')
    } else {
      con <- file(file,open='w')
    }
    tmp <- FALSE
  } else if( any( class(file) == 'connection' ) ) {
    con <- file
  } else {
    stop('file must be a character string or connection')
  }
  if( tmp && isOpen(con) ) {
    R2txt.vars$closecon <- FALSE
  } else {
    R2txt.vars$closecon <- TRUE
    if(tmp){
        if(append) {
            open(con, open='a')
        } else {
            open(con, open='w')
        }
    }
  }
  R2txt.vars$vis <- visible.only
  R2txt.vars$cmd <- commands
  R2txt.vars$res <- results
  R2txt.vars$con <- con
  R2txt.vars$first <- TRUE

  if(results) {
      R2txt.vars$outcon <- textConnection(NULL, open='w')
      sink(R2txt.vars$outcon, split=TRUE)
  }

  if( !missing(cmdfile) ) {
    tmp <- TRUE
    if(is.character(cmdfile)) {
      con2 <- file(cmdfile, open='w')
      tmp <- FALSE
    } else if( any( class(cmdfile) == 'connection' ) ) {
      con2 <- cmdfile
    }
    if( tmp && isOpen(con2) ) {
      R2txt.vars$closecon2 <- FALSE
    } else {
      R2txt.vars$closecon2 <- TRUE
      if(tmp) {
          open(con2, open='w')
      }
    }
    R2txt.vars$con2 <- con2
    R2txt.vars$cmdfile <- TRUE
  } else {
    R2txt.vars$cmdfile <- FALSE
  }

  R2txt.vars$prompt <- unlist(options('prompt'))
  R2txt.vars$continue <- unlist(options('continue'))

  options(prompt= paste('txt',R2txt.vars$prompt,sep=''),
          continue= paste('txt',R2txt.vars$continue,sep='') )

  cat('Output being copied to text file,\nuse txtStop to end\n')
  addTaskCallback(R2txt, name='r2txt')
  invisible(NULL)
}

txtStop <- function() {
  removeTaskCallback('r2txt')
  if( R2txt.vars$closecon ) {
    close( R2txt.vars$con )
  }
  if( R2txt.vars$cmdfile && R2txt.vars$closecon2 ) {
    close( R2txt.vars$con2 )
  }
  options( prompt=R2txt.vars$prompt,
           continue=R2txt.vars$continue )
  if(R2txt.vars$res) {
      sink()
      close(R2txt.vars$outcon)
  }
  evalq( rm(list=ls()), envir=R2txt.vars )
  invisible(NULL)
}

txtComment <- function(txt,cmdtxt) {
    R2txt.vars$first <- TRUE
    if(!missing(txt)) {
        cat("\n",txt,"\n\n", file=R2txt.vars$con)
    }
    if(!missing(cmdtxt)) {
        cat("# ",cmdtxt,"\n", file=R2txt.vars$con2)
    }
}

txtSkip <- function(expr) {
    R2txt.vars$first <- TRUE
    expr
}


######### etxt extended or enscriptable

R2etxt <- function(cmd,res,s,vis) {
  if(R2txt.vars$first) {
      R2txt.vars$first <- FALSE
      if( R2txt.vars$res ) {
          sink()
          close(R2txt.vars$outcon)
          R2txt.vars$outcon <- textConnection(NULL, open='w')
          sink(R2txt.vars$outcon, split=TRUE)
      }
  } else {

      if( R2txt.vars$cmd ) {
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', "\n", cmdline)
          cmdline <- gsub('}', "\n}", cmdline)
          writeChar("",R2txt.vars$con)
          cat(R2txt.vars$cmdbg,file=R2txt.vars$con)
          writeChar("",R2txt.vars$con)
          cat(R2txt.vars$cmdcol,file=R2txt.vars$con)
          cat(R2txt.vars$prompt, cmdline, "\n", sep='',
              file=R2txt.vars$con)
      }
      if( R2txt.vars$cmdfile ) {
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', "\n ", cmdline)
          cmdline <- gsub('}', "\n}", cmdline)
          cat(cmdline,"\n", file=R2txt.vars$con2)
      }
      if( R2txt.vars$res ) {
          tmp <- textConnectionValue(R2txt.vars$outcon)
          if(length(tmp)) {
              writeChar("",R2txt.vars$con)
              cat(R2txt.vars$resbg, file=R2txt.vars$con)
              writeChar("",R2txt.vars$con)
              cat(R2txt.vars$rescol,file=R2txt.vars$con)
              cat(tmp,sep='\n',file=R2txt.vars$con)
              sink()
              close(R2txt.vars$outcon)
              R2txt.vars$outcon <- textConnection(NULL, open='w')
              sink(R2txt.vars$outcon, split=TRUE)
          }
      }
  }

  TRUE
}

etxtStart <- function(dir=tempfile('etxt'), file='transcript.txt',
                      commands=TRUE, results=TRUE, append=FALSE,
                      cmdbg='white',cmdcol='red', resbg='white',
                      rescol='navy',combg='cyan',comcol='black',
                      cmdfile, visible.only=TRUE) {

  if( !file_test("-d", dir) ) {
      dir.create(dir)
  }
  tmp <- TRUE
  if(is.character(file)){
    file2 <- file.path(dir,file)
    if(append){
      con <- file(file2,open='a')
    } else {
      con <- file(file2,open='w')
    }
    tmp <- FALSE
    R2txt.vars$file2 <- file2
  } else if( any( class(file) == 'connection' ) ) {
    con <- file
  } else {
    stop('file must be a character string or connection')
  }
  if(tmp && isOpen(con)) {
    R2txt.vars$closecon <- FALSE
  } else {
    R2txt.vars$closecon <- TRUE
    if(tmp) {
        if(append) {
            open(con, open='a')
        } else {
            open(con, open='w')
        }
    }
  }
  R2txt.vars$dir <- dir
  R2txt.vars$vis <- visible.only
  R2txt.vars$cmd <- commands
  R2txt.vars$res <- results
  R2txt.vars$con <- con
  R2txt.vars$first <- TRUE

  tmp <- round( col2rgb( c(cmdbg,cmdcol,resbg,rescol,combg,comcol) )/255,
               3)
  tmp2 <- paste( rep(c('bgcolor{','color{'),3),
                 tmp[1,], ' ', tmp[2,], ' ', tmp[3,], '}',
                 sep='' )

  R2txt.vars$cmdbg  <- tmp2[1]
  R2txt.vars$cmdcol <- tmp2[2]
  R2txt.vars$resbg  <- tmp2[3]
  R2txt.vars$rescol <- tmp2[4]
  R2txt.vars$combg  <- tmp2[5]
  R2txt.vars$comcol <- tmp2[6]

  if(results) {
      R2txt.vars$outcon <- textConnection(NULL, open='w')
      sink(R2txt.vars$outcon, split=TRUE)
  }
  tmp3 <- TRUE
  if( !missing(cmdfile) ) {
    if(is.character(cmdfile)) {
      con2 <- file(cmdfile, open='w')
      tmp3 <- FALSE
    } else if( any( class(cmdfile) == 'connection' ) ) {
      con2 <- cmdfile
    }
    if( tmp3 && isOpen(con2) ) {
      R2txt.vars$closecon2 <- FALSE
    } else {
      R2txt.vars$closecon2 <- TRUE
      if(tmp3) {
          open(con2, open='w')
      }
    }
    R2txt.vars$con2 <- con2
    R2txt.vars$cmdfile <- TRUE
  } else {
    R2txt.vars$cmdfile <- FALSE
  }

  R2txt.vars$prompt <- unlist(options('prompt'))
  R2txt.vars$continue <- unlist(options('continue'))

  options(prompt= paste('etxt',R2txt.vars$prompt,sep=''),
          continue= paste('etxt',R2txt.vars$continue,sep='') )

  cat('Output being copied to text file,\nuse etxtStop to end\n')
  addTaskCallback(R2etxt, name='r2etxt')
  invisible(NULL)
}

etxtStop <- function() {
  removeTaskCallback('r2etxt')
  if( R2txt.vars$closecon ) {
    close( R2txt.vars$con )
  }
  if( R2txt.vars$cmdfile && R2txt.vars$closecon2 ) {
    close( R2txt.vars$con2 )
  }
  options( prompt=R2txt.vars$prompt,
           continue=R2txt.vars$continue )
  if(R2txt.vars$res) {
      sink()
      close(R2txt.vars$outcon)
  }

  if( 'file2' %in% names(R2txt.vars) ) {
      out <- R2txt.vars$file2
  } else {
      out <- invisible(NULL)
  }
  evalq( rm(list=ls()), envir=R2txt.vars )
  out
}

etxtComment <- function(txt,cmdtxt) {
    R2txt.vars$first <- TRUE
    if(!missing(txt)) {
        writeChar("",R2txt.vars$con)
        cat(R2txt.vars$combg,file=R2txt.vars$con)
        writeChar("",R2txt.vars$con)
        cat(R2txt.vars$comcol,file=R2txt.vars$con)
        cat("\n",txt,"\n\n", file=R2txt.vars$con)
    }
    if(!missing(cmdtxt)) {
        cat("# ",cmdtxt,"\n", file=R2txt.vars$con2)
    }
}

etxtSkip <- function(expr) {
    R2txt.vars$first <- TRUE
    expr
}

etxtPlot <- function(file=paste(tempfile('plot',R2txt.vars$dir),'.eps',sep=''),
                     width=4, height=4) {
    dev.copy2eps(file=file, height=height, width=width)
    writeChar("",R2txt.vars$con)
    cat('epsf{',file,'}\n', sep='', file=R2txt.vars$con)
    R2txt.vars$first <- TRUE
    invisible(NULL)
}

#### version for sending output to MSword

R2wdtxt <- function(cmd,res,s,vis) {

  requireNamespace('R2wd', quietly = TRUE)
  if(R2txt.vars$first) {
      R2txt.vars$first <- FALSE
      if( R2txt.vars$res ) {
          sink()
          close(R2txt.vars$outcon)
          R2txt.vars$outcon <- textConnection(NULL, open='w')
          sink(R2txt.vars$outcon, split=TRUE)
      }
  } else {

      if( R2txt.vars$cmd ){
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', paste("\n",R2txt.vars$continue, sep=''),
                          cmdline)
          cmdline <- gsub('}', paste("\n",R2txt.vars$continue,"}", sep=''),
                          cmdline)
          R2wd::wdVerbatim( paste(R2txt.vars$prompt, cmdline, sep=''),
                       fontsize=R2txt.vars$fontsize )
      }
      if( R2txt.vars$cmdfile ) {
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', "\n ", cmdline)
          cmdline <- gsub('}', "\n}", cmdline)
          cat(cmdline,"\n", file=R2txt.vars$con2)
      }

      if( R2txt.vars$res ) {
          tmp <- textConnectionValue(R2txt.vars$outcon)
          if(length(tmp)) {
              R2wd::wdVerbatim(paste(tmp,sep='\n'), fontsize=R2txt.vars$fontsize)
              sink()
              close(R2txt.vars$outcon)
              R2txt.vars$outcon <- textConnection(NULL, open='w')
              sink(R2txt.vars$outcon, split=TRUE)
          }
      }
  }

  TRUE
}

wdtxtStart <- function(commands=TRUE, results=TRUE, fontsize=9,
                     cmdfile, visible.only=TRUE) {

  if( !requireNamespace('R2wd', quietly = TRUE) ) stop('the R2wd package is required')

  R2wd::wdGet()

  R2txt.vars$vis <- visible.only
  R2txt.vars$cmd <- commands
  R2txt.vars$res <- results
  R2txt.vars$first <- TRUE
  R2txt.vars$fontsize <- fontsize

  if(results) {
      R2txt.vars$outcon <- textConnection(NULL, open='w')
      sink(R2txt.vars$outcon, split=TRUE)
  }

  if( !missing(cmdfile) ) {
    tmp <- TRUE
    if(is.character(cmdfile)) {
      con2 <- file(cmdfile, open='w')
      tmp <- FALSE
    } else if( any( class(cmdfile) == 'connection' ) ) {
      con2 <- cmdfile
    }
    if( tmp && isOpen(con2) ) {
      R2txt.vars$closecon2 <- FALSE
    } else {
      R2txt.vars$closecon2 <- TRUE
      if(tmp) {
          open(con2, open='w')
      }
    }
    R2txt.vars$con2 <- con2
    R2txt.vars$cmdfile <- TRUE
  } else {
    R2txt.vars$cmdfile <- FALSE
  }

  R2txt.vars$prompt <- unlist(options('prompt'))
  R2txt.vars$continue <- unlist(options('continue'))

  options(prompt= paste('wdTxt',R2txt.vars$prompt,sep=''),
          continue= paste('wdTxt',R2txt.vars$continue,sep='') )

  cat('Output being copied to text file,\nuse wdtxtStop to end\n')
  addTaskCallback(R2wdtxt, name='r2wdtxt')
  invisible(NULL)
}

wdtxtStop <- function() {
  removeTaskCallback('r2wdtxt')

  if( R2txt.vars$cmdfile && R2txt.vars$closecon2 ) {
    close( R2txt.vars$con2 )
  }
  options( prompt=R2txt.vars$prompt,
           continue=R2txt.vars$continue )
  if(R2txt.vars$res) {
      sink()
      close(R2txt.vars$outcon)
  }
  evalq( rm(list=ls()), envir=R2txt.vars )
  invisible(NULL)
}

wdtxtComment <- function(txt,cmdtxt) {
    requireNamespace('R2wd', quietly = TRUE)
    R2txt.vars$first <- TRUE
    if(!missing(txt)) {
        R2wd::wdParagraph()
        R2wd::wdBody(txt)
        R2wd::wdParagraph()
    }
    if(!missing(cmdtxt)) {
        cat("# ",cmdtxt,"\n", file=R2txt.vars$con2)
    }
}

wdtxtSkip <- function(expr) {

    R2txt.vars$first <- TRUE
    expr
}

wdtxtPlot <- function(height=5, width=5, pointsize=10) {
    requireNamespace('R2wd', quietly = TRUE)
    R2txt.vars$first <- TRUE

    tmp <- recordPlot()
    R2wd::wdPlot(tmp, plotfun=replayPlot, height=height, width=width,
           pointsize=pointsize)
}


########## mdtxt Use markdown

R2mdtxt <- function(cmd,res,s,vis) {
      if(R2txt.vars$first) {
      R2txt.vars$first <- FALSE
      if( R2txt.vars$res ) {
          sink()
          close(R2txt.vars$outcon)
          R2txt.vars$outcon <- textConnection(NULL, open='w')
          sink(R2txt.vars$outcon, split=TRUE)
      }
  } else {

      if( R2txt.vars$cmd ){
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', "\n", cmdline)
          cmdline <- gsub('}', "\n}", cmdline)
          cat("```r\n",file=R2txt.vars$con)
          cat(R2txt.vars$prompt, cmdline, "\n", sep='',
              file=R2txt.vars$con)
          cat("```\n\n", file=R2txt.vars$con)
      }
      if( R2txt.vars$cmdfile ) {
          cmdline <- deparse(cmd)
          cmdline <- gsub('    ', "\n ", cmdline)
          cmdline <- gsub('}', "\n}", cmdline)
          cat(cmdline,"\n", file=R2txt.vars$con2)
      }
      if( R2txt.vars$res ) {
          tmp <- textConnectionValue(R2txt.vars$outcon)
          if(length(tmp)) {
              cmdline <- deparse(cmd)
              if( grepl("^\\s*pand(er|oc)", cmdline)[1] ) {
                  cat("\n",file=R2txt.vars$con)
                  cat(tmp,sep='\n',file=R2txt.vars$con)
                  cat("\n\n",file=R2txt.vars$con)
                  sink()
                  close(R2txt.vars$outcon)
                  R2txt.vars$outcon <- textConnection(NULL, open='w')
                  sink(R2txt.vars$outcon, split=TRUE)
              } else {
                  cat("```\n",file=R2txt.vars$con)
                  cat(tmp,sep='\n',file=R2txt.vars$con)
                  cat("```\n\n",file=R2txt.vars$con)
                  sink()
                  close(R2txt.vars$outcon)
                  R2txt.vars$outcon <- textConnection(NULL, open='w')
                  sink(R2txt.vars$outcon, split=TRUE)
              }
          }
      }
  }

  TRUE
}

mdtxtStart <- function(dir=tempfile('mdtxt'), file='transcript.md',
                      commands=TRUE, results=TRUE, append=FALSE,
                      cmdfile, visible.only=TRUE) {

  if( !file_test("-d", dir) ) {
      dir.create(dir)
  }
  tmp <- TRUE
  if(is.character(file)){
    file2 <- file.path(dir,file)
    if(append){
      con <- file(file2,open='a')
    } else {
      con <- file(file2,open='w')
    }
    tmp <- FALSE
    R2txt.vars$file2 <- file2
  } else if( any( class(file) == 'connection' ) ) {
    con <- file
  } else {
    stop('file must be a character string or connection')
  }
  if(tmp && isOpen(con)) {
    R2txt.vars$closecon <- FALSE
  } else {
    R2txt.vars$closecon <- TRUE
    if(tmp) {
        if(append) {
            open(con, open='a')
        } else {
            open(con, open='w')
        }
    }
  }
  R2txt.vars$dir <- dir
  R2txt.vars$vis <- visible.only
  R2txt.vars$cmd <- commands
  R2txt.vars$res <- results
  R2txt.vars$con <- con
  R2txt.vars$first <- TRUE

  if(results) {
      R2txt.vars$outcon <- textConnection(NULL, open='w')
      sink(R2txt.vars$outcon, split=TRUE)
  }
  tmp3 <- TRUE
  if( !missing(cmdfile) ) {
    if(is.character(cmdfile)) {
      con2 <- file(cmdfile, open='w')
      tmp3 <- FALSE
    } else if( any( class(cmdfile) == 'connection' ) ) {
      con2 <- cmdfile
    }
    if( tmp3 && isOpen(con2) ) {
      R2txt.vars$closecon2 <- FALSE
    } else {
      R2txt.vars$closecon2 <- TRUE
      if(tmp3) {
          open(con2, open='w')
      }
    }
    R2txt.vars$con2 <- con2
    R2txt.vars$cmdfile <- TRUE
  } else {
    R2txt.vars$cmdfile <- FALSE
  }

  R2txt.vars$prompt <- unlist(options('prompt'))
  R2txt.vars$continue <- unlist(options('continue'))

  options(prompt= paste('mdtxt',R2txt.vars$prompt,sep=''),
          continue= paste('mdtxt',R2txt.vars$continue,sep='') )

  cat('Output being copied to text file,\nuse mdtxtStop to end\n')
  addTaskCallback(R2mdtxt, name='r2mdtxt')
  invisible(NULL)
}

mdtxtStop <- function() {
  removeTaskCallback('r2mdtxt')
  if( R2txt.vars$closecon ) {
    close( R2txt.vars$con )
  }
  if( R2txt.vars$cmdfile && R2txt.vars$closecon2 ) {
    close( R2txt.vars$con2 )
  }
  options( prompt=R2txt.vars$prompt,
           continue=R2txt.vars$continue )
  if(R2txt.vars$res) {
      sink()
      close(R2txt.vars$outcon)
  }

  if( 'file2' %in% names(R2txt.vars) ) {
      out <- R2txt.vars$file2
  } else {
      out <- invisible(NULL)
  }
  evalq( rm(list=ls()), envir=R2txt.vars )
  out
}

mdtxtComment <- function(txt,cmdtxt) {
    R2txt.vars$first <- TRUE
    if(!missing(txt)) {
        cat("\n",txt,"\n\n", file=R2txt.vars$con)
    }
    if(!missing(cmdtxt)) {
        cat("# ",cmdtxt,"\n", file=R2txt.vars$con2)
    }
}

mdtxtSkip <- function(expr) {
    R2txt.vars$first <- TRUE
    expr
}

mdtxtPlot <- function(file=tempfile('plot',R2txt.vars$dir,'.png'),
                     width=4, height=4) {
    file <- gsub("\\\\","/",file)
    dev.copy(png, file=file, height=height, width=width, units='in',
             res=300)
    dev.off()
    cat('![plot ',file,'](',file,') \\\n\n', sep='', file=R2txt.vars$con)
    R2txt.vars$first <- TRUE
    invisible(NULL)
}


