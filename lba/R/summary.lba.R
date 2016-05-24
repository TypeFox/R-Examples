summary.lba <- function(object,
                        digits=2L,
                        ...){

 cat('|-------------------------------------------|\n')
 cat('|         COMPOSITION DATA MATRIX           |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[1]],digits=digits))
 cat('\n') 

 cat('|-------------------------------------------|\n')
 cat('|              EXPECTED BUDGET              |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[2]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|              RESIDUAL MATRIX              |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[3]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|      UNIDENTIFIED MIXING PARAMETERS       |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[4]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|       UNIDENTIFIED LATENT BUDGETS         |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[5]],digits=digits))
 cat('\n')

 if(object$what == 'outer'){
  cat('|-------------------------------------------|\n')
  cat('|     OUTER EXTREME MIXING PARAMETERS       |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[6]],digits=digits))
  cat('\n')

  cat('|-------------------------------------------|\n')
  cat('|       OUTER EXTREME LATENT BUDGETS        |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[7]],digits=digits))
  cat('\n')
 } else {
  cat('|-------------------------------------------|\n')
  cat('|     INNER EXTREME MIXING PARAMETERS       |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[6]],digits=digits))
  cat('\n')

  cat('|-------------------------------------------|\n')
  cat('|       INNER EXTREME LATENT BUDGETS        |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[7]],digits=digits))
  cat('\n')
 }

 cat('|-------------------------------------------|\n')
 cat('|         RESCALED LATENT BUDGET            |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[8]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|             BUDGET PROPORTIONS            |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[9]],digits=digits))
 cat('\n')

 ifelse(any(class(object)=='lba.ls'),
        vof <- 'LS',
        vof <- '-LOGLF')

 func <- formatC(round(object[[10]],digits=digits),
                 format='f',
                 digits=digits)
 cat('|-------------------------------------------|\n')
 cat(paste('| VALUE OF THE',vof,'FUNCTION:',sep=' '),
     func,
     paste(rep(' ',(45-(29+nchar(vof)+nchar(func)))),collapse=''), '|\n')
 cat('| NUMBER OF UNIDENTIFIED ITERACTIONS:',round(object[[11]],digits=digits),
     paste(rep(' ',(45-(41+nchar(object[[11]])))),collapse=''),'|\n')

 if(ncol(object[[6]])> 2){

   cat('| NUMBER OF IDENTIFIED ITERACTIONS:',round(object[[12]],digits=digits),
       paste(rep(' ',(45-(38+nchar(object[[12]])))),collapse=''),'|\n') 
  
 }
 cat('|-------------------------------------------|\n') 

}
