summary.lba.ls.logit <- function(object,
                                 digits=2L,
                                 ...){
 logitA <- getCall(object)$logitA
 logitB <- getCall(object)$logitB

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
 cat('|             MIXING PARAMETERS             |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[4]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|              LATENT BUDGETS               |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[5]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|           RESCALED LATENT BUDGETS         |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[6]],digits=digits))
 cat('\n')

 cat('|-------------------------------------------|\n')
 cat('|             BUDGET PROPORTIONS            |\n')
 cat('|-------------------------------------------|\n')
 print(round(object[[7]],digits=digits))
 cat('\n')

 if(!is.null(logitA) & !is.null(logitB)){

  cat('|-------------------------------------------|\n')
  cat('|            LOGIT PARAMETERS OMSK          |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[10]],digits=digits))
  cat('\n')

  cat('|-------------------------------------------|\n')
  cat('|            LOGIT PARAMETERS PSITK         |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[11]],digits=digits))
  cat('\n')

 } else

 #===========================================================================
 #Multinomial logit on mixing parameters only
 #===========================================================================

 if(!is.null(logitA) & is.null(logitB)){

  cat('|-------------------------------------------|\n')
  cat('|            LOGIT PARAMETERS OMSK          |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[10]],digits=digits))
  cat('\n')

 } else {

  #===========================================================================
  #Multinomial logit on latent components only
  #===========================================================================

  cat('|-------------------------------------------|\n')
  cat('|            LOGIT PARAMETERS PSITK         |\n')
  cat('|-------------------------------------------|\n')
  print(round(object[[10]],digits=digits))
  cat('\n')

 }

 func <- formatC(round(object[[8]],digits=digits),
                 format='f',
                 digits=digits)
 cat('|-------------------------------------------|\n')
 cat(paste('| VALUE OF THE LS FUNCTION:',sep=' '),
     func,
     paste(rep(' ',(45-(31+nchar(func)))),collapse=''), '|\n')
 cat('| NUMBER OF ITERACTIONS:',round(object[[9]],digits=digits),
     paste(rep(' ',(45-(28+nchar(object[[9]])))),collapse=''),'|\n')
 cat('|-------------------------------------------|\n') 

}
