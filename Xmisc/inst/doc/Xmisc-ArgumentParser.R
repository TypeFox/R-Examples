### R code from vignette source 'Xmisc-ArgumentParser.Rnw'

###################################################
### code chunk number 1: Xmisc-ArgumentParser.Rnw:35-36
###################################################
  cat(as.character(packageVersion('Xmisc')))


###################################################
### code chunk number 2: Xmisc-ArgumentParser.Rnw:40-41
###################################################
  cat(unlist(strsplit(packageDescription('Xmisc')[['Date']],' '))[1])


###################################################
### code chunk number 3: Xmisc-ArgumentParser.Rnw:150-152 (eval = FALSE)
###################################################
##   ## install Xmisc
##   install.packages("Xmisc")


###################################################
### code chunk number 4: Xmisc-ArgumentParser.Rnw:163-165
###################################################
  require(Xmisc)
  parser <- ArgumentParser$new()  


###################################################
### code chunk number 5: Xmisc-ArgumentParser.Rnw:171-196
###################################################
## add a character object
  parser$add_argument(
    '--a_name',type='character',
    help='A a_name.'
  )  
## add an integer object with default
  parser$add_argument(
    '--a_int',type='integer',
    default=1,    
    help='A integer.'
  )

## add a numeric object with default
  parser$add_argument(
    '--a_num',type='numeric',
    default=1,    
    help='A number.'
  )

## add a logical object with default
  parser$add_argument(
    '--if.test',type='logical',
    default=FALSE,
    help='Whether it is a test?!'
  )


###################################################
### code chunk number 6: Xmisc-ArgumentParser.Rnw:202-203
###################################################
  parser$add_usage('Xmisc-ArgumentParser.R [options]')


###################################################
### code chunk number 7: Xmisc-ArgumentParser.Rnw:210-212
###################################################
  parser$add_description(
    'An executable R script parsing arguments from Unix-like command line.')


###################################################
### code chunk number 8: Xmisc-ArgumentParser.Rnw:219-230
###################################################
  parser$add_argument(
    '--h',type='logical',
    action='store_true',
    help='Print the help page'
  )

  parser$add_argument(
    '--help',type='logical',
    action='store_true',
    help='Print the help page'
  )


###################################################
### code chunk number 9: Xmisc-ArgumentParser.Rnw:239-240
###################################################
  parser$get_args()


###################################################
### code chunk number 10: Xmisc-ArgumentParser.Rnw:249-250
###################################################
  parser$helpme()


###################################################
### code chunk number 11: Xmisc-ArgumentParser.Rnw:262-266
###################################################
system.file('bin', 'Xmisc-ArgumentParser.R', package='Xmisc', mustWork=TRUE)

## Or,
Xmisc::get_executable('Xmisc','Xmisc-ArgumentParser.R')


###################################################
### code chunk number 12: Xmisc-ArgumentParser.Rnw:427-433 (eval = FALSE)
###################################################
## ## add a required character object
##   parser$add_argument(
##     '--a_name',type='character',
##     required=TRUE, 
##     help='A a_name.'
##   )  


