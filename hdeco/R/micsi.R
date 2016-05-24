"micsi" <-
function () {

  cat("\nDefine the Decomposition Path\n")
  CIM <- readline("Title: ")
  KI <- list()
  EXEK<-readline("One image => 1   ,  Several images => 2  :")
  if(EXEK=="1") EXEK <- 0 else EXEK<-1
  for (i in 1:999){
    VN <- readline("\nName of the variable (X/Y/Z/press <ENTER> to finish): ")
    if(nchar(VN)==0){
      attr(KI,"cim") <- CIM
      return(KI)
    }
    if(VN=="x") VN<-"X"
    if(VN=="y") VN<-"Y"
    if(VN=="z") VN<-"Z"
    cat("\nWhich one(s)? (number(s)):\n")
    MIK <- scan(file="",what=0,quiet=T)
    KUKK <- list(VARIAB=VN,WHICH=MIK,EXE=EXEK)
    KI <- c(KI,KUKK)
  }
}

