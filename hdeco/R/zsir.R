"zsir" <-
function (KIK=c("Z"==substring(names(.QND),1,1)),fnev="huhu",EGYKE=1) {

  # KIK IS A T/F VECTOR WHERE Z VARIABLES ARE TRUE
  # fnev IS THE OUTPUT BASE FILENAME
  # EGYKE IS THE NUMBER OF INPUT IMAGES

  # BUILD .ZS OBJECT BY CALLING zs()
  assign(".ZS",NULL,pos=1)

  KIK <- c("Z"==substring(names(.QND),1,1))
  zs(KIK,EGYKE)
  
  H <- length(.ZS)

  # PROVIDE SEPARATOR LINE IN OUTPUT FILE
  write("\n--------------------------------------------------------",ncol=33,file=fnev,append=T)

  if(H < 1) {
    cat("\nNo Go")
    write("\nThere is no valid Mutual Info decomposition. ",ncol=33,file=fnev,append=T)
    return()
  }

  cimki <- T
  if(EGYKE == 1) {
    cat("\nSingle image colour decomposition")
    for (h in 1:H) { 
      if(!is.na(sum(.ZS[[h]]))){
        if(cimki) {
          write("\nMutual Info can be decomposed as sum over all Color Categories:\nP(ColorCategory)* H|ColorCategory",ncol=33,file=fnev,append=T)
        }
        cimki <- F
        write(c("\nStep",h),ncol=33,file=fnev,append=T)
        kii(.ZS[[h]][,2],.ZS[[h]][,1],fnev)
      }
    }
  }
  else {
    cat("\nDouble image colour decomposition")
    for (h in 1:H) {
      if(!is.na(sum(.ZS[[h]]))){
        if(cimki) {
          write("\nMutual Info can be decomposed as sum over all Color Categories:\nP(ColorCategory)* H|ColorCategory",ncol=33,file=fnev,append=T)
        }
        cimki <- F
        write(c("\nStep",h-1),ncol=33,file=fnev,append=T)
        kii(.ZS[[h]][,2],.ZS[[h]][,1],fnev)
      }
    }
  }

  if(cimki) {
    write("\nThere is no valid Mutual Info decomposition. ",ncol=33,file=fnev,append=T)
  }

}

