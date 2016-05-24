matricear <- function(string, maxdist = 1){
  V2 <- c()
  if(maxdist==1){
    V1<-deletion1(string)
    V1<-as.data.table(V1)
    V1<-V1[,V2 := string]
  }
  else if(maxdist==2){
    V1<-deletion2(string)
    V1<-as.data.table(V1)
    V1<-V1[,V2 := string]
  }
  else if(maxdist==3){
    V1<-deletion3(string)
    V1<-as.data.table(V1)
    V1<-V1[,V2 := string]
  }
  else{stop("Argument 'maxdist' out of range (1:3)")}
  V1
}
