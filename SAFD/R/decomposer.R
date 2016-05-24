decomposer <-
function(X){
 #function decomposes a given fuzzy set X in the basis-represantation 
 ok<-checking(X)
 if(ok==1){
  nl<-nrow(X)/2
   vc<-mean(X$x[nl:(nl+1)])
   a<-X$x[1:nl]
   yl<-vc-a
   cl<-yl-c(yl[2:nl],0)

   b<-X$x[(nl+1):(2*nl)]
   yr<-b-vc
   cr<-yr-c(0,yr[1:(nl-1)])
 
   Coordinates<-data.frame(coor=c(cl,vc,cr),alpha=c(X$alpha[1:nl],1,X$alpha[(nl+1):(2*nl)]))
  invisible(Coordinates)
 }
}
