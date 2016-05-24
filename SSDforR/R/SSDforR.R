SSDforR <-
function(){
  #print(ls("package:SSDforR"))
  f<-ls("package:SSDforR")
  f <- f[ f != "undo" ] 
  f <- f[ f != "insert" ] 
  print(f)
}
