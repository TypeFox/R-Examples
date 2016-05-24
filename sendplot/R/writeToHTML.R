
# header v1
# no hyperlinks 
writeToHTML1 <- function(obj, DFs, iType){

  cdat = DFs$cdat
  ndat = DFs$ndat

  if(iType == "circle") writeCircle.1(DFs, cdat, ndat, obj) 
  if(iType == "rect") writeRect.1(DFs, cdat, ndat)
  if(iType == "poly") writePoly.1(DFs, cdat, ndat, obj)
  
}



# header v2
# hyperlink activity
writeToHTML2 <- function(obj, DFs, iType){
  
  
  cdat = DFs$cdat
  ndat = DFs$ndat
  
  if(iType == "circle") writeCircle.2(DFs, cdat, ndat, obj) 
  if(iType == "rect") writeRect.2(DFs, cdat, ndat, obj)
  if(iType == "poly") writePoly.2(DFs, cdat, ndat, obj)


}




