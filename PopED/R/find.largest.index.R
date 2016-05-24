find.largest.index <- function (func.str="sfg",lab="bpop",mat=F,mat.row=T) {
  if(is.function(func.str)){
    txt <- capture.output(func.str)
  } else {
    txt <- capture.output(eval(parse(text=func.str)))
  }
  txt <- grep(paste("^[^\\#]*",lab,"\\[",sep=""),txt,value=T)
  ind <- 0
  if(length(txt)!=0 && !mat)  ind <- gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
  if(length(txt)!=0 && mat && mat.row)  ind <- gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*,.*?\\].*",sep=""),"\\1",txt)
  if(length(txt)!=0 && mat && !mat.row)  ind <- gsub(paste("^[^\\#]*",lab,"\\[.*?,\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
  
  max(as.numeric(ind))
}
# 
#  find.largest.index("sfg","bpop")
#  find.largest.index("sfg","b")
#find.largest.index("sfg","bocc",mat=T,mat.row=T)
#  find.largest.index("sfg","x")
#  find.largest.index("sfg","a")
# 
# txt <- capture.output(eval(parse(text="sfg")))
# txt <- grep(paste("^[^\\#]*bpop","\\[",sep=""),txt,value=T)
# txt
# ind <- gsub(paste("^[^\\#]*","bpop","\\[","(\\d+)\\].*",sep=""),"\\1",txt)
# max(as.numeric(ind))
