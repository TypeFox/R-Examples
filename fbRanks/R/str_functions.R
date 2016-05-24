str_remove.nonascii.2=function(string,sub=NULL){
  #sub is what to replace the non-ascii character with
  x=string
  asc <- iconv(x, "latin1", "ASCII")
  ind <- is.na(asc) | asc != x
    tmp.fun=function(x){
      tmp=iconv(x, "latin1", "ASCII", sub = "byte")
      for(i in 1:length(tmp)){
      starts=str_locate_all(tmp[i],"[<]")[[1]][,1]
      ends=str_locate_all(tmp[i],"[>]")[[1]][,1]
      if(is.null(sub)){ tmp[i]=str_remove(tmp[i], starts, ends)
      }else{ tmp[i]=str_remove(tmp[i], starts, ends, sub=sub) }
      }
      tmp
      }
  for(x.i in 1:length(x)){
    if(ind[x.i]){ x[x.i]=tmp.fun(x[x.i]) }
  } 
  x
}

str_remove.nonascii=function(string,sub=""){
  #sub is what to replace the non-ascii character with
  x=string
  asc <- iconv(x, "latin1", "ASCII")
  ind <- is.na(asc) | asc != x
  for(x.i in 1:length(x)){
    if(ind[x.i]){ x[x.i]=iconv( x[x.i], "latin1", "ASCII", sub = sub ) }
  } 
  x
}

str_remove=function(string, start = 1L, end = -1L, sub=""){
    rep.char=c("&","@","#","~","%") #surely one of these will be missing
    for(i in rep.char){
      if(all(!str_detect(string, i))){ repc=i; break }
    }
    if(length(string)!=1){
      start = rep(start,length.out=length(string))
      end = rep(end,length.out=length(string))
    }
    for(i in 1:length(start)){
    str_sub(string[i], start[i], end[i]) = paste(rep(repc,end[i]-start[i]+1),collapse="")
    }
    str_replace_all(string,repc,sub)
  }

str_proper=function(string) {
  for(i in 1:length(string)){
    x=tolower(string[i])
  s <- strsplit(x, " ")[[1]]
  string[i]=paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
  }
  string
}

str_strip.white=function(string, sub=" ") 
{
  pattern <- "^\\s+|\\s+$"
  string <- str_replace_all(string, pattern, "")
  pattern <- "\\s+"
  str_replace_all(string, pattern, sub)
}
