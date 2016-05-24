ex_to_string <- function(ex_file,comment_dontrun=T){
  
  pattern_count <- function(pattern,string){
    pat_match <- gregexpr(pattern,string)[[1]][1]!=-1
    length(pat_match[pat_match==T])
  }
  # could also use stringr::str_count(ex_file[i],"\\}")
  
  file_lines <- readLines(ex_file)
  dont_runs <- grep("\\\\dontrun\\{",file_lines)
  for(i in dont_runs){
    open_br <- 0
    it <- 0
    while(open_br>0 || it==0){
      open_br_line <- pattern_count("\\{",file_lines[i])
      closed_br_line <- pattern_count("\\}",file_lines[i])
      open_br <- open_br + open_br_line - closed_br_line
      if(it==0 || comment_dontrun) file_lines[i] <- paste("#",file_lines[i])
      i <-  i+1
      it <- it + 1 
    }
    file_lines[i-1] <- paste("#",file_lines[i-1])
  }
  #eval(parse(text=paste(file_lines,sep="\n")),parent.frame(n=2))
  return(paste(file_lines,sep="\n"))
}