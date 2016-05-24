# Automatically generated from all.nw using noweb
#$Log: printBanner.q,v $
#Revision 1.4  2007/01/23 21:00:27  sinnwell
#rm ending newline \n.  Users can space if desired.
#
#Revision 1.3  2005/02/04 20:57:18  sinnwell
#banner.width now based on options()$width
#char.perline based on banner.width
#
#Revision 1.2  2004/06/25 15:56:48  sinnwell
#now compatible with R, changed end when a line is done
#
#Revision 1.1  2004/02/26 21:34:55  sinnwell
#Initial revision
#

printBanner <- function(str, banner.width=options()$width, char.perline=.75*banner.width, border = "="){

# char.perline was calculated taking the floor of banner.width/3

  vec <- str
  new<-NULL
  onespace<-FALSE
  for(i in 1:nchar(vec)){
    if (substring(vec,i,i)==' ' && onespace==FALSE){
      onespace<-TRUE
      new<-paste(new,substring(vec,i,i),sep="")}
    else if (substring(vec,i,i)==' ' && onespace==TRUE)
      {onespace<-TRUE}
    else{
      onespace<-FALSE
      new<-paste(new,substring(vec,i,i),sep="")}
  }
  
  where.blank<-NULL
  indx <- 1
  
  for(i in 1:nchar(new)){
    if((substring(new,i,i)==' ')){
      where.blank[indx]<-i
      indx <- indx+1
    }
  }
  

# Determine the position in the where.blank vector to insert the Nth character position of "new"
  j<-length(where.blank)+1

# Add the Nth character position of the "new" string to the where.blank vector.
  where.blank[j]<-nchar(new)
  
  begin<-1
  end<-max(where.blank[where.blank<=char.perline])

# If end.ok equals NA then the char.perline is less than the position of the 1st blank.
  end.ok <- is.na(end) 

# Calculate a new char.perline. 
  if (end.ok==TRUE){ 
    char.perline <- floor(banner.width/2)
    end<-max(where.blank[where.blank<=char.perline])
  }

  cat(paste(rep(border, banner.width), collapse = ""),"\n")

  repeat {
    titleline<-substring(new,begin,end)
    n <- nchar(titleline)
    if(n < banner.width)
      {
        n.remain <- banner.width - n
        n.left <- floor(n.remain/2)
        n.right <- n.remain - n.left
        for(i in 1:n.left) titleline <- paste(" ",titleline,sep="")
        for(i in 1:n.right) titleline <- paste(titleline," ",sep="")
        n <- nchar(titleline)
      }
     cat(titleline,"\n")
    begin<-end+1
    end.old <- end
   # Next line has a problem when used in R.  Use print.banner.R until fixed.
   # Does max with an NA argument
    tmp <- where.blank[(end.old<where.blank) & (where.blank<=end.old+char.perline+1)]
    if(length(tmp)) end <- max(tmp)
    else break
   
#   end<-max(where.blank[(end.old<where.blank)&(where.blank<=end.old+char.perline+1)])
#   end.ok <- is.na(end)
#   if (end.ok==TRUE)
#      break
  }
  
  cat(paste(rep(border, banner.width), collapse = ""), "\n")
  invisible()
  
}

