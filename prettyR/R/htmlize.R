CreateIndexFile<-function(HTMLbase,HTMLdir,title="R listing") {
 basecon<-file(paste(HTMLdir,"/",HTMLbase,".html",sep="",collapse=""),"w")
 cat("<html><head><title>",title,"</title></head>\n",file=basecon)
 cat(" <frameset cols=200,* border=1>\n",file=basecon)
 cat(" <frame src=",paste(HTMLbase,"_nav.html",sep="",collapse=""),
  file=basecon,sep="")
 cat("  name=\"nav\" scrolling=yes>\n",file=basecon)
 cat(" <frame src=",paste(HTMLbase,"_list.html",sep="",collapse=""),
  file=basecon,sep="")
 cat("  name=\"list\" scrolling=yes>\n </frameset>\n</html>",file=basecon)
 close(basecon)
}

AddNavItem<-function(Rcommand,navcon,listname,navindex) {
 navitem<-ifelse(nzchar(Rcommand)>20,
  paste(paste(unlist(strsplit(Rcommand,""))[1:18],sep="",collapse=""),
   "...",sep="",collapse=""),
  Rcommand)
 nametag<-paste("ni",navindex,sep="",collapse="")
 cat(" <a href=\"",listname,"#",nametag,"\" target=\"list\">\n",
  file=navcon,sep="")
 cat("  ",navitem,"<br>\n </a>\n",file=navcon,sep="")
 return(nametag)
}

BeginNav<-function(navcon,bgcolor="#dddddd") {
 cat("<html>\n <body bgcolor=\"",bgcolor,"\">\n",file=navcon,sep="")
}

StartList<-function(listcon,title="R listing",bgcolor="#dddddd",useCSS=NULL) {
 cat("<html>\n <head>\n  <title>",title,"  </title>\n",file=listcon)
 if(!is.null(useCSS))
  cat(" <link rel=\"stylesheet\" type=\"text/css\" href=\"",
   useCSS,"\">\n",sep="",file=listcon)
 cat(" </head>\n",file=listcon)
 cat(" <body bgcolor=\"",bgcolor,"\">\n",file=listcon,sep="")
 cat("  <center><h1>",title,"</h1></center>\n",file=listcon)
}

EndHTML<-function(con,ending=NULL) {
 if(!is.null(ending)) cat(ending,file=con)
 cat(" </body>\n</html>\n",file=con)
}

htmlize<-function(Rfile,HTMLbase,HTMLdir,title,
 bgcolor="#dddddd",echo=FALSE,do.nav=FALSE,useCSS=NULL,...) {

 if(missing(Rfile)) stop("Minimal usage: htmlize(Rfile,...)")
 # Is Rfile there?
 if (!file.exists(Rfile)) {
  stopstring<-paste("Can't find",Rfile,collapse="")
  stop(stopstring)
 }
 # read in the commands
 Rcon<-file(Rfile,"r")
 Rcommands<-readLines(Rcon)
 close(Rcon)
 # if there is no HTML base name, use the Rfile name
 if(missing(HTMLbase)) {
  HTMLbase<-unlist(strsplit(basename(Rfile),"\\."))
  HTMLbase<-HTMLbase[1:(length(HTMLbase)-1)]
 }
 if(missing(title)) {
  if(charmatch("#title~",Rcommands[1],0))
   title<-strsplit(Rcommands[1],"~")[[1]][2]
  else title<-paste("Listing of",HTMLbase)
 }
 # If there is no HTML directory, use the path on the Rfile
 if(missing(HTMLdir)) {
  if(length(grep("/",Rfile))) HTMLdir<-unlist(strsplit(Rfile,"/"))
  else HTMLdir<-unlist(strsplit(Rfile,"\\\\"))
  HTMLdir<-
  ifelse(length(HTMLdir)==1,".",file.path(HTMLdir[-length(HTMLdir)]))
 }
 if(do.nav) {
  CreateIndexFile(HTMLbase,HTMLdir,title)
  navcon<-file(paste(HTMLdir,"/",HTMLbase,"_nav.html",sep="",collapse=""),"w")
  BeginNav(navcon,bgcolor=bgcolor)
  #navname<-paste(HTMLbase,"_nav.html",sep="",collapse="")
  listname<-paste(HTMLdir,"/",HTMLbase,"_list.html",sep="",collapse="")
 }
 else {
 listname<-paste(HTMLdir,"/",HTMLbase,".html",sep="",collapse="")
 }
 listcon<-file(listname,"w")
 StartList(listcon,title=title,bgcolor=bgcolor,useCSS=useCSS)
 listname<-paste(HTMLbase,"_list.html",sep="",collapse="")
 sink(listcon)
 on.exit({
  sink(NULL,type="output");
  sink(NULL,type="message");
  if(do.nav) close(navcon);
  close(listcon)
 })
 forbidden <- c("connection","fifo","file","sink")
 GraphicDevices<-c("bitmap","bmp","jpeg","png","tiff")
 thiscommand<-""
 navindex<-1
 for(i in 1:length(Rcommands)) {
  if(echo) cat(Rcommands[i],"<br>\n",file=listcon)
  # check for a comment character
  commentpos<-grep("#",Rcommands[i],fixed=TRUE)
  # cut out the comment
  if(length(commentpos)) Rcommands[i]<-strsplit(Rcommands[i],"#")[[1]][1]
  commexp<-try(parse(text=Rcommands[i]),silent=TRUE)
  # if this line is a complete command and thiscommand is not
  if(is.expression(commexp) && nzchar(thiscommand) > 0)
   # add a semicolon in case we're in braces
   Rcommands[i]<-paste(Rcommands[i],";",sep="",collapse="")
  thiscommand<-paste(thiscommand,Rcommands[i],sep="",collapse="")
  commexp<-try(parse(text=thiscommand),silent=TRUE)
  # If try doesn't yield an expression, it must be an incomplete line
  if(is.expression(commexp)) {
   if(do.nav) {
    nametag<-AddNavItem(thiscommand,navcon,listname,navindex)
    cat("<a name=\"",nametag,"\">\n",file=listcon,sep="")
   }
   fname<-strsplit(thiscommand,"\\(")[[1]][1]
   dont<-fname %in% forbidden
   if(dont) thiscommand<-paste("#",thiscommand,sep="",collapse=" ")
   navindex<-navindex+1
   # if the function is assigning its value, peel off the assignment
   if(!is.na(charmatch("<-",fname))) fname<-strsplit(fname,"-")[[1]][2]
   # does this command open a graphics device?
   if(fname %in% GraphicDevices) {
    hasquote<-"\"" %in% unlist(strsplit(thiscommand,""))
    gclist<-unlist(strsplit(thiscommand,ifelse(hasquote,"[\"]","[(,)]")))
    if(!hasquote) gclist[2]<-get(gclist[2])
    cat("<img src=\"",gclist[2],"\"><p>\n",sep="")
    if(hasquote) thiscommand<-paste(gclist[1],"\"",HTMLdir,"/",gclist[2],"\"",
     paste(gclist[3:length(gclist)],sep="",collapse=""),sep="",collapse="")
    else thiscommand<-paste(gclist[1],"(\"",HTMLdir,"/",gclist[2],"\",",
     paste(gclist[3:length(gclist)],sep="",collapse=""),")",sep="",collapse="")
   }
   # If it's not forbidden, evaluate the command
   if(!dont) {
    cat("<pre>\n",file=listcon)
    eval(parse(text=thiscommand))
    cat("</pre>\n",file=listcon)
   }
   # Get ready for a new command
   thiscommand<-""
  }
 }
 EndHTML(listcon)
 if(do.nav) EndHTML(navcon)
}
