######################################################################
#
# fileio.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 5/2/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines relating to file I/O.
#
# Contents:
#   read.dot
#   read.nos
#   write.dl
#   write.nos
#
######################################################################


#read.dot - Import a file in Graphviz .dot format.  This code was contributed
#by Matthijs den Besten.
read.dot <- function(...) {
  lines <- readLines(...);
  body <- lines[grep("->", lines, fixed=TRUE)];

  nodePairs <- sub('^[[:space:]]+\"', '\"',
                   sub('\"[;[:space:]]+$', '\"',
                       unlist(strsplit(body, "->"))));
  nodeLists <- split(nodePairs,1:length(nodePairs)%%2);

  nodes <- unique(nodePairs);
  edges <- data.frame(orig=nodeLists[[2]], dest=nodeLists[[1]]);

  n <- length(nodes);
  graph <- matrix(0, n, n, dimnames=list(nodes, nodes));
                                        #for(i in 1:nrow(edges)) {
                                        #  edge <- edges[i,];
                                        #  graph[edge$orig,edge$dest] <- 1;
                                        #} // Did not work as intended.
  for(node in nodes) {
    graph[node,nodes%in%edges$dest[edges$orig==node]] <- 1;
  }
                                       
  return(graph);
}


#read.nos - Read an input file in Neo-OrgStat format.  At this time, only the 
#graph stack is read; any coloring information is ignored.
read.nos<-function(file,return.as.edgelist=FALSE){
   #Get the formatting information
   f<-sapply(readLines(file,n=2),strsplit," ")
   #Parse the formatting information
   m<-as.numeric((f[[1]])[1])
   n<-as.numeric((f[[2]])[1])
   o<-as.numeric((f[[2]])[2])
   #Read the input data
   dat<-scan(file,skip=3)
   #Unpack the data in the proper order
   gstack<-array(dim=c(m,n,o))
   for(i in 1:m)
      for(j in 1:n)
         for(k in 1:o)
            gstack[i,j,k]<-dat[(i-1)*n*o+(j-1)*o+k]   
   #Return the stack
   if(return.as.edgelist)
     as.edgelist.sna(gstack)
   else
     gstack
}


#write.dl - Write a graph or graph stack in DL format
write.dl<-function(x,file,vertex.lab=NULL,matrix.lab=NULL){
  x<-as.sociomatrix.sna(x)
  if(is.list(x))
    stop("DL format requires all graphs to be of identical order.")
  if(is.matrix(x))
    x<-array(x,dim=c(1,NROW(x),NCOL(x)))
  m<-dim(x)[1]
  n<-dim(x)[2]
  #Write the DL header
  cat("DL n = ",n,", nm = ",m,", format = edgelist1\n",sep="",file=file)
  #Write the labels
  if(is.null(vertex.lab))
    vertex.lab<-dimnames(x)[[2]]
  if(is.null(vertex.lab))
    vertex.lab<-1:n
  if(is.character(vertex.lab))
    vertex.lab<-paste("\"",vertex.lab,"\"",sep="")
  cat("labels:\n",file=file,append=TRUE)
  cat(paste(vertex.lab,collapse=","),"\n",sep="",file=file,append=TRUE)
  if(is.null(matrix.lab))
    matrix.lab<-dimnames(x)[[1]]
  if(is.null(matrix.lab))
    matrix.lab<-1:m
  if(is.character(matrix.lab))
    matrix.lab<-paste("\"",matrix.lab,"\"",sep="")
  cat("matrix labels:\n",file=file,append=TRUE)
  cat(paste(matrix.lab,sep="",collapse=","),"\n",sep="",file=file, append=TRUE)
  #Write the data
  cat("data:\n",file=file,append=TRUE)
  for(i in 1:m){
    edges<-x[i,,]                  #Obtain the matrix of edges
    edges[is.na(edges)]<-0
    edges<-edges!=0
    rn<-row(x[i,,])[edges]         #Get rows, columns, values
    cn<-col(x[i,,])[edges]
    val<-x[i,,][edges]
    if(sum(edges>0)){
      for(j in 1:length(rn))         #Write the edges
        cat(rn[j],cn[j],val[j],"\n",file=file,append=TRUE)
    }
    if(i<m)
      cat("!\n",file=file,append=TRUE)
  }
}


#write.nos - Write a graph or graph stack in Neo-OrgStat format
write.nos<-function(x,file,row.col=NULL,col.col=NULL){
  if(is.list(x)||(class(x)=="network"))
    x<-as.sociomatrix.sna(x)
  if(is.list(x))
    stop("NOS format requires all graphs to be of identical order.")
  if(is.matrix(x))
    x<-array(x,dim=c(1,NROW(x),NCOL(x)))
  m<-dim(x)[1]
  n<-dim(x)[2]
  o<-dim(x)[3]
  #Write NOS header
  cat(m,"\n",n," ",o,"\n",sep="",file=file)
  if(is.null(row.col))
    row.col<-rep(0,n)
  if(is.character(row.col))
    row.col<-paste("\"",row.col,"\"",sep="")
  if(is.null(col.col))
    col.col<-rep(0,o)
  if(is.character(col.col))
    col.col<-paste("\"",col.col,"\"",sep="")
  cat(paste(c(row.col,col.col),collapse=" "),"\n",sep="",file=file,append=TRUE)
  #Write the data
  for(i in 1:m){
    for(j in 1:n)
      cat(paste(x[i,j,],collapse=" "),"\n",sep="",file=file,append=TRUE)
  }
}
