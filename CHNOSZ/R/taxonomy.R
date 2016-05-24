# CHNOSZ/taxonomy.R
# functions to work with NCBI taxonomy files
# 20100311 jmd
# 20100613 getnames and getnodes now only use the 'scan' function
#   to read files (instead of 'read.table')

## read the entire taxonomy nodes file
getnodes <- function(taxdir) {
  # we only read the first five tab-separated fields
  n <- scan(paste(taxdir,"nodes.dmp",sep="/"),what=as.list(character(5)),
    quote="",flush=TRUE,sep="\t")
  # the information we want is in fields 1, 3 and 5
  # (fields 2 and 4 contain the vertical bar)
  id <- as.numeric(n[[1]])
  parent <- as.numeric(n[[3]])
  rank <- n[[5]]
  nodes <- data.frame(id,parent,rank)
  return(nodes)
}

## get the taxonomic rank of the selected node
getrank <- function(id,taxdir,nodes=NULL) { 
  if(is.null(nodes)) nodes <- getnodes(taxdir)
  return(as.character(nodes$rank[match(id,nodes$id)]))
}

# get the parent id of the selected node
# or higher ranking ancestor if specified 
parent <- function(id,taxdir,rank=NULL,nodes=NULL) {
  if(is.null(nodes)) nodes <- getnodes(taxdir)
  id.out <- numeric()
  for(i in 1:length(id)) {
    myid <- id[i]
    if(is.null(rank)) {
      # the immediate parent
      id.out <- c(id.out,nodes$parent[nodes$id==myid])
    } else {
      # parent, grandparent, etc.
      myrank <- getrank(myid,taxdir,nodes)
      if(is.na(myrank)) myid <- NA
      else while(myrank != rank) {
        myid <- parent(myid,taxdir,nodes=nodes)
        myrank <- getrank(myid,taxdir,nodes)
        if(myid==1) break
      }
      id.out <- c(id.out,myid)
    }
  }
  return(id.out)
}

## get the ids of all of the parents of the selected node
allparents <- function(id,taxdir,nodes=NULL) {
  if(is.null(nodes)) nodes <- getnodes(taxdir)
  thisid <- id
  while(thisid != 1) {
    thisid <- parent(thisid,taxdir=taxdir,nodes=nodes)
    id <- c(id,thisid)
  }
  return(id)
}

## read the entire taxonomy names file
getnames <- function(taxdir) {
  n <- scan(paste(taxdir,"names.dmp",sep="/"),
    what=as.list(character(7)),flush=TRUE,quote="",sep="\t")
  id <- as.numeric(n[[1]])
  name <- n[[3]]
  type <- n[[7]]
  names <- data.frame(id,name,type)
  return(names)
}

## get the scientific name(s) of the selected node(s)
sciname <- function(id,taxdir,names=NULL) {
  if(is.null(names)) names <- getnames(taxdir)
  # we're only interested in the scientific names 
  isci <- which(names$type=="scientific name")
  names <- names[isci,]
  # now identify the entries of interest
  names <- names[match(id,names$id),]
  n <- as.character(names$name)
  if(length(id)==1) n <- s2c(n,sep="\t")[1]
  return(n)
}
