##' @title Network construction
##' 
##' @description Construction a network from experimental data or integrated PPI database.
##' 
##' @param input A data frame containing the experimental data.
##' @param local.net Logical value, indicating whether to construct a network from experimental data (if \code{TRUE}) or not (if \code{FLASE}). Default value is \code{FALSE}.
##' @param node.attribute A data frame containing node attributes. Default value is \code{NULL}.
##' @param db Integrated PPI database, either \code{Biogrid} or \code{\link{HPRD}}.
##' @param species This parameter indicates the biological species to which analyzable PPI data is related; currently \code{human} for "Homo sapiens" and \code{ath} for "Arabidopsis thaliana" are available.
##' @param ID.type The ID type of the biological genes or proteins, possible values are \code{Entrez gene} and \code{Gene symbol} when \code{db} is \code{Biogrid}, or \code{Gene symbol} when \code{db} is \code{\link{HPRD}}.
##' @param hierarchy This parameter indicates how many hierarchy are included in the network, currently it can be \code{0}, \code{1} or \code{2}. Default value is \code{1}.
##' @return A network in igraph format.
##' @seealso \code{\link{construct_local}}, \code{\link{construct_nlocal}}
##' @export
##' @examples
##' ## Construction a local network.
##' local<-data.frame(1:5,2:6)
##' attribute<-data.frame(1:6,c(2.2,5.3,1.2,4.5,6.2,0.6))
##' net<-construction(input=local,local.net=TRUE,node.attribute=attribute)
##' ## Construction a network from the human HPRD database.
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)

construction<-function(input,local.net=FALSE,node.attribute=NULL,
                       db=c("Biogrid","HPRD"),species=c("human","ath"),
                       ID.type=c("Gene symbol","Entrez Gene"),
                       hierarchy=1)
{
  db <- match.arg(db)
    
  species <- match.arg(species)
  ID.type <- match.arg(ID.type)
    
  if(local.net){
    graph<-construct_local(input=input,node.attribute=node.attribute)
  }else{
    species<-match.arg(species)
    graph<-construct_nlocal(input=input,db=db,species=species,
                            ID.type=ID.type,hierarchy=hierarchy)
  }

  return(graph)
}

##' @title Local network construction
##' 
##' @description Construction a network from experimental data.
##' 
##' @param input A data frame containing the experimental data.
##' @param node.attribute A data frame containing node attributes.
##' @return A network in igraph format.
##' @seealso \code{\link{construction}}, \code{\link{construct_nlocal}}.
##' @export
##' @examples
##' local<-data.frame(1:5,2:6)
##' attribute<-data.frame(c(0.2,0.3,0.2,0.5,0.1))
##' net<-construct_local(input=local,node.attribute=attribute)
construct_local<-function(input,node.attribute)
{
  if(missing(input)){
    stop("Not a data frame")
  }
  if(!is.data.frame(input)){
    stop("Not a data frame")
  }
  
  graph<-NULL
  
  if(nrow(input)>1){
    graph<-graph.data.frame(input,directed=FALSE)
    ##  add node attribute
    if(!is.null(node.attribute)){
      graph<-add.vertex.attri(graph=graph,input=node.attribute)
    }
  }else{
    stop("No edges are imported")
  }
  
  return(graph)
}

##' @title Non-local network construction
##' 
##' @description Construction a network from integrated PPI database.
##' 
##' @param input A data frame containing the experimental data.
##' @param db Integrated PPI database, either \code{Biogrid} or \code{\link{HPRD}}.
##' @param species This parameter indicates the biological species to which analyzable PPI data is related; currently \code{human} for "Homo sapiens" and \code{ath} for "Arabidopsis thaliana" are available.
##' @param ID.type The ID type of the biological genes or proteins, possible values are \code{Entrez gene} and \code{Gene symbol} when \code{db} is \code{Biogrid}, or \code{Gene symbol} when \code{db} is \code{\link{HPRD}}.
##' @param hierarchy This parameter indicates how many hierarchy are included in the network, currently it can be \code{0}, \code{1} or \code{2}. Default value is \code{1}.
##' @return A network in igraph format.
##' @seealso \code{\link{construction}}, \code{\link{construct_local}}.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construct_nlocal(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
construct_nlocal<-function(input,db=c("Biogrid","HPRD"),species=c("human","ath"),
                           ID.type=c("Gene symbol","Entrez Gene"),
                           hierarchy=1)
{
  db <- match.arg(db)
  
  species <- match.arg(species)
  ID.type <- match.arg(ID.type)
  
  if(!missing(input)){
    if( !is.data.frame(input)){
      stop("Not a data frame")
    }
  }
  
  ##import the data
  if(db=="HPRD"){
    HPRD<-NULL
    data("HPRD", envir=environment())
    net<-graph.data.frame(HPRD[,c("Interactor.1.Gene.symbol","Interactor.2.Gene.symbol")],
                            directed=FALSE)       
  }else if(db=="Biogrid"){
    if(species=="human"){
      human<-NULL
    }else if(species=="ath"){
      ath<-NULL
    }
    data(list=species, envir=environment())
    if(ID.type=="Entrez Gene"){
      if(species=="human"){
        net<-graph.data.frame(human[,c("Entrez.Gene.Interactor.A","Entrez.Gene.Interactor.B")],
                              directed=FALSE)
      }else if(species=="ath"){
        net<-graph.data.frame(ath[,c("Entrez.Gene.Interactor.A","Entrez.Gene.Interactor.B")],
                              directed=FALSE)
      }
    }else if(ID.type=="Gene symbol"){
      if(species=="human"){
        net<-graph.data.frame(human[,c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")],
                              directed=FALSE)
      }else if(species=="ath"){
        net<-graph.data.frame(ath[,c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")],
                              directed=FALSE)
      }
    }
  }
  
  if(missing(input)){
    return(net)
  }
  
  ##  match
  index<-match(input[,1],V(net)$name)
  index<-index[!is.na(index)]
  ##  create a  sub network
  graph<-induced.subgraph(graph=net,unlist(neighborhood(graph=net,order=hierarchy,nodes=index)))
  ##  vertex.hierarchy
  V(graph)$vertex.hierarchy<-rep(hierarchy,vcount(graph))
  if(hierarchy>=1){
    for(i in (hierarchy-1):0){
      index<-unlist(neighborhood(graph=graph,order=i,nodes=which(V(graph)$name %in% input[,1])))
      V(graph)$vertex.hierarchy[index]<-i
    }
  }
  ##  add vertex attributes
  graph<-add.vertex.attri(graph=graph,input=input)
  return(graph)
}

## Add attributes to the vertex. 
add.vertex.attri<-function(graph=NULL,input)
{
  if(!is.igraph(graph)){
    stop("Not an igraph object")
  }
  if(!is.data.frame(input)){
    input<-as.data.frame(input)
  }
  coln<-ncol(input)
  if(coln>1){
    colname<-colnames(input)
    index<-match(V(graph)$name,input[,1],nomatch=0)
    for(i in 2:coln){
      graph<-set.vertex.attribute(graph,colname[i],index=as.character(input[index,1]),
                                  value=input[index,i])
    }
  }
  return(graph)
}
