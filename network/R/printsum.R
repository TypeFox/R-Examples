######################################################################
#
# printsum.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 7/05/11
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various routines for printing/summarizing 
# network class objects.
#
# Contents:
#
#   print.network
#   print.summary.network
#   summary.character
#   summary.network
#
######################################################################


# Printing for network class objects.
#
print.network<-function(x, matrix.type=which.matrix.type(x), 
                      mixingmatrices=FALSE, na.omit=TRUE, print.adj=FALSE, ...)
{
    cat(" Network attributes:\n")
    for(i in 1:length(x$gal)){
      if (names(x$gal)[i]=="n"){
          attributeName<-"vertices"
          attributeValue<-x$gal[[i]]
      } else {
          attributeName<-names(x$gal)[i]
          attributeValue<-x$gal[[i]]
      }
      if(is.network(attributeValue)){
        if(attributeName=="design"){
          cat(" ",attributeName,"=\n")
          cat("      total missing =",network.edgecount(attributeValue),"\n")
          cat("    percent missing =",network.density(attributeValue),"\n")
        }else{
          cat(" ",attributeName,":\n",sep="")
          if(is.discrete(attributeValue)){
            assign(paste(" ",attributeName),attributeValue)
            print(table(get(paste(" ",attributeName))))
            if(mixingmatrices){
              cat("\n","mixing matrix for ",attributeName,":\n",sep="")
              print(mixingmatrix(x,attributeName))
            }
         }else{
            print(summary(attributeValue))
         }
        }
      }else{
        if(attributeName!="mnext"){
           if(is.discrete(attributeValue)){
             #assign(paste(" ",attributeName),attributeValue)
             #print(table(get(paste(" ",attributeName))))
             print(table(attributeValue,dnn=paste('  ',attributeName,':',sep='')))
          }else{
             # for short attributes, just print out the values
             if((class(attributeValue)%in%c("factor","character","numeric", "logical","integer","double","NULL"))&&(length(attributeValue) < 10)){
               # handle NULL case because cat won't print NULL
               if (is.null(attributeValue)){
                 cat(" ",attributeName,"= NULL\n")
               } else {
                 cat(" ",attributeName,"=",attributeValue,"\n")
               }
             } else{
               # special handling for classes where summary would give messy or non-useful output
               # don't print summary for net obs period or active attributes
               if (attributeName=='net.obs.period' || grepl('.active$',attributeName) ){
                 cat("  ",attributeName,": (not shown)\n", sep="")
               } else if (class(attributeValue)%in%c("matrix")){
                 cat("  ",attributeName,": ",nrow(attributeValue),"x",ncol(attributeValue)," matrix\n", sep="")
                 
               } else {
               # default to printing out the summary for the attribute   
                 cat("  ",attributeName,":\n", sep="")
                 print(summary(attributeValue))
               }
             }
          }
        }
      }
    }
    cat("  total edges=",network.edgecount(x,na.omit=FALSE),"\n")
    cat("    missing edges=",network.naedgecount(x),"\n")
    cat("    non-missing edges=",network.edgecount(x,na.omit=TRUE),"\n")
    vna<-list.vertex.attributes(x)
    if(na.omit){
      vna<-vna[vna!="na"]
    }
    if(length(vna)==0){
      cat("\n","No vertex attributes","\n",sep="")
    }else{
      cat("\n","Vertex attribute names:","\n")
      cat("   ",vna,"\n")
    }
    # Print list of edge attributes, but only if there are not very many edges
    # because list.edge.attributes is expensive on large nets
    if(length(x$mel)<=1000){
      ena<-list.edge.attributes(x)
      if(na.omit){
        ena<-ena[ena!='na']
      }
      if(length(ena)==0){
        cat("\n","No edge attributes","\n",sep="")
      }else{
        cat("\n","Edge attribute names:","\n")
        cat("   ",ena,"\n")
      }
    } else {
      cat("\n","Edge attribute names not shown","\n")
    }
    
    
    #Print the adjacency structure, if desired
    if(print.adj){
      if(is.multiplex(x)&&(matrix.type=="adjacency"))
        matrix.type<-"edgelist"
      if(is.hyper(x))
        matrix.type<-"incidence"
      cat("\n",matrix.type,"matrix:\n")    
      if(network.edgecount(x)>0){
        mat<-as.matrix.network(x,matrix.type=matrix.type)
        attr(mat,"n")<-NULL            #Get rid of any extra attributes
        attr(mat,"vnames")<-NULL
        attr(mat,"bipartite")<-NULL
        print(mat)
      }else
        cat("Empty Graph\n")
    }
    invisible(x)
}


#Print method for summary.character
print.summary.character <- function(x, max.print=10, ...){
  x<-table(x)
  nam<-names(x)
  x<-as.vector(x)
  names(x)<-nam
  if(length(x) <= max.print){
    print(x)
  }else{
    ord<-order(as.vector(x),decreasing=TRUE)
    cat(paste("   the ",max.print," most common values are:\n",sep=""))
    print(x[ord][1:max.print])
  }
  invisible(x)
}


#Print method for summary.network
print.summary.network<-function(x, ...){
    #Pull any extra goodies from summary.network (stored in gal)
    na.omit<-x%n%"summary.na.omit"
    mixingmatrices<-x%n%"summary.mixingmatrices"
    print.adj<-x%n%"summary.print.adj"
    #Print the network-level attributes
    class(x)<-"network"
    cat("Network attributes:\n")
    for(i in 1:length(x$gal)){
      if (names(x$gal)[i]=="n"){
        attributeName<-"vertices"
        attributeValue<-x$gal[[i]]
      } else {
        attributeName<-names(x$gal)[i]
        attributeValue<-x$gal[[i]]
      }
      if(!(attributeName%in%c("mnext","summary.na.omit", "summary.mixingmatrices","summary.print.adj"))){
        if(is.network(attributeValue)){
          if(attributeName=="design"){
            cat(" ",attributeName,"=\n")
            cat("   total missing = ",network.edgecount(attributeValue),"\n", sep="")
            cat("   percent missing  =",network.density(attributeValue),"\n", sep="")
          }else{
            cat(" ",attributeName,"=\n")
            print(attributeValue)
          }
        }else{
          if(is.discrete(attributeValue)){
            assign(paste(" ",attributeName),attributeValue)
            print(table(get(paste(" ",attributeName))))
            if(mixingmatrices){
              cat("\n","mixing matrix for ",attributeName,":\n",sep="")
              print(mixingmatrix(x,attributeName))
            }
          }else{
            if((class(attributeValue)%in%c("factor","character","numeric", "logical","integer","double"))&& (length(attributeValue) < 10)){
              cat("  ",attributeName," = ",attributeValue,"\n",sep="")
            }else{
              cat("  ",attributeName,":\n", sep="")
              print(summary(attributeValue))
            }
          }
        }
      }
    }
    cat(" total edges =",network.edgecount(x,na.omit=FALSE),"\n")
    cat("   missing edges =",network.naedgecount(x),"\n")
    cat("   non-missing edges =",network.edgecount(x,na.omit=TRUE),"\n")
    cat(" density =",network.density(x),"\n")

    #Print the network-level attributes
    van<-list.vertex.attributes(x)
    if(na.omit){
      van<-van[van!="na"]
    }
    if(length(van)==0){
      cat("\n","No vertex attributes","\n",sep="")
    }else{
      cat("\nVertex attributes:\n")
      for (i in (1:length(van))){ 
        if(van[i]=="vertex.names"){
          cat("  vertex.names:\n")
          cat("   character valued attribute\n")
          cat("   ",sum(!is.na(network.vertex.names(x)))," valid vertex names\n",sep="")
        }else{
        cat("\n ",van[i],":\n",sep="")
        aaval<-get.vertex.attribute(x,van[i],unlist=FALSE)
        aaclass<-unique(sapply(aaval,class))
        aaclass<-aaclass[aaclass!="NULL"]
        if(length(aaclass)>1){
          cat("   mixed class attribute\n")
          cat("   ",sum(!sapply(aaval,is.null)),"values\n")
        }else if(aaclass%in%c("logical","numeric","character","list")){
          cat("   ",aaclass," valued attribute\n",sep="")
          aalen<-sapply(aaval,length)
          if(all(aalen<=1)&&(aaclass!="list")){
            cat("   attribute summary:\n")
            print(summary(unlist(aaval)))
            if(is.discrete(unlist(aaval))&&mixingmatrices){
              cat("   mixing matrix:\n")
              print(mixingmatrix(x,van[i]))
            }
          }else{
            cat("   uneven attribute lengths; length distribution is\n")
            print(table(aalen))
          }
        }else{
          cat("   ",aaclass," valued attribute\n",sep="")
          cat("   ",length(aaval)," values\n",sep="")
        }
      }
    }
  }

    #Print the edge-level attributes
    ean <- list.edge.attributes(x)
    if(na.omit){
      ean<-ean[ean!="na"]
    }
    if(length(ean)==0){
      cat("\n","No edge attributes","\n",sep="")
    }else{
      cat("\nEdge attributes:\n")
      for (i in (1:length(ean))){ 
        cat("\n ",ean[i],":\n",sep="")
        eaval<-get.edge.attribute(x$mel,ean[i],unlist=FALSE)
        eaclass<-unique(sapply(eaval,class))
        eaclass<-eaclass[eaclass!="NULL"]
        if(length(eaclass)>1){
          cat("   mixed class attribute\n")
          cat("   ",sum(!sapply(eaval,is.null)),"values\n")
        }else if(eaclass%in%c("logical","numeric","character","list")){
          cat("   ",eaclass," valued attribute\n",sep="")
          ealen<-sapply(eaval,length)
          if(all(ealen<=1)&&(eaclass!="list")){
            cat("   attribute summary:\n")
            print(summary(unlist(eaval)))
          }else{
            cat("   uneven attribute lengths; length distribution is\n")
            print(table(ealen))
          }
        }else{
          cat("   ",eaclass," valued attribute\n",sep="")
          cat("   ",length(eaval),"values\n",sep="")
        }
      }
    }
    
    #Print the adjacency structure
    if(print.adj){
      matrix.type=which.matrix.type(x)
      if(is.multiplex(x)&&(matrix.type=="adjacency"))
        matrix.type<-"edgelist"
      if(is.hyper(x))
        matrix.type<-"incidence"
      cat("\nNetwork ",matrix.type," matrix:\n",sep="")    
      if(network.edgecount(x)>0){
        mat<-as.matrix.network(x,matrix.type=matrix.type)
        attr(mat,"n")<-NULL            #Get rid of any extra attributes
        attr(mat,"vnames")<-NULL
        attr(mat,"bipartite")<-NULL
        print(mat)
      }else
        cat("Empty Graph\n")
    }
    invisible(x)
}


#An internal routine to handle summaries of characters
summary.character <- function(object, ...){
  class(object)<-c("summary.character",class(object))
  object
}


# Summaries of network objects
#
summary.network<-function(object, na.omit=TRUE, mixingmatrices=FALSE, print.adj=TRUE, ...){
  #Add printing parameters as network objects, and change the class
  object%n%"summary.na.omit"<-na.omit
  object%n%"summary.mixingmatrices"<-mixingmatrices
  object%n%"summary.print.adj"<-print.adj
  class(object)<-"summary.network"
  #Return the object
  object
}

