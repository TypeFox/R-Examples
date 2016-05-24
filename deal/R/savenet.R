## savenet.R --- 
## Author          : Claus Dethlefsen
## Created On      : Thu Sep 26 15:19:02 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Sep 28 13:34:12 2006
## Update Count    : 97
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bottcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

savenet <- function(nw,con=file("default.net")) {
    ## save network to .net file that can be read by eg. Hugin

    open(con,"w")
    cat("%Created by deal,",date(),"\n",file=con) # create empty file
    cat("%deal is Copyright (C) 2002-2006  Susanne Gammelgaard Bottcher, Claus Dethlefsen\n",file=con)
    cat(rep("%",60),"\n\n",sep="",file=con)
    ## ########################################
    ## Global information
    ## ########################################
    
    cat("net\n",file=con)
    cat("{\n",file=con)
    
    cat("\tnode_size = (40 40);\n",file=con)
    
    cat("}\n\n",file=con)
    
    ## ########################################
    ## DEFINE NODES
    ## ########################################
    for (i in 1:nw$n) { ## for each node
        nd <- nw$nodes[[i]]
        cat(
            nd$type,
            "node",
            nd$name,
            "\n",
            file=con)
        cat("{\n",file=con)
        if (nd$type=="discrete") {
            cat("\tstates = (",
                paste("\"",nd$levelnames,"\"",sep=""),
                ");\n",
                file=con)
        }
        cat("\tlabel = \"", nd$name,"\";\n",sep="",file=con)
        
        cat("\tposition = (",
            nd$position,
            ");\n",file=con)
        
        
        cat("}\n\n",file=con)
        
    }
    ## ########################################
    ## DEFINE POTENTIALS
    ## ########################################
    
    for (i in 1:nw$n) {
        nd <- nw$nodes[[i]]
        
        cat("\npotential (",
            nd$name, file=con)
        
        if (length(nd$parents)>0) {
            cat(" | ",file=con)
            ##            for (j in nd$parents)
            ##                cat(nw$nodes[[j]]$name," ",file=con)
            ## apparently, discrete parents must appear before cont.
            for (j in intersect(nd$parents,nw$discrete))
                cat(nw$nodes[[j]]$name," ",file=con)
            for (j in intersect(nd$parents,nw$continuous))
                cat(nw$nodes[[j]]$name," ",file=con)
            
        }
        cat(" )\n",file=con)
        
        cat("{\n",file=con)
        
        ## # parameters defining local distribution

        ## ##################################################################
        ## discrete nodes
        ## ##################################################################
        if (nd$type=="discrete") {

            cat("\tdata=(",file=con)
            
            
            ## the distribution of nd|parents in row-major layout
            if (length(nd$parents)>0) {

                cat("\n\t",file=con)
                
                if (FALSE) {
                    cat("nd$prob\n")
                    print(nd$prob)
                }
                
                dpar <- intersect(nd$parents,nw$discrete)
                ## Determine the discrete parents and their dimensions
                ## include (node$levels) as first component
                ## Dim <- c(nd$levels)
                Dim <- c()                
                for (i in dpar) {
                    Dim <- c(Dim, nw$nodes[[i]]$levels)
                }
                TD <- prod(Dim)
                
                ## dan alle teksterne i den rigtige rækkefolge
                lablist <- c()
                for (i in 1:TD) {
                    cf <- findex( i, Dim, FALSE)
##                    label <- nd$levelnames[cf[1,1]]
                    label <- ""
##                    for (j in 1:(ncol(cf)-1)) {
                    for (j in 1:ncol(cf)) {
##                        label <- paste(label, nw$nodes[[dpar[j]]]$levelnames[cf[1,j+1]]
                        label <- paste(label,
                                       nw$nodes[[dpar[j]]]$levelnames[cf[1,j]],
                                       sep=":")
                    }
                    lablist <- c(lablist,label)
                }
                
                ##we need to transform our column major mode
                ##to row major mode (used in .net files)
                
                cmajor <- array(1:TD,Dim)
                rmajor <- array(NA,Dim)
                for (i in 1:TD) {
                    lD <- length(Dim)
                    cf <- findex(i,Dim,config=FALSE)
                    a  <- c(1,cumprod(Dim[lD:1]))[lD:1]
                    idx <- sum(a*(cf-1))+1
                    rmajor[cf] <- idx
                }
                ##                        rmajor <- array(1:TD,Dim[nDim])
                ##                        rmajor <- aperm(rmajor,nDim)
                if (FALSE) {
                    cat("cmajor:\n");print(cmajor)
                    cat("rmajor:\n");print(rmajor)
                }
            
                ## write distribution for each config of disc. parents
                for (j in 1:TD) {
                    
                    ## transform from cmajor to rmajor
                    if (length(dpar)>1)
                        i <- cmajor[rmajor==j]
                    else
                        i <- j

                    cf <- findex(i,Dim,config=FALSE)
                    cfm <- cbind(1:nd$levels,
                                 matrix( rep(cf,nd$levels),
                                        ncol=length(cf), byrow=TRUE))
                    idx <- findex(cfm,c(nd$levels,Dim),config=TRUE)
                    if (FALSE) {
                    cat("cf=\n");print(cf)
                    cat("cfm=\n");print(cfm)                    
                    cat("Dim=\n");print(Dim)
                    cat("c(nd$levels,Dim)=\n");print(c(nd$levels,Dim))
                    cat("idx=\n");print(idx)
                }
                    cat(nd$prob[idx],"\t",file=con)
                    cat("\t%",lablist[i],"\n\t",file=con)
                }
                
            }
            else { 
                cat(nd$prob,file=con)
            }
            cat(");\n",file=con)
            
            
        } # discrete
    
####################################################################
    ## continuous nodes
####################################################################
    if (nd$type=="continuous") {
        cat("\tdata=(\n\t",file=con)
        
        
        ## skal skelne mellem kont. og disk. forældre
        ## hvis rene kont. forældre er prob en vektor
        if (length(nd$parents)>0) { # we have parents!
            if (length(intersect(nd$parents,nw$discrete))>0) {
                ## we have discrete parents
                
                dpar <- intersect(nd$parents,nw$discrete)
                ## Determine the discrete parents and their dimensions
                Dim <- c()
                for (i in dpar) {
                    Dim <- c(Dim, nw$nodes[[i]]$levels)
                    }
                TD <- prod(Dim)
                
                ## dan alle teksterne i den rigtige rækkefolge
                lablist <- c()
                for (i in 1:TD) {
                    cf <- findex( i, Dim, FALSE)
                    label <- ""
                    for (j in 1:ncol(cf)) {
                        label <- paste(label, nw$nodes[[dpar[j]]]$levelnames[cf[1,j]]
                                       ,sep=":")
                    }
                    lablist <- c(lablist,label)
                    }
                
                    ##we need to transform our column major mode
                    ##to row major mode (used in .net files)
                    if (length(dpar)>1) {
                        nDim <- c(2,1)
                        if (length(dpar)>2)
#                            nDim <- c(nDim,3:length(Dim[-c(1,2)]))
                            nDim <- c(nDim,3:length(Dim))

                        if (FALSE) {
                        cat("Dim:",Dim,"\n")
                        cat("nDim:",nDim,"\n")
                    }


                        ## ny strategi: udfyld rmajor paa en anden maade                        
                        cmajor <- array(1:TD,Dim)
                        rmajor <- array(NA,Dim)
                        for (i in 1:TD) {
                            lD <- length(Dim)
                            cf <- findex(i,Dim,config=FALSE)
                            a  <- c(1,cumprod(Dim[lD:1]))[lD:1]
                            idx <- sum(a*(cf-1))+1
                            rmajor[cf] <- idx
                        }
#                        rmajor <- array(1:TD,Dim[nDim])
#                        rmajor <- aperm(rmajor,nDim)
                        if (FALSE) {
                            cat("cmajor:\n");print(cmajor)
                            cat("rmajor:\n");print(rmajor)                            }
                    }
                    ## write distribution for each config of disc. parents
                    for (j in 1:TD) {

                        ## transform from cmajor to rmajor
                        if (length(dpar)>1)
                            i <- cmajor[rmajor==j]
                        else
                            i <- j
                        
                        cat("\tnormal ( ",nd$prob[i,2],file=con)
                        
                        if (length(nd$prob[i,])>2) { #cont.parents
                            for (j in 1:(length(nd$prob[i,])-2)) {
                                if (nd$prob[i,j+2]>=0)
                                    cat("+",file=con)
                                cat(nd$prob[i,j+2],"*",
                                    nw$nodes[[(intersect(nd$parents,nw$continuous))[j]]]$name,file=con)
                            }
                        }
                        ## print remark in file with the config of disc.par.
                        cat(", ",nd$prob[i,1],")","\t%",lablist[i],"\n",file=con)
                    }
                }
                else {
                    cat("normal ( ",nd$prob[2],file=con)
                    
                    for (j in 1:(length(nd$prob)-2)) {
                        if (nd$prob[j+2]>=0)
                            cat("+",file=con)                  
                        cat(nd$prob[j+2],"*",
                            nw$nodes[[(intersect(nd$parents,nw$continuous))[j]]]$name,file=con)
                            }
                    cat(", ",nd$prob[1],")\n",file=con)                    
                }
            }
            else {
                cat("normal ( ",nd$prob[2],", ",nd$prob[1],")\n",file=con)
            }
            
            cat("\t);\n",file=con)
        }
        
        cat("}\n",file=con)
                
    }
    
#    cat("Connection",filename,"created\n")
    close(con)
    invisible()
}
