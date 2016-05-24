# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 28-05-2014
# last modification: 10-10-2014
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# compile some Rdata files into one Rdata file
#           OR
# compile a list of bootSPLS objects into one bootSPLS object

compile.bootsPLS.object=function(bootsPLS.list,path,pattern,file,save.file)
{
    
    # --------if compiling from a list of bootSPLS objects:
    # bootSPLS.list
    
    
    # --------if compiling from Rdata files and looking for the Rdata files:
    # path: where to look for the Rdata files
    # pattern: pattern of the Rdata files (as in list.files()
    
    # --------if compiling from Rdata files and the files are given:
    # file= vector of the files to be compiled
    
    # -------- in any case
    # H is the number of component
    # save.file=TRUE, should the concatenation of result be saved? if yes, the file is saved as pattern.Rdata
    
    if(missing(pattern) & missing(file) & missing(bootsPLS.list)) stop("missing arguments")
    if((!missing(pattern))&(!missing(file))&(!missing(bootsPLS.list)))  warning("Argument `bootSPLS.list' is used")
    
    if((!missing(pattern) | !missing(file)) & (missing(bootsPLS.list)))  from="file"
    if(!missing(bootsPLS.list)) from="object"

    
    if((!missing(pattern))&(!missing(file))) warning("Argument `file' is used")
    
    if((!missing(pattern))&(missing(file)))
    {
        file=list.files(pattern=pattern,path=path)
        cat("the files used are", file, "\n")
        
        if(length(file)==0) stop("no files found with the path/pattern")
        nbr.objects=length(file)
    }
    
    
    if(from=="file")
    {
        bootsPLS.list=list()
        for(i in 1:length(file))
        {
            load(paste(path,file[i],sep=""))
            out=list(ClassifResult=ClassifResult,loadings.X=loadings.X,selection.variable=selection.variable,
            frequency=frequency,nbr.var=nbr.var,learning.sample=learning.sample,prediction=prediction,data=data,nzv=nzv)
            class(out)="bootsPLS"
            
            bootsPLS.list[[i]]=out
        }
        
    }
    

    nbr.objects=length(bootsPLS.list)
    
    #first quick check
    for(i in 1:nbr.objects) if(class(bootsPLS.list[[i]])!="bootsPLS") stop("bootsPLS.list does not only contain `bootsPLS' object")


    for(i in 1:nbr.objects)
    {
        print(i)


        # number of replications
        nbr.replication=nrow(bootsPLS.list[[i]]$nbr.var)
        
        if(i==1)
        {
            
            p=ncol(bootsPLS.list[[i]]$data$X) #assume same data for all the files, should check that somewhere
            n=nrow(bootsPLS.list[[i]]$data$X)
            H=ncol(bootsPLS.list[[i]]$nbr.var)
            nlevelY=nlevels(bootsPLS.list[[i]]$data$Y)
            
            ClassifResult=bootsPLS.list[[i]]$ClassifResult
            loadings.X=bootsPLS.list[[i]]$loadings.X
            selection.variable=bootsPLS.list[[i]]$selection.variable
            nbr.var=bootsPLS.list[[i]]$nbr.var
            learning.sample=bootsPLS.list[[i]]$learning.sample
            prediction=bootsPLS.list[[i]]$prediction
            X=bootsPLS.list[[i]]$data$X
            Y=bootsPLS.list[[i]]$data$Y
            method=bootsPLS.list[[i]]$data$method
            data=list(X=X,Y=Y,method=method)
            nzv=bootsPLS.list[[i]]$nzv
            
            nbr.tot=nbr.replication
            
        }else{
            
            nbr.var.temp=bootsPLS.list[[i]]$nbr.var
            nbr.var=rbind(nbr.var,nbr.var.temp)

            ClassifResult.temp=array(0,c(nlevelY,nlevelY,H,nbr.tot+nbr.replication))
            ClassifResult.temp[,,,1:nbr.tot]=ClassifResult
            ClassifResult.temp[,,,(nbr.tot+1):(nbr.tot+nbr.replication)]=bootsPLS.list[[i]]$ClassifResult
            ClassifResult=ClassifResult.temp
            
            loadings.X.temp=array(0,c(p,H,nbr.tot+nbr.replication))
            loadings.X.temp[,,1:nbr.tot]=loadings.X
            loadings.X.temp[,,(nbr.tot+1):(nbr.tot+nbr.replication)]=bootsPLS.list[[i]]$loadings.X
            loadings.X=loadings.X.temp
            
            selection.variable.temp=array(0,c(H,p,nbr.tot+nbr.replication))
            selection.variable.temp[,,1:nbr.tot]=selection.variable
            selection.variable.temp[,,(nbr.tot+1):(nbr.tot+nbr.replication)]=bootsPLS.list[[i]]$selection.variable
            selection.variable=selection.variable.temp
            
            
            learning.sample.temp=matrix(0,nrow=n,ncol=nbr.tot+nbr.replication) #record which sample are in the learning set
            learning.sample.temp[,1:nbr.tot]=learning.sample
            learning.sample.temp[,(nbr.tot+1):(nbr.tot+nbr.replication)]=bootsPLS.list[[i]]$learning.sample
            learning.sample=learning.sample.temp
            
            prediction.temp=array(0,c(n,nbr.tot+nbr.replication,H)) #record the class associated to each sample (either in learning or test set)
            prediction.temp[,1:nbr.tot,]=prediction
            prediction.temp[,(nbr.tot+1):(nbr.tot+nbr.replication),]=bootsPLS.list[[i]]$prediction
            prediction=prediction.temp
            
            nbr.tot=nbr.tot+nbr.replication
            
        }
        
        
    }
    
    rownames(ClassifResult)=levels(Y)
    colnames(ClassifResult)=colnames(bootsPLS.list[[1]]$ClassifResult)
    dimnames(ClassifResult)[[3]]=dimnames(bootsPLS.list[[1]]$ClassifResult)[[3]]
    dimnames(ClassifResult)[[4]]=paste("iteration.",1:nbr.tot,sep="")
    dimnames(selection.variable)[[2]]=dimnames(bootsPLS.list[[1]]$selection.variable)[[2]]
    dimnames(loadings.X)[[1]]=dimnames(bootsPLS.list[[1]]$selection.variable)[[2]]


    
    #calculation of the frequency of selection for each variable, on each components
    frequency=matrix(0,nrow=H,ncol=dim(loadings.X)[1])
    for(j in 1:nbr.tot)
    {
        for(k in 1:H)
        {a=which(loadings.X[,k,j]!=0)
            frequency[k,a]=frequency[k,a]+1 #add 1 everytime the gene is selected
        }
    }
    frequency=frequency/nbr.tot #get the probability of selection (percentage of times each gene is selected, per component
    colnames(frequency)=colnames(bootsPLS.list[[1]]$frequency)
    

    
    if(!missing(save.file))
    save(ClassifResult,loadings.X,selection.variable,frequency,nbr.var,learning.sample,prediction,data=data,nzv,file=save.file)


    out=list(ClassifResult=ClassifResult,loadings.X=loadings.X,selection.variable=selection.variable,frequency=frequency,
    nbr.var=nbr.var,learning.sample=learning.sample,prediction=prediction,data=data,nzv=nzv)
    structure(out,class="bootsPLS")


}
