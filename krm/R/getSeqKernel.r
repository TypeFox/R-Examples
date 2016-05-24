getSeqKernel <- function(sequences, kern.type=c("mm", "prop", "mi"), tau, call.C=TRUE, seq.start=NULL, seq.end=NULL) {
    
    kern.type <- match.arg(kern.type)    
    
    if (is.character(sequences)) {
        # sequences is a file name
        seq.file.name=sequences
        seq.alignment <- readFastaFile(seq.file.name)
        
    } else if (is.list(sequences)) {
        # sequences is a list of strings
        if (is.null(names(sequences))) names(sequences)=1:length(sequences)
        seq.alignment = sequences
        # get substr
        if (!is.null(seq.start) & !is.null(seq.end)) {
            for (i in 1:length(seq.alignment)) {
                seq.alignment[[i]]=substr(seq.alignment[[i]], seq.start, seq.end)
            }
        }
    } else stop("sequences input incorrect, should be either sequence file name or a list of strings")
    
    N=length(seq.alignment)
    I=length(strsplit(seq.alignment[[1]], split="")[[1]])
    
#   sum(strsplit(c1.data$seq[[1]],split="")[[1]]!=strsplit(c1.data$seq[[2]],split="")[[1]])
    
    #### 20 amino acids plus gap ####
    if (kern.type == "mm") {
        
        aa.list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
        K=length(aa.list)
    
#       seq.array <- array(0, dim=c(N,I,K), dimnames=list(seq.alignment[[2]], NULL, aa.list))
        seq.array <- array(0, dim=c(N,I,K), dimnames=list(names(seq.alignment), NULL, aa.list))
        seq.matrix <- matrix(0, nrow=N, ncol=(I*K))
        for (n in 1:N) {
            for (i in 1:I) {
#               seq.array[n,i,strsplit(seq.alignment$seq[[n]],split="")[[1]][i]]=1
                seq.array[n,i,strsplit(seq.alignment[[n]],split="")[[1]][i]]=1
                seq.matrix[n,(((i - 1) * K + 1):((i - 1) * K + K))]=seq.array[n,i,] 
            }
        }
        rm(seq.array)
        
        mismatch.dist <- matrix(0, nrow=N, ncol=N)
        for (j in 1:N) {
            for (k in 1:N) {
                mismatch.dist[j,k]=sum((seq.matrix[j,] - seq.matrix[k,])^2)
            }
        }
        return(exp(-tau*(mismatch.dist/2)))
        
    } else if (kern.type == "prop") {    
    
#       factanal(krm::aa.prop.list[,4:13], factors=5) suggests 5 factors(properties) are enough
        K=7       
        prop.final <- matrix(0, nrow=20, ncol=K,dimnames=list(krm::aa.prop.list[,1],c("Surface_Area_Chothia","Hydrophobicity_Engleman","Refractivity_Jones","Polarity_Jones","Residue_Volume_Zamayatin","Bulkiness_Jones","Hydrophilicity_Hopp")))
        prop.final[,"Surface_Area_Chothia"]=krm::aa.prop.list[,"Surface_Area_Chothia"]/max(abs(krm::aa.prop.list[,"Surface_Area_Chothia"]))
        prop.final[,"Refractivity_Jones"]=krm::aa.prop.list[,"Refractivity_Jones"]/max(abs(krm::aa.prop.list[,"Refractivity_Jones"]))
        prop.final[,"Hydrophobicity_Engleman"]=krm::aa.prop.list[,"Hydrophobicity_Engleman"]/max(abs(krm::aa.prop.list[,"Hydrophobicity_Engleman"]))
        prop.final[,"Polarity_Jones"]=krm::aa.prop.list[,"Polarity_Jones"]/max(abs(krm::aa.prop.list[,"Polarity_Jones"]))
        prop.final[,"Residue_Volume_Zamayatin"]=krm::aa.prop.list[,"Residue_Volume_Zamayatin"]/max(abs(krm::aa.prop.list[,"Residue_Volume_Zamayatin"]))
        prop.final[,"Bulkiness_Jones"]=krm::aa.prop.list[,"Bulkiness_Jones"]/max(abs(krm::aa.prop.list[,"Bulkiness_Jones"]))
        prop.final[,"Hydrophilicity_Hopp"]=krm::aa.prop.list[,"Hydrophilicity_Hopp"]/max(abs(krm::aa.prop.list[,"Hydrophilicity_Hopp"]))
        seq.matrix <- matrix(0, nrow=N, ncol=(I*K))
        for (n in 1:N) {
            for (i in 1:I) {
#               if(toupper(strsplit(seq.alignment$seq[[n]],split="")[[1]][i]) == "-") {
                if(strsplit(seq.alignment[[n]],split="")[[1]][i] == "-") {    
                    seq.matrix[n,(((i - 1) * K + 1):((i - 1) * K + K))]=rep(10000, times=K)  ### putting a spurious large value for a gap, the max distance otherwise should be 1
                } else {
#                   seq.matrix[n,(((i - 1) * K + 1):((i - 1) * K + K))]=prop.final[toupper(strsplit(seq.alignment$seq[[n]],split="")[[1]][i]),]
                    seq.matrix[n,(((i - 1) * K + 1):((i - 1) * K + K))]=prop.final[strsplit(seq.alignment[[n]],split="")[[1]][i],]
                }
            }
        }
        
        ###  Calculate distance here
        prop.dist <- matrix(0, nrow=N, ncol=N)
        for (j in 1:N) {
            for (k in 1:N) {
                gap.length=0
                gap.init=FALSE
                for (i in 1:(I*K)) {
                     if(((seq.matrix[j,i] == 10000)&&(seq.matrix[k,i] != 10000))||((seq.matrix[k,i] == 10000)&&(seq.matrix[j,i] != 10000))) {
                         gap.init=TRUE
                         gap.length=gap.length + 1
                         if ((gap.init == TRUE) && (gap.length > 1)) {
                             prop.dist[j,k]=prop.dist[j,k] + 1^2 ### using a gap.extension penalty of 1 
                         } else {
                             prop.dist[j,k]=prop.dist[j,k] + 3^2 ### using a gap.initiation penalty of 3 
                         } 
                     } else {
                         gap.init=FALSE
                         prop.dist[j,k]=prop.dist[j,k] + ((seq.matrix[j,i] - seq.matrix[k,i])^2)
                     }
                }
            }
        }
        return(exp(-tau*prop.dist))        
        
    } else if (kern.type == "mi") {
                       
        if (!call.C) {
            
            aaPrior=krm::cloud9
            dat <- string2arabic(seq.alignment)
            N=nrow(dat)
            
            kernel.dist <- matrix(0, nrow=N, ncol=N)
            for (j in 1:(N-1)) {
                for (k in (j+1):N) {
                    num.loglik=hmmMargLlik(dat[c(j,k),,drop=FALSE], aaPrior, 1)
                    denom.loglik=0.5*( hmmMargLlik(dat[c(j,j),,drop=FALSE], aaPrior, 1) + hmmMargLlik(dat[c(k,k),,drop=FALSE], aaPrior, 1) )
                    kernel.dist[j,k]=exp(tau*(num.loglik - denom.loglik))               
                }
            }
            for (j in 1:(N-1)) for (k in (j+1):N) kernel.dist[k,j]=kernel.dist[j,k]
            for (j in 1:N) kernel.dist[j,j]=1
            
        } else {
            
            N=length(seq.alignment)
            kernel.dist <- matrix(0.0, nrow=N, ncol=N)
            if (!is.double(kernel.dist)) kernel.dist <- as.double(kernel.dist)
            
            #aux=.C("MI_kernel", dataFileCh=as.character(seq.file.name), dataFileLen=nchar(seq.file.name), tau=as.double(tau), K=kernel.dist)
            aux=.C("MI_kernel_str", seq=concatList(seq.alignment,sep=""), nSeq=N, seqLength=I, tau=as.double(tau), K=kernel.dist)
            kernel.dist=aux$K
        
        } # end if !call.C
            
        return(kernel.dist)
        
    } # end if kern.type
    
    
} # end function
