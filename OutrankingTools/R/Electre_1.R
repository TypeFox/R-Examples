Electre_1 <-
function(performanceMatrix,alternatives,criteria,criteriaWeights,minmaxcriteria,concordance_threshold=1,discordance_threshold=0){
        
        cat("\014")
        ####################################################################################################################################################
        #                                                                                                                                                  #
        # Copyright Michel Prombo, 2014                                                                                                                    #
        # Module : electre III                                                                                                                             #
        # Version 1.0                                                                                                                                      #
        #                                                                                                                                                  #
        # Contributors:                                                                                                                                    #
        #   Michel Prombo <michel.prombo@statec.etat.lu>                                                                                                   #                                                                                              #
        #                                                                                                                                                  #    
        # This software, ELECTRE III , is a package for the R OR software which  allows to use MCDA algorithms and methods.                                #
        #                                                                                                                                                  #
        # This software is governed by the terms and conditions of the Free and open-source licenses. You can                                              #
        # use, modify and/ or redistribute the software under the terms of .......                                                                         #
        #                                                                                                                                                  #
        #                                                                                                                                                  #
        # As a counterpart to the access to the source code and rights to copy,modify and redistribute granted by the license, users are provided only     #
        # with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability.  #
        #                                                                                                                                                  #        
        # In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the       #
        # software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also         #
        # therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore         #
        # encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or #
        # data to be ensured and, more generally, to use and operate it in the # same conditions as regards security.                                      #
        #                                                                                                                                                  #    
        # The fact that you are presently reading this means that you have had knowledge of the terms and conditions of the Free and open-source licenses  #
        # and that you accept its terms.                                                                                                                   #
        #                                                                                                                                                  #
        ####################################################################################################################################################
        #require("igraph")  
        
        
        ## *******************  check validity of the objects manipulated by the current function *************************** ##
        
        
        ## *******************************************************  checking the input data  ***********************************************************#
        
        # data is filtered, check for some data consistency
        
        # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
        
        if (is.null(dim(performanceMatrix))) 
            stop("less than 2 criteria or 2 alternatives")    
        
        ## check the input data
        
        if (!((is.matrix(performanceMatrix) || (is.data.frame(performanceMatrix))))) 
            stop("wrong performanceMatrix, should be a matrix or a data frame")
        
        if (!(is.vector(alternatives)))
            stop("alternatives should be a vector")
        if (!is.character(alternatives))	
            stop("alternatives should be a character vector")
        if(!(nrow(performanceMatrix)==length(alternatives)))
            stop("length of alternatives should be checked") 	
        
        
        if (!(is.vector(criteria)))
            stop("criteria should be a vector")
        if (!is.character(criteria))	
            stop("criteria should be a character vector") 
        if(!(ncol(performanceMatrix)==length(criteria)))
            stop("length of criteria should be checked") 	
        
        if (!(is.vector(minmaxcriteria)))
            stop("minmaxcriteria should be a vector")
        minmaxcriteria=tolower(minmaxcriteria)
        if(!(ncol(performanceMatrix)==length(minmaxcriteria)))
            stop("length of minmaxcriteria should be checked") 
        n=length(minmaxcriteria)	
        for (i in 1:n){
            if(!((minmaxcriteria[i]=='min') ||(minmaxcriteria[i]=='max'))){
                stop(" Vector minmaxcriteria must contain 'max' or 'min' ")
            }
        }
        
        if (!(is.vector(criteriaWeights)))
            stop("criteriaWeights should be a vector")
        if (!is.numeric(criteriaWeights))	
            stop("criteriaWeights should be a numeric vector") 
        if(!(ncol(performanceMatrix)==length(criteriaWeights)))
            stop("length of criteriaWeights should be checked") 
        
        ##     *************************************  End of checking the validity of the "inputs" **********************************************************####    
        
        
        #  *************************************************************   Variables transformation *********************************************************   #	 
        pm  	<- performanceMatrix
        av  	<-  as.list(alternatives)
        cv  	<- as.list(criteria)
        vp  	<- criteriaWeights
        mmv  <- minmaxcriteria
        s1 	<- concordance_threshold
        s2 	<- discordance_threshold
        
        if (!(is.na(match("min", mmv)))) {
            for(j in 1:ncol(pm)){
                valmax=max(pm[,j])
                if (mmv[j]=="min") {
                    for (i in 1:nrow(pm)){
                        pm[i,j]=valmax-pm[i,j]
                    }
                }
            }
            
        }
        
        
        #  ****************************************           Calcul matrice des concordance                ***********************************   #	  
        c1=0
        c2=0
        
        i=1
        n=nrow(pm)
        mc <- matrix (rep(0, n*n), n, n)
        diag(mc)=1
        for (i in 1:nrow(pm)){
            
            k=i+1
            while (k<=n){
                for (j in 1:ncol(pm)){
                    
                    tp=pm[i,j]-pm[k,j]
                    if (tp>= 0) {
                        c1=c1+vp[j]
                    }
                    if (tp<= 0) {
                        c2=c2+vp[j]
                    }
                }
                mc[k,i]=c1/sum(vp)
                mc[i,k]=c2/sum(vp)
                c1=0
                c2=0
                k=k+1
            }
            
            concordance=mc
        }
        
        #  *******************************                Calcul matrice des discordances              ***************************************   #
        
        max1=0
        max2=0
        
        md <- matrix (rep(0, n*n), n, n)
        
        sigma=0						# sigma is  the maximum difference  gi(b)-gi(a) on any criterion i
        for (s in 1:ncol(pm)){
            sigtemp=max(pm[,s])-min(pm[,s])
            if (sigtemp>= sigma){
                sigma=sigtemp
            }
        }
        
        for (i in 1:nrow(pm)){
            
            k=i+1
            while (k<=n){
                for (j in 1:ncol(pm)){	  
                    tp1=pm[k,j]-pm[i,j]
                    if (tp1>= max1) {
                        max1=tp1
                    }
                }
                for (j in 1:ncol(pm)){	
                    tp2=pm[i,j]-pm[k,j]
                    if (tp2>= max2) {
                        max2=tp2
                    }
                }
                md[k,i]=max1/sigma
                md[i,k]=max2/sigma
                max1=0
                max2=0
                k=k+1
            }
            
            discordance =md
        }
        
        
        matfiltre <- matrix (rep(0, n*n), n, n)
        
        for (i in 1:n){
            for (j in 1:n){
                if (i != j){
                    if((mc[i,j]>=s1)&&(md[i,j]<=s2)){
                        matfiltre[i,j]=1
                    }
                }
            }
            
            filtre=matfiltre
        }
        
        
        rownames(matfiltre)=av
        colnames(matfiltre)=av
        
        matgraph=t(matfiltre)
        g1<-graph.adjacency(matgraph); 
        g2=plot(g1)
        
        
        mc <- round(matrix(mc,ncol=n,nrow=n),digits=4)
        md <- round(matrix(md,ncol=n,nrow=n),digits=4)    
        rownames(pm)=av
        colnames(pm)=cv
        rownames(mc)=av
        colnames(mc)=av
        rownames(md)=av
        colnames(md)=av
        rownames(matfiltre)=av
        colnames(matfiltre)=av
        
        
        
        
        
        
        # prepare the output
        sink("result.txt")
        cat("------------------------------------------------------------------------------  Performance table    -------------------------------------------------------------------------","\n")
        cat(" ","\n")
        cat(" ","\n")
        print(pm)
        cat(" ","\n")
        cat(" ","\n")
        cat("--------------------------------------------------------------------------  Concordance  matrix    ---------------------------------------------------------------------------","\n")
        cat(" ","\n")
        cat(" ","\n")
        print(t(mc))
        cat(" ","\n")
        cat(" ","\n")
        cat("---------------------------------------------------------------------------  Discordance  matrix    ---------------------------------------------------------------------------","\n")
        cat(" ","\n")
        cat(" ","\n")
        print(t(md))
        cat(" ","\n")
        cat(" ","\n")
        cat("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------","\n")
        sink()
        
        mc=t(mc)
        md=t(md)  
        out <-  list( "Performance Matrix"=pm,"Concordance Matrix"=mc,"Discordance Matrix" =md)       
    }
