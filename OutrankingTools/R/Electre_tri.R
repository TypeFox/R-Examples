Electre_tri <-
function(performanceMatrix,alternatives,profiles,profiles_names,criteria,minmaxcriteria,criteriaWeights,IndifferenceThresholds,PreferenceThresholds,VetoThresholds,lambda=NULL){
        
        cat("\014") 
        ####################################################################################################################################################
        #                                                                                                                                                  #
        # Copyright Michel Prombo, 2014                                                                                                                    #
        # Module : electre tri                                                                                                                             #
        # ELECTRE TRI assigns alternatives to predefined categories. The assignment of an alternative results from the comparison with the profiles        # 
        #defining the limits of the categories.                                                                                                            #
        # Version 1.0                                                                                                                                      # 
        #                                                                                                                                                  #
        # Contributors:                                                                                                                                    #
        #   Michel Prombo <michel.prombo@statec.etat.lu>            																	  #
        #   Kevin Prombo-Rosamont <kevinrosamont@ymail.com>  					               											  #
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
        
        ## *******************************************************  checking the input data  *****************************************************************#
        
        # data is filtered, check for some data consistency
        
        # if there are less than 2 criteria or 2 alternatives, there is no problem to resolve here
        
        if (is.null(dim(performanceMatrix))) 
            stop("less than 2 criteria or 2 alternatives")		
        
        
        if (!((is.matrix(performanceMatrix) || (is.data.frame(performanceMatrix))))) 
            stop("wrong performanceMatrix, should be a matrix or a data frame")
        
        if (!((is.matrix(profiles) || (is.data.frame(profiles))))) 
            stop("wrong profiles, should be a matrix or a data frame")
        if (is.null(dim(profiles))) 
            stop("less than 2 criteria or 2 profiles")	
        if (!(is.vector(profiles_names)))
            stop("profiles_names should be a vector")
        
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
        
        
        if (!(is.vector(IndifferenceThresholds)))
            stop("IndifferenceThresholds should be a vector")
        if(!(ncol(performanceMatrix)==length(IndifferenceThresholds)))
            stop("length of IndifferenceThresholds should be checked") 
        
        if (!(is.vector(PreferenceThresholds)))
            stop("PreferenceThresholds should be a vector")
        if(!(ncol(performanceMatrix)==length(PreferenceThresholds)))
            stop("length of PreferenceThresholds should be checked") 
        
        if (!(is.vector(VetoThresholds)))
            stop("VetoThresholds should be a vector")
        if(!(ncol(performanceMatrix)==length(VetoThresholds)))
            stop("length of VetoThresholds should be checked") 
        
        
        
        ##     *************************************  End of checking the validity of the "inputs" **********************************************************####
        
        pm=performanceMatrix
        pr_names=profiles_names
        nrow_pm=nrow(pm)
        al <- alternatives
        pr <- profiles
        nrow_pr=nrow(pr)
        cr <- criteria
        length_cr=length(cr)
        vp <- criteriaWeights
        if (sum(vp)>1){
            vp <- vp/sum(vp)
        }
        
        mmv <- minmaxcriteria
        p_t=PreferenceThresholds
        q_t=IndifferenceThresholds
        v_t=VetoThresholds
        
        if (is.null(lambda)){
            lambda <-0.75
        }
        
        
        partialConcordance_al_pr_gj  <-  matrix (rep(0, nrow_pm*length_cr), nrow_pm, length_cr)
        partialConcordance_pr_al_gj  <-  matrix (rep(0, nrow_pm*length_cr), nrow_pm, length_cr)
        partialDiscordance_al_pr_gj  <-  matrix (rep(0, nrow_pm*length_cr), nrow_pm, length_cr)
        partialDiscordance_pr_al_gj  <-  matrix (rep(0, nrow_pm*length_cr), nrow_pm, length_cr)
        globalconcordance      	    <-  matrix (rep(0, nrow_pm*nrow_pr*2), nrow_pm, nrow_pr*2)
        credibility     		    <-  matrix (rep(0, nrow_pm*nrow_pr*2), nrow_pm, nrow_pr*2)
        relations  	              <-  matrix (rep(0, nrow_pm*nrow_pr*2), nrow_pm, nrow_pr*2)
        Pessimistic 	     	    <-  matrix (rep(0, nrow_pm*(nrow_pr+1)), nrow_pm, nrow_pr+1)
        Optimistic 	     	    <-  matrix (rep(0, nrow_pm*(nrow_pr+1)), nrow_pm, nrow_pr+1)
        
        names_glo_con=list()
        for (i in 1:nrow_pr){
            v1=paste("c(","ai",",","b",i,")",sep="")
            names_glo_con <-union(names_glo_con,v1)
            v2=paste("c(","b",i,",","ai",")",sep="")
            names_glo_con <-union(names_glo_con,v2)
        }
        
        names_relations=list()
        for (i in 0:(2*nrow_pr-1)){
            v1=paste("b",i,sep="")
            names_relations <-union(names_relations,v1)
        }
        
        names_categories=list()
        for (i in 1:(nrow_pr+1)){
            v1=paste("C",i,sep="")
            names_categories <-union(names_categories,v1)
        }
        
        colnames(globalconcordance)=names_glo_con
        rownames(globalconcordance)=al
        colnames(credibility)=names_glo_con
        rownames(credibility)=al
        # colnames(relations)=names_glo_con
        rownames(relations)=al
        colnames(Pessimistic)=names_categories
        rownames(Pessimistic)=al
        colnames(Optimistic)=names_categories
        rownames(Optimistic)=al
        
        rownames(partialConcordance_al_pr_gj)=al
        colnames(partialConcordance_al_pr_gj)=cr
        rownames(partialConcordance_pr_al_gj)=al
        colnames(partialConcordance_pr_al_gj)=cr
        
        rownames(partialDiscordance_al_pr_gj)=al
        colnames(partialDiscordance_al_pr_gj)=cr
        rownames(partialDiscordance_pr_al_gj)=al
        colnames(partialDiscordance_pr_al_gj)=cr
        
        credibility_al_pr  <-  matrix (rep(0, nrow_pm%*%nrow_pr), nrow_pm, nrow_pr)
        credibility_pr_al  <-  matrix (rep(0, nrow_pr*nrow_pm), nrow_pr, nrow_pm)
        # globalconcordance_al_pr=0
        # globalconcordance_pr_al=0
        
        
        #   1. Computation of partial concordance indices  cj(a,b)  and cj(b,a)
        
        h=0
        k=1
        while (k<=nrow_pr){
            
            for (i in 1: nrow_pm){
                # Columns' iterations	
                
                for (j in 1:length_cr){
                    if (mmv[j]=="max"){
                        c_a_b = (pm[i,j]-pr[k,j]+p_t[j])/(p_t[j]-q_t[j])
                    }
                    if (mmv[j]=="min"){
                        c_a_b = (pr[k,j]-pm[i,j]+p_t[j])/(p_t[j]-q_t[j])	
                    }
                    c_a_b=round(c_a_b,digits=4)
                    if ( c_a_b < 0){
                        c_a_b <- 0
                    }
                    if ( c_a_b < 0){
                        partialConcordance_al_pr_gj[i,j]=0
                    }
                    if ( c_a_b >= 1 ){
                        partialConcordance_al_pr_gj[i,j]=1
                    } else {
                        partialConcordance_al_pr_gj[i,j]= c_a_b
                    }										
                }				
                
            }
            h=h+1
            globalconcordance[,h]=partialConcordance_al_pr_gj%*%vp
            cat("----------------------------------------------------------","\n")  
            cat("concordance partielle ","ai","_","b",k,sep="","\n")
            cat("----------------------------------------------------------","\n")  								
            print(partialConcordance_al_pr_gj)
            cat("----------------------------------------------------------","\n")
            
            for (i in 1: nrow_pm){
                # Columns' iterations	
                
                for (j in 1:length_cr){
                    if (mmv[j]=="max"){							
                        c_a_b = (pr[k,j]-pm[i,j]+p_t[j])/(p_t[j]-q_t[j])
                    }
                    if (mmv[j]=="min"){
                        c_a_b = (pm[i,j]-pr[k,j]+p_t[j])/(p_t[j]-q_t[j])									
                    }									
                    c_a_b=round(c_a_b,digits=4)
                    if ( c_a_b < 0){
                        c_a_b <- 0
                    }
                    if ( c_a_b < 0){
                        partialConcordance_pr_al_gj[i,j]=0
                    }
                    if ( c_a_b >= 1 ){
                        partialConcordance_pr_al_gj[i,j]=1
                    } else {
                        partialConcordance_pr_al_gj[i,j]= c_a_b
                    }										
                }				
                
            }
            h=h+1
            globalconcordance[,h]=partialConcordance_pr_al_gj%*%vp
            cat("----------------------------------------------------------","\n")  
            cat("concordance partielle ","b",k,"_","ai",sep="","\n")
            cat("----------------------------------------------------------","\n")  								
            print(partialConcordance_pr_al_gj)
            cat("----------------------------------------------------------","\n")
            k=k+1
        } 
        
        cat("----------------------------------------------------------","\n")
        
        # 2. Computation of the partial discordance indices dj(a,b) and dj(b,a).
        
        l=0
        k=1
        while (k<=nrow_pr){
            for (i in 1: nrow_pm){
                # Columns' iterations	
                
                for (j in 1:length_cr){
                    if ((mmv[j]=="max")){
                        d_a_b = (pr[k,j]-pm[i,j]-p_t[j])/(v_t[j]-p_t[j])
                    }
                    if ((mmv[j]=="min")){
                        d_a_b = (pm[i,j]-pr[k,j]-p_t[j])/(v_t[j]-p_t[j])								
                    }
                    
                    d_a_b=round(d_a_b,digits=4)
                    if ( d_a_b < 0){
                        d_a_b <- 0
                    }
                    if ( d_a_b < 0){
                        partialDiscordance_al_pr_gj[i,j]=0
                    }
                    if ( d_a_b >= 1 ){
                        partialDiscordance_al_pr_gj[i,j]=1
                    } else {
                        partialDiscordance_al_pr_gj[i,j]= d_a_b
                    }										
                }				
                
            }
            l=l+1
            for (h in 1: nrow_pm){
                if (max(partialDiscordance_al_pr_gj[h,]==1)){
                    credibility[h,l]=0
                }else if (max(partialDiscordance_al_pr_gj[h,])<globalconcordance[h,l]){
                    credibility[h,l]=globalconcordance[h,l]
                } else if (max(partialDiscordance_al_pr_gj[h,]) > globalconcordance[h,l]){
                    credibility[h,l]=((1-max(partialDiscordance_al_pr_gj[h,]))/(1-globalconcordance[h,l]))*globalconcordance[h,l]
                }
            }
            
            cat("----------------------------------------------------------","\n")  
            cat("partial discordance ","ai","_","b",k,sep="","\n")
            cat("----------------------------------------------------------","\n")  								
            print(partialDiscordance_al_pr_gj)
            cat("----------------------------------------------------------","\n")
            
            for (i in 1: nrow_pm){
                # Columns' iterations	
                
                for (j in 1:length_cr){
                    if ((mmv[j]=="max")){							
                        d_b_a = (pm[i,j]-pr[k,j]-p_t[j])/(v_t[j]-p_t[j])
                    }
                    if ((mmv[j]=="min")){
                        d_b_a = (pr[k,j]-pm[i,j]-p_t[j])/(v_t[j]-p_t[j])
                    }								
                    d_b_a = round(d_b_a,digits=4)
                    if ( d_b_a < 0){
                        d_b_a <- 0
                    }
                    if ( d_b_a < 0){
                        partialDiscordance_pr_al_gj[i,j]=0
                    }
                    if ( d_b_a >= 1 ){
                        partialDiscordance_pr_al_gj[i,j]=1
                    } else {
                        partialDiscordance_pr_al_gj[i,j]= d_b_a
                    }										
                }				
                
            }
            l=l+1
            for (h in 1: nrow_pm){
                if (max(partialDiscordance_pr_al_gj[h,]==1)){
                    credibility[h,l]=0
                }else if (max(partialDiscordance_pr_al_gj[h,])< globalconcordance[h,l]){
                    credibility[h,l]=globalconcordance[h,l]
                } else if (max(partialDiscordance_pr_al_gj[h,]) > globalconcordance[h,l]){
                    credibility[h,l]=((1-max(partialDiscordance_pr_al_gj[h,]))/(1-globalconcordance[h,l]))*globalconcordance[h,l]
                }
            }
            
            cat("----------------------------------------------------------","\n")  
            cat("partial discordance ","b",k,"_","ai",sep="","\n")
            cat("----------------------------------------------------------","\n")  								
            print(partialDiscordance_pr_al_gj)
            cat("----------------------------------------------------------","\n") 
            
            k=k+1
        } 
        
        cat("----------------------------------------------------------","\n")
        cat("----------------  global concordance indexes  ------------","\n")
        print(globalconcordance)
        cat("----------------------------------------------------------","\n")
        
        cat("------------------  Credibility indexes  -----------------","\n")
        print(credibility)
        cat("----------------------------------------------------------","\n")
        
        relations=data.frame(relations)
        colnames(relations)=names_relations
        for (i in 1:nrow_pm){
            # Columns' iterations	    
            for (j in 0:(nrow_pr+1)){
                if (j==0) {
                    relations[i,1]=">"
                } else if (j==nrow_pr+1) {
                    relations[i,nrow_pr+2]="<"
                }else {
                    if (credibility[i,2*j-1]>= lambda){
                        if (credibility[i,2*j]>=lambda){
                            relations[i,j+1]="I"
                        } else {  
                            relations[i,j+1]=">"
                        }
                    } else if (credibility[i,2*j] >=lambda){
                        relations[i,j+1]="<"
                    } else {
                        relations[i,j+1]="R"
                    }
                }
            }																
        }
        cat("------------------  Relations  symbols  -----------------","\n")	
        print(relations)
        cat("----------------------------------------------------------","\n")
        
        
        # pessimistic assignment
        # a) compare alternatives ("a") successively to "b(i)" , for i=p,p-1, ..., 0,
        # b) let "b(h)" = the first profile such that "a" outranks "b(h).",
        # affect "a" to the category C(h+1).
        
        
        # Pessimistic decisions ... 	
        for (i in 1:nrow_pm){
            # Columns' iterations
            ok=0
            for (j in nrow_pr:1){
                if ( (credibility[i,2*j-1]>=0.75) && (credibility[i,2*j]<0.75)) {
                    if (ok==0) {
                        Pessimistic[i,j+1]=1
                        ok=1
                    }
                }
            }
            if (ok==0){
                Pessimistic[i,j]=1
            }														
        }
        cat("------------------  Pessimistic decisions  -----------------","\n")	
        print(Pessimistic)
        cat("----------------------------------------------------------","\n")
        
        # optimistic assignment,
        # a) compare alternatives ("a") successively to "b(i)" ,for  i=1, 2, ..., p+1,
        # b) let "b(h)" = the first profile such that "b(h)" outranks ("a"),
        # affect "a" to the category C(h).
        
        # Optimistic decisions ... 
        for (i in 1:nrow_pm){
            # Columns' iterations
            ok=0
            for (j in 1:nrow_pr){
                if ( (credibility[i,2*j-1]<0.75) && (credibility[i,2*j]>=0.75)) {
                    if (ok==0) {
                        Optimistic[i,j]=1
                        ok=1
                    }
                }
            }
            if (ok==0){
                Optimistic[i,j+1]=1
            }							
        }
        
        cat("------------------  Optimistic decisions  -----------------","\n")	
        print(Optimistic)
        cat("----------------------------------------------------------","\n")
    }
