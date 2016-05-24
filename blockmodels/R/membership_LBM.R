
# membership definition

setRefClass("LBM",
    fields = list(
        Z1="matrix",
        Z2="matrix"
    ),
    methods = list(
        initialize = function(network_size=NULL,classif=NULL,from_cc=NULL)
        {
            if(length(classif)==0)
            {
                if(length(network_size)>0)
                {
                    Z1 <<- matrix(1,nrow=network_size[1],ncol=1)
                    Z2 <<- matrix(1,nrow=network_size[2],ncol=1)
                }
                else
                {
                    Z1 <<- from_cc$Z1
                    Z2 <<- from_cc$Z2
                }
            }
            else
            {
                fclassif1 <- factor(classif[[1]])
                classif1 <- as.numeric(fclassif1)
                
                fclassif2 <- factor(classif[[2]])
                classif2 <- as.numeric(fclassif2)
                
                Q1 <- length(levels(fclassif1))
                Q2 <- length(levels(fclassif2))

                Z1 <<- matrix(0,nrow=length(classif1),ncol=Q1)
                Z2 <<- matrix(0,nrow=length(classif2),ncol=Q2)

                for(i in 1:length(classif1))
                {
                    Z1[i,classif1[i]] <<- 1
                }
                
                for(i in 1:length(classif2))
                {
                    Z2[i,classif2[i]] <<- 1
                }
            }
        },
        digest = function()
        {
            digest::digest(
                list(
                    order_round_matrix(Z1),
                    order_round_matrix(Z2)
                ),algo='sha256')
        },
        show_short = function()
        {
            paste(ncol(Z1),"row groups,",ncol(Z2),"col groups")
        },
        show = function()
        {
            cat("LBM membership\n")
            cat(paste("    Groups:",ncol(Z1),"row groups,",ncol(Z2),"col groups\n"))
            cat(paste("    Nodes:",nrow(Z1),"row nodes,",nrow(Z2),"col nodes\n"))
            cat("    Usefull fields and methods:\n")
            cat("        $Z1 : matrix of row nodes memberships\n")
            cat("        $Z2 : matrix of col nodes memberships\n")
            cat("        $plot() : plot the memberships\n")
        },
        to_cc = function()
        {
            list(Z1=Z1,Z2=Z2)
        },
        map = function()
        {
            list(
                C1 = apply(Z1,1,which.max),
                C2 = apply(Z2,1,which.max)
            )
        },
        ICL_penalty = function()
        {
            (dim(Z1)[2]-1)*log(dim(Z1)[1]) + (dim(Z2)[2]-1)*log(dim(Z2)[1])
        },
        merges = function()
        {
            result <- list()
            Q1 <- dim(Z1)[2]
            if(Q1>1)
            {
                for(k1 in 1:(Q1-1))
                {
                    for(k2 in (k1+1):Q1)
                    {
                        Zn<-as.matrix(Z1[,-k2])
                        Zn[,k1]<-Z1[,k1]+Z1[,k2]
                        result <- c(result,list(getRefClass('LBM')(from_cc=list(Z1=Zn,Z2=Z2))))
                    }
                }
            }
            
            Q2 <- dim(Z2)[2]
            if(Q2>1)
            {
                for(k1 in 1:(Q2-1))
                {
                    for(k2 in (k1+1):Q2)
                    {
                        Zn<-as.matrix(Z2[,-k2])
                        Zn[,k1]<-Z2[,k1]+Z2[,k2]
                        result <- c(result,list(getRefClass('LBM')(from_cc=list(Z1=Z1,Z2=Zn))))
                    }
                }
            }

            return(result)
        },
        plot = function()
        {
            par(mfrow=c(1,2))

            rn1<-rownames(Z1)
            if(is.null(rn1))
            {
                rn1<-1:nrow(Z1)
            }
            ordering1 <- order(.self$map()$C1)
            matrixplot(as.matrix(Z1[ordering1,]),rowlabels=rn1[ordering1])

            rn2<-rownames(Z2)
            if(is.null(rn2))
            {
                rn2<-1:nrow(Z2)
            }
            ordering2 <- order(.self$map()$C2)
            matrixplot(as.matrix(Z2[ordering2,]),rowlabels=rn2[ordering2])
        }   

    )
)
                

