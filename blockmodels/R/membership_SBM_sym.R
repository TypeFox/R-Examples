
# membership definition

setRefClass("SBM_sym",
    fields = list(
        Z="matrix"
    ),
    methods = list(
        initialize = function(network_size=FALSE,classif=FALSE,from_cc=FALSE)
        {
            if(!classif[1])
            {
                if(network_size[1])
                {
                    Z <<- matrix(1,nrow=network_size[1],ncol=1)
                }
                else
                {
                    Z <<- from_cc$Z
                }
            }
            else
            {
                fclassif <- factor(classif)
                classif <- as.numeric(fclassif)
                Q <- length(levels(fclassif))
                Z <<- matrix(0,nrow=length(classif),ncol=Q)
                for(i in 1:length(classif))
                {
                    Z[i,classif[i]] <<- 1
                }
            }
        },
        digest = function()
        {
            digest::digest(order_round_matrix(Z),algo='sha256')
        },
        show_short = function()
        {
            paste(ncol(Z),"groups")
        },
        show = function()
        {
            cat("SBM_sym membership\n")
            cat("    Groups:",paste(ncol(Z),"groups\n"))
            cat("    Nodes:",paste(nrow(Z),"nodes\n"))
            cat("    Usefull fields and methods:\n")
            cat("        $Z : matrix of nodes memberships\n")
            cat("        $plot() : plot the memberships\n")
        },
        to_cc = function()
        {
            list(Z=Z)
        },
        map = function()
        {
            list(
                C=apply(Z,1,which.max)
            )
        },
        ICL_penalty = function()
        {
            (dim(Z)[2]-1)*log(dim(Z)[1])
        },
        merges = function()
        {
            result <- list()
            Q <- dim(Z)[2]
            for(k1 in 1:(Q-1))
            {
                for(k2 in (k1+1):Q)
                {
                    Z2<-Z[,-k2]
                    Z2[,k1]<-Z[,k1]+Z[,k2]
                    result <- c(result,list(getRefClass('SBM_sym')(from_cc=list(Z=Z2))))
                }
            }
            return(result)
        },
        plot = function()
        {
            rn<-rownames(Z)
            if(is.null(rn))
            {
                rn<-1:nrow(Z)
            }
            ordering <- order(.self$map()$C)
            matrixplot(as.matrix(Z[ordering,]),rowlabels=rn[ordering])
        }


    )
)
                

