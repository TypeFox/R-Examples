
setRefClass("scalar_model",
    contains = "model",
    fields = list(
        adj = "matrix"
    ),
    methods = list(
        postinit = function()
        {
            callSuper()
            if(membership_name=="SBM" || membership_name=="SBM_sym")
            {
                if(nrow(adj)!=ncol(adj))
                {
                    stop(paste("The adjacency matrix does not have the same number of rows and columns.",
                               "Vocatus periculosum ad sanitatem est."))
                }
            }

            if(membership_name=="SBM_sym")
            {
                if(isSymmetric(adj))
                {
                    adj <<- (adj+t(adj))/2
                }
                else
                {
                    stop("Adjacency matrix is not symmetric. You need more sleep.")
                }
            }
        },
        number_of_nodes = function() { dim(adj) },
        show_network = function()
        {
            paste(nrow(adj),"x",ncol(adj),"scalar network")
        },  
        network_to_cc = function() { list(adjacency = adj) },
        split_membership_model = function(Q)
        {
            membership <- memberships[[Q]]
            error <- .self$residual(Q)

            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                result <- list()
                for(q in 1:Q)
                {
                    sub_classif <- coordinates_split(
                        cbind(error, t(error)),
                        membership$Z[,q]
                        )
                    Z <- cbind(membership$Z,membership$Z[,q])
                    Z[,q] <- Z[,q]*sub_classif
                    Z[,Q+1] <- Z[,Q+1]*(1-sub_classif)
                    result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z=Z))
                        ))
                }
                return(result)
            }
            if(membership_name == "LBM")
            {

                Q1<-dim(membership$Z1)[2]
                Q2<-dim(membership$Z2)[2]

                result <- list()
                for(q in 1:Q1)
                {
                    sub_classif <- coordinates_split(
                        cbind(error),
                        membership$Z1[,q]
                        )
                    Z1 <- cbind(membership$Z1,membership$Z1[,q])
                    Z1[,q] <- Z1[,q]*sub_classif
                    Z1[,Q1+1] <- Z1[,Q1+1]*(1-sub_classif)
                    result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z1=Z1,Z2=membership$Z2))
                        ))
                }
                
                for(q in 1:Q2)
                {
                    sub_classif <- coordinates_split(
                        cbind(t(error)),
                        membership$Z2[,q]
                        )
                    Z2 <- cbind(membership$Z2,membership$Z2[,q])
                    Z2[,q] <- Z2[,q]*sub_classif
                    Z2[,Q2+1] <- Z2[,Q2+1]*(1-sub_classif)
                    result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z1=membership$Z1,Z2=Z2))
                        ))
                }
                return(result)
            }
        },
        data_number = function()
        {
            if(membership_name=="SBM")
            {
                return(dim(adj)[1]*(dim(adj)[1]-1))
            }
            if(membership_name=="SBM_sym")
            {
                return(dim(adj)[1]*(dim(adj)[1]-1)/2)
            }
            else
            {
                return(dim(adj)[1]*(dim(adj)[2]))
            }
        },
        precompute = function()
        {
            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                if(length(precomputed)>0)
                {
                    return()
                }
                else
                {
                    if(length(ICL)!=0)
                    {
                        say(1,"Computation of eigen decomposition used for initalizations")
                        
                        error <- .self$residual(1)
                        W<- error %*% t(error)
                        W<-1/(1+exp(-W/sd(W)))
                        D<- diag(1/sqrt(rowSums(W)))
                        L<- D %*% W %*% D

                        precomputed$eigen <<- eigen(L, symmetric=TRUE)
                        
                        cat("\n")
                    }
                }
            }
            if(membership_name == "LBM")
            {
                if(length(precomputed)>0)
                {
                    return()
                }
                else
                {
                    if(length(ICL)!=0)
                    {
                        say(1,"Computation of eigen decomposition used for initalizations")
                        error <- .self$residual(2)
                        
                        say(2,"for rows")
                        
                        W1<- error %*% t(error)
                        W1<-1/(1+exp(-W1/sd(W1)))
                        D1<- diag(1/sqrt(rowSums(W1)))
                        L1<- D1 %*% W1 %*% D1

                        precomputed$eigen1 <<- eigen(L1, symmetric=TRUE)
                        
                        say(2,"for cols")
                        
                        W2<- t(error) %*% error
                        W2<-1/(1+exp(-W2/sd(W2)))
                        D2<- diag(1/sqrt(rowSums(W2)))
                        L2<- D2 %*% W2 %*% D2

                        precomputed$eigen2 <<- eigen(L2, symmetric=TRUE)
                        
                        cat("\n")
                    }
                }
            }
        },
        provide_init = function(Q)
        {
            .self$precompute()
            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                return(
                    list(
                        getRefClass(membership_name)(
                            classif=blockmodelskmeans(
                                as.matrix(precomputed$eigen$vectors[,1:Q]),
                                Q
                            )
                        )
                    )
                )
            }
            if(membership_name == "LBM")
            {
                result <- list()
                for(Q1 in 1:(Q-1))
                {
                    Q2<-Q-Q1
                    if(Q1<=nrow(adj) && Q2<=ncol(adj))
                    {
                        result[[Q1]] <- getRefClass(membership_name)(
                            classif=list(
                                blockmodelskmeans(as.matrix(precomputed$eigen1$vectors[,1:Q1]),Q1),
                                blockmodelskmeans(as.matrix(precomputed$eigen2$vectors[,1:Q2]),Q2)
                            )
                        )
                    }
                    
                }
                result <- result[!sapply(result,is.null)]

                return(result)
            }
        },
        plot_obs_pred = function(Q)
        {
            pred <- .self$prediction(Q)

            par(mfrow=c(2,1))
            if(membership_name == "LBM")
            {
                order1 <- order(memberships[[Q]]$map()$C1)
                order2 <- order(memberships[[Q]]$map()$C2)
            }
            else
            {
                order1 <- order(memberships[[Q]]$map()$C)
                order2 <- order1
            }

            rn<-rownames(adj)
            cn<-colnames(adj)

            if(is.null(rn))
            {
                rn<-1:nrow(adj)
            }
            if(is.null(cn))
            {
                cn<-1:ncol(adj)
            }

            matrixplot(.self$plot_transform(adj[order1,order2]),rowlabels=rn[order1],collabels=cn[order2])
            
            matrixplot(.self$plot_transform(pred[order1,order2]),rowlabels=rn[order1],collabels=cn[order2])
        },
        plot_transform = function(x){x},
        residual = function(Q)
        {
            adj-.self$prediction(Q)
        }

            
    )
)


