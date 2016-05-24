
setRefClass("model",
    fields = list(
        model_name = "character",       # e.g. "bernoulli"
        membership_name = "character",  # e.g. "SBM"

        autosave = "character", # autosave filename

        # digests are used here
        digest_already_tried = "list",     # initialization already tried
        digest_already_quality_computed = "list", 
                                           # initialization with quality already
                                           # computed each value is a list with
                                           # two members, first the digest,
                                           # second, the value
        digest_already_splitted = "list",
        digest_already_merged = "list",

        exploration_factor = "numeric",
        explore_max = "numeric",
        explore_min = "numeric",
        memberships = "list",           # found memberships,
        model_parameters = "list",         # found model parameters
        PL = "numeric",                 # Pseudo liklihood of found models
        H = "numeric",                  # Entropy of found models
        ICL = "numeric",                # ICL of found models
        precomputed = "list",
        last_reinitialization_effort = "numeric",
        verbosity = "numeric",
        plotting = "character",
        plotting_data = "list",

        # profiling
        profiling = "numeric",
        profiling_active = "logical",
        profiling_t = "numeric",

        ncores = "numeric"
    ),
    methods = list(
        postinit = function()
        {
            if(membership_name != 'SBM'
               &&
               membership_name != 'SBM_sym'
               &&
               membership_name != 'LBM'
               )
            {
                stop(paste('Membership',membership_name,'unknown. Are you drunk?'))
            }
            if(length(verbosity)==0)
            {
                verbosity<<-6
            }

            if(length(profiling_active)==0)
            {
                profiling_active<<-FALSE
            }

            if(length(autosave)==0)
            {
                autosave<<-""
            }

            if(length(exploration_factor)==0)
            {
                exploration_factor <<- 1.5
            }

            if(length(explore_max)==0)
            {
                explore_max <<- Inf
            }

            if(length(explore_min)==0)
            {
                explore_min <<- 4
            }

            if(length(ncores)==0)
            {
                ncores <<- detectCores()
            }
        },
        save_now = function()
        {
            if(nzchar(autosave)>0)
            {
                saveRDS(.self,file=autosave)
            }
        },
        tic = function()
        {
            if(profiling_active)
            {
                profiling_t<<-cumtime()
            }
        },
        toc = function(field)
        {
            if(profiling_active)
            {
                t2 <- cumtime()
                if(is.na(profiling[field]))
                {
                    profiling[field] <<- 0
                }
                profiling[field] <<- profiling[field] + t2 - profiling_t
                profiling_t <<-t2
            }
        },
        say = function(level,...)
        {
            if(level<=verbosity)
            {
                if(level>1)
                {
                    for(i in 2:level)
                    {
                        cat("    ")
                    }
                }
                cat("-> ")
                cat(paste(...))
                cat("\n")
            }
        },
        estimate = function(reinitialization_effort=1)
        {
            if(!any(last_reinitialization_effort==reinitialization_effort))
            {
                digest_already_splitted <<- list()
                digest_already_merged <<- list()
                last_reinitialization_effort<<-reinitialization_effort
                changing_effort <- TRUE
            }
            else
            {
                changing_effort <- FALSE
            }
            
            if(length(memberships)==0)
            {
                if(membership_name=="LBM")
                {
                    say(1,"Estimation for 2 groups (1+1)")
                    do_with_inits(
                        list(getRefClass(membership_name)(
                            network_size=.self$number_of_nodes())),
                        2,reinitialization_effort)
                }
                else
                {
                    say(1,"Estimation for 1 groups")
                    do_with_inits(
                        list(getRefClass(membership_name)(
                            network_size=.self$number_of_nodes())),
                        1,reinitialization_effort)

                }
            }

            .self$precompute()

            l<-TRUE
            n<-1
            while(l)
            {
                say(1,"Pass",n)

                say(2,"With ascending number of groups")
                ra<-.self$estim_ascend(reinitialization_effort,changing_effort)

                say(2,"With descending number of groups")
                rb<-.self$estim_descend(reinitialization_effort)

                l<-ra||rb
                n<-n+1
                changing_effort<-FALSE
            }
        },

        estim_ascend = function(reinitialization_effort,changing_effort)
        {
            
            if(membership_name=="LBM")
            {
                Q <- 2
            }
            else
            {
                Q <- 1
            }
            
            Q_without_ICL <- max(length(ICL),explore_min)
            Q_stop <- explore_max

            ret<-FALSE
            while(Q<Q_stop && (which.max(ICL)*exploration_factor>length(ICL) || Q<Q_without_ICL))
            {
                Q<-Q+1

                say(3,"For",Q,"groups")
                say(4,'Selecting initialization')

                if(Q>length(ICL) || changing_effort)
                {
                    say(5,"Init from spectral clustering")

                    tic()
                     
                    inits <- .self$provide_init(Q)

                    toc('init_SC')
                    
                }
                else
                {
                    inits <- list()
                }

                say(5,"Init from splitting groups from",Q-1,"groups")
               
                tic() 
                
                inits <- c(inits,.self$split_membership(Q-1))

                toc('init_split')

                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    say(4,"already done")
                    if(Q>length(ICL))
                    {
                        break
                    }
                }
            }
            return(ret)
        },

        estim_descend = function(reinitialization_effort)
        {
            ret<-FALSE
            if(membership_name=="LBM")
            {
                Qmin <- 2
            }
            else
            {
                Qmin <- 1
            }
            if(length(ICL)<=Qmin+1)
            {
                return(FALSE)
            }

            for(Q in seq(length(ICL)-1,Qmin+1))
            {
                say(3,"For",Q,"groups")
                say(4,"Selecting intializations")
                say(5,"Init from merging groups from",Q+1,"groups")

                tic()

                inits <- merge_membership(memberships[[Q+1]])

                toc('init_merges')                

                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    say(4,"Already done")
                }

            }
            return(ret)
        },

        do_with_inits = function(inits,Q,reinitialization_effort)
        {
            say(5,length(inits),"initializations provided")

            tic()
            
            filter<-sapply(
                inits,
                function(init)
                {
                    d<-init$digest()
                    !any(
                        sapply(
                            digest_already_tried,
                            function(x)
                            {
                                x==d
                            }
                        )
                    )
                }
            )

            toc('estimation_already_tried')            

            nb_init_max <- floor(1+4*reinitialization_effort*sqrt(Q))

            say(5,length(inits)-sum(filter),"initializations already used")
            
            if(length(inits)>nb_init_max)
            {
                say(5,'Computing intializations quality')
                quality<-.self$membership_init_quality(inits)
                seuil <- (-sort(-quality))[nb_init_max]
                filter <- filter & (quality >= seuil)
            }

            inits <- inits[filter]

            say(4,"Estimation with",length(inits),"initializations")

            ret<-FALSE

            
            if(length(inits)>0)
            {

                tic()

                results<-parallel_lapply(
                    inits,
                    .self$do_one_estim,
                    mc.cores=ncores,
                    verbose=(verbosity>4))
            
                toc('estimation_run')

                good <- FALSE

                ICLs <- sapply(results, function(r){
                            r$PL - .5*(r$model$n_parameters *
                                        log(.self$data_number())
                                +
                                getRefClass(membership_name)(
                                        from_cc=r$membership)$ICL_penalty())
                      })

                if(length(ICL)>=Q)
                {

                    if(ICL[Q]<max(ICLs))
                    {
                        good <- TRUE
                    }
                }
                else
                {
                    good <- TRUE
                }
               
                toc('estimation_computation_ICL')

                if(membership_name=="SBM" || membership_name=="SBM_sym")
                {
                    xQ <- sapply(
                        results,
                        function(r)
                        {
                            ncol(getRefClass(membership_name)(from_cc=r$membership)$Z)
                        }
                    )
                }
                else
                {
                    xQ <- sapply(
                        results,
                        function(r)
                        {
                            membership<-getRefClass(membership_name)(from_cc=r$membership)
                            return(sqrt(ncol(membership$Z1)*ncol(membership$Z2)))
                        }
                    )
                }
                plotting_data$xQ <<- c(plotting_data$xQ,xQ)
                plotting_data$ICL <<- c(plotting_data$ICL,ICLs)

                toc('updating_plotting_data')

                digest_already_tried <<- c(digest_already_tried,
                                    lapply(inits,function(x){x$digest()}))
                
                toc('estimation_adding_already_tried')

                if(good)
                {
                    
                    say(5,"Better ICL criterion found")
                    say(6,"new ICL:",max(ICLs))
                    say(6,"old ICL:",ICL[Q])
                    

                    kmax<-which.max(ICLs)

                    r<-results[[kmax]]
                    memberships[[Q]] <<-
                        getRefClass(membership_name)(from_cc=r$membership)

                    if(membership_name=="LBM")
                    {
                        say(5,memberships[[Q]]$show_short())
                    }

                    model_parameters[[Q]] <<- r$model
                    PL[Q] <<- r$PL
                    H[Q] <<- r$H
                    ICL[Q] <<- r$PL - .5*(r$model$n_parameters *
                        log(.self$data_number()) + memberships[[Q]]$ICL_penalty())
                    ret<-TRUE

                    toc('estimation_saving_goods')



                
                }
                else
                {
                    say(5,"Useless, no better ICL criterion found")
                    say(6,"better ICL found:",max(ICLs))
                    say(6,"old ICL:",ICL[Q])
                }
                
                tic()

                plot_type<-0
                if(length(plotting)==0)
                {
                    plot_type<-1
                }
                else
                {
                    if(nzchar(plotting))
                    {
                        plot_type<-2
                    }
                }

                if(plot_type>0)
                {
                    if(plot_type==2)
                    {
                        pdf(plotting)
                    }
                    if(membership_name=="LBM")
                    {
                        xlab<-"gemetrical_mean (Q1,Q2)"
                        xQ<-sapply(memberships[2:length(ICL)],function(m){
                                   sqrt(ncol(m$Z1)*ncol(m$Z2))})
                        yICL<-ICL[2:length(ICL)]
                    }
                    else
                    {
                        xlab<-"Q"
                        xQ<-1:length(ICL)
                        yICL<-ICL
                    }
                    par(bty="l")
                    if(max(xQ)<=4)
                    {
                        ylim <- range(yICL)
                    }
                    else
                    {
                        ylim <- range(yICL[xQ>=2])
                    }

                    plot(
                        x=plotting_data$xQ,
                        y=plotting_data$ICL,
                        xlab=xlab,
                        ylab="ICL",
                        ylim=ylim,
                        pch=19
                    )
                    o<-order(xQ)
                    points(x=xQ[o],y=yICL[o],pch=19,col="red")
                    if(plot_type==2)
                    {
                        dev.off()
                    }
                }

            }

            .self$save_now()

            return(ret)
        },

        membership_init_quality = function(inits)
        {
            tic()

            quals <- sapply(
                inits,
                function(init)
                {
                    qual <- digest_already_quality_computed[[init$digest()]]
                    if(is.null(qual))
                    {
                        return(NA)
                    }
                    else
                    {
                        return(qual)
                    }
                }
            )
           
            toc('quality_already_computed') 
            
            
            if(any(is.na(quals)))
            {
                inits<-inits[is.na(quals)]

                naquals<- simplify2array(
                    parallel_lapply(
                        inits,
                        function(membership_init)
                        {
                            

                            r <- dispatcher(membership_name,
                                            membership_init$to_cc(),
                                            model_name,
                                            .self$network_to_cc(),
                                            FALSE)

                            local_ICL <- r$PL - .5*(r$model$n_parameters *
                                                log(.self$data_number())
                                        +
                                        getRefClass(membership_name)(
                                                from_cc=r$membership)$ICL_penalty())
                        },
                        mc.cores=ncores,
                        verbose=(verbosity>4)
                    )
                )

                for(i in 1:length(inits))
                {
                    digest_already_quality_computed[[inits[[i]]$digest()]] <<- naquals[i]
                }

                quals[is.na(quals)] <- naquals
           
                toc('quality_computation') 
            }
            return(quals)
        },

        do_one_estim = function(membership_init)
        {
            return(
                dispatcher(membership_name,
                    membership_init$to_cc(),
                    model_name,
                    .self$network_to_cc(),
                    TRUE)
                )
        },
        split_membership = function(Q)
        {
            d<-memberships[[Q]]$digest()
            if(any(sapply(digest_already_splitted,function(x){x==d})))
            {
                return(list())
            }
            else
            {
                splitted_membership<-.self$split_membership_model(Q)
                digest_already_splitted <<- c(digest_already_splitted,list(d))
                return(splitted_membership)
            }
        },
        merge_membership = function(membership)
        {
            d<-membership$digest()
            if(any(sapply(digest_already_merged,function(x){x==d})))
            {
                return(list())
            }
            else
            {
                merged_membership<-membership$merges()
                digest_already_merged <<- c(digest_already_merged,list(d))
                return(merged_membership)
            }
        },
        precompute = function() {},
        plot_obs_pred = function(Q) {},
        plot_parameters = function(Q) {},
        show = function()
        {
            cat("blockmodels object\n")
            cat(paste("    model:",model_name,"\n"))
            cat(paste("    membership:",membership_name,"\n"))
            cat(paste("    network:",.self$show_network(),"\n"))
            if(length(ICL)>0)
            {
                cat(paste("    maximum of ICL:",memberships[[which.max(ICL)]]$show_short(),"\n"))
                cat("    Most usefull fields and methods:\n")
                cat("        The following fields are indexed by the number of groups:\n")
                cat("            $ICL : vector of ICL\n")
                cat("            $PL : vector of pseudo log liklihood\n")
                cat("            $memberships : list of memberships founds by estimation\n")
                cat("                           each membership is represented object\n")
                cat("            $model_parameters : models parameters founds by estimation\n")
                cat("        Estimation methods:\n")
                cat("            $estimate(reinitalization_effort=1) : to run again estimation with a\n")
                cat("                                                  higher reinitalization effort\n")
                cat("        Plotting methods:\n")
                cat("            $plot_obs_pred(Q) : to plot the obeserved and predicted network for Q groups\n")
                cat("            $plot_parameters(Q) : to plot the model_parameters for Q groups\n")
                cat("            Please note that each membership object have a plotting pethod\n")
            }
            else
            {
                cat("    Estimation not done.\n")
                cat("    Run $estimate(). You can specify a reinitialization effort, by default 1.\n")
            }

        }
    )
)

