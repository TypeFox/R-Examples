landscape.sample.stages <- function(Rland,ns=NULL,svec=NULL)
    {
        if (is.null(svec))
            {
                Rland
            } else { #there are stages to check out
                stgs <- 0:((Rland$intparam$stages * Rland$intparam$habitats)-1)
                if (sum(!(svec%in%stgs))>0) stop ("you have specified demographic stages that do no occur in this landscape")
                Rland$individuals <- Rland$individuals[Rland$individuals[,1]%in%svec,]
                if (!is.null(ns))
                    {
                        for (s in svec)
                            {
                                inds <- Rland$individuals[Rland$individuals[,1]==s,]
                                Rland$individuals <-  Rland$individuals[Rland$individuals[,1]!=s,]
                                inddim <- dim(inds)[1]
                                if (inddim>0)
                                    if (inddim<ns)
                                        Rland$individuals <- rbind(Rland$individuals,inds)
                                    else
                                        Rland$individuals <- rbind(Rland$individuals,inds[sample(1:inddim,ns,replace=F),])
                            }
                        Rland$individuals <- Rland$individuals[order(Rland$individuals[,1],Rland$individuals[,4]),]
                    }
                Rland
            }
    }

landscape.sample.pops <- function(Rland,ns=NULL,pvec=NULL)
    {
        if (is.null(pvec))
            {
                Rland
            } else { #there are pops
                pops <- 1:Rland$intparam$habitats
                if (sum(!(pvec%in%pops))>0) stop ("you have specified populations that can not occur in this landscape")
                Rland$individuals <- Rland$individuals[landscape.populations(Rland)%in%pvec,]
                if (!is.null(ns))
                    {
                        for (p in pvec)
                            {
                                inds <- Rland$individuals[landscape.populations(Rland)==p,]
                                Rland$individuals <-  Rland$individuals[landscape.populations(Rland)!=p,]
                                inddim <- dim(inds)[1]
                                if (!is.null(inddim))
                                    {
                                        if (inddim<ns)
                                            Rland$individuals <- rbind(Rland$individuals,inds)
                                        else
                                            Rland$individuals <- rbind(Rland$individuals,inds[sample(1:inddim,ns,replace=F),])
                                    }
                            
                            }
                        Rland$individuals <- Rland$individuals[order(Rland$individuals[,1],Rland$individuals[,4]),]
                    }
                Rland
            }
    }



landscape.sample <- function (Rland, np = NULL, ns = NULL, pvec = NULL, svec=NULL) 
{
    if ((!is.null(svec))&(!is.null(pvec))) stop ("either pvec or svec should be specified, not both")
    if ((!is.null(np))|(!is.null(pvec)))
        {
            if (is.null(pvec))
                {
                    if (Rland$intparam$habitats<np){stop("can't sample more populations than exist in landscape")}
                    pvec <- sample(1:Rland$intparam$habitats,np,F)
                } else {
                        pvec <- unique(pvec)
                    }
            if ((!is.null(np))&(!is.null(pvec)))
                if (np!=length(pvec)) stop ("if you specify both np and pvec, the length of pvec has to equal np")
            retland <- landscape.sample.pops(Rland,ns=ns,pvec=pvec)
        } else 
            if (!is.null(svec))
                {
                    retland <- landscape.sample.stages(Rland,ns=ns,svec=svec)
                } else
                    if (!is.null(ns))
                        {
                            retland <- landscape.sample.pops(Rland,ns=ns,pvec=1:Rland$intparam$habitats)
                        } else retland <- Rland
    retland
}
