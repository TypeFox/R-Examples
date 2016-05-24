#in development code
#[TBC - NUMBER] functions 

#cond handlers

#do properly



###########################
###########################
#
###########################
###########################






groupsAndZcasesPanelHandler <- function(panel=NULL, ..., plot = TRUE, process = TRUE){

    #wrapper for groupsPanelHandler and zcasesPanelHandler

    temp <- names(list(...))

    if(!"groups" %in% temp & !"zcases" %in% temp)
        panel(...) else if(!"groups" %in% temp)
            zcasesPanelHandler(..., panel = panel, plot = plot, process = process) else if(!"zcases" %in% temp)
                groupsPanelHandler(..., panel = panel, plot = plot, process = process) else {

#this is currently not fully working
#here we know we have both groups and zcases

                    #reset args
                    extra.args <- list(...)
                    extra.args$cond.args <- unique(c(extra.args$group.args, 
                                                     extra.args$zcase.args))
                    extra.args$conds <- paste(extra.args$groups,                    
                                              extra.args$zcases,
                                              sep = "#by#")

##could rationalise the next two bits?
##

                    if("group.args" %in% names(extra.args))
                        for(i in extra.args$group.args)
                            extra.args[[i]] <- rep(extra.args[[i]], length(extra.args$zcase.ids))

                    if("zcase.args" %in% names(extra.args))
                        for(i in extra.args$zcase.args)
                            extra.args[[i]] <- rep(extra.args[[i]], each = length(extra.args$zcase.ids))
                    
#chop out groups, etc
extra.args <- extra.args[!names(extra.args) %in% c("groups", "group.ids", "zcases", "zcase.ids")]

                    do.call(condsPanelHandler, listUpdate(extra.args, 
                                                          list(panel=panel, plot=plot, process=process, 
                                                               reset.conds = c("groups","zcases"))))
                }
}










#####################################




zcasesPanelHandler <- function(..., zcases = NULL, panel = NULL, by.zcase = NULL,
         process = TRUE, plot = TRUE){

     #restructure 
     #so we can use groupsHandler code

     extra.args <- list(...)

     if(!"cond.args" %in% names(extra.args))
         extra.args$cond.args <- extra.args$zcase.args

     if(is.null(extra.args$panel.elements)){
                          temp <- lapply(extra.args, function(x) length(x) == length(extra.args$x))
                          extra.args$panel.elements <- names(extra.args)[temp]
                      }

     extra.args$cond.ids <- if(is.null(extra.args$zcase.ids))
                                if(is.factor(zcases)) levels(zcases) else unique(zcases) else 
                                    extra.args$zcase.ids
     
     extra.args <- listUpdate(extra.args,
                              list(panel = panel, conds = zcases,
                                   by.cond = by.zcase, reset.conds = "zcases", 
                                   process = process, plot = plot))

     #run
     do.call(condsPanelHandler, extra.args)
}








###########################################

groupsPanelHandler <- function(..., groups = NULL, panel = NULL, by.group = NULL,
         process = TRUE, plot = TRUE){

     #restructure 
     #so we can use groupsHandler code

     extra.args <- list(...)

     if(!"cond.args" %in% names(extra.args))
         extra.args$cond.args <- extra.args$group.args

     if(is.null(extra.args$panel.elements)){
                          temp <- lapply(extra.args, function(x) length(x) == length(extra.args$x))
                          extra.args$panel.elements <- names(extra.args)[temp]
                      }

     extra.args$cond.ids <- if(is.null(extra.args$group.ids))
                                if(is.factor(groups)) levels(groups) else unique(groups) else 
                                    extra.args$group.ids
     
     extra.args <- listUpdate(extra.args,
                              list(panel = panel, conds = groups,
                                   by.cond = by.group, reset.conds = "groups", 
                                   process = process, plot = plot))

     #run
     do.call(condsPanelHandler, extra.args)
}









#############################################



condsPanelHandler <- function(..., conds = NULL, panel = NULL, by.cond = NULL,
         process = TRUE, plot = TRUE){


##think about this
##process and plot need to be passed on




    #extra.args <- list(...)

    extra.args <- listUpdate(list(...), list(process=process, plot=plot))



    if(is.null(panel) && is.null(by.cond)) return(NULL)
    if(!process && !plot) return(NULL)

    if(is.null(conds)) return(do.call(panel, extra.args))

    if(is.null(by.cond))
        by.cond <- list(panel)

    cond.ids <- if(is.null(extra.args$cond.ids))
                     if(is.factor(conds)) levels(conds) else unique(conds) else 
                         extra.args$cond.ids


#might want to rename this
#might be fine
    panel.elements <- if(is.null(extra.args$panel.elements)){

                          #if not known 
                          #assume we want all that are the same length as x
                          #could be bad for col.regions
                          temp <- lapply(extra.args, function(x) length(x) == length(extra.args$x))
                          names(extra.args)[temp]

                      } else extra.args$panel.elements


  
     n <- rep(1:length(by.cond), length.out = length(cond.ids))


##could the next bit be reduced/simplify

     if(process & !plot){

         #the process method 
         out <- lapply(1:length(n), function(i){

             #process and plot
             #or plot only
             #or null case

             temp <- extra.args[names(extra.args) %in% panel.elements]
             temp <- lapply(temp, function(x) x[conds == cond.ids[i]])

             temp <- listUpdate(extra.args, temp)
             for(j in temp$cond.args)
                 temp[[j]] <- temp[[j]][i]

#            temp[names(temp) %in% c(temp$cond.args] <- temp[names(temp) %in% c(temp$cond.args][i]

             out2 <- if(is.list(by.cond[[n[i]]])) do.call(panel, listUpdate(temp, by.cond[[n[i]]])) else 
                        if(is.function(by.cond[[n[i]]])) do.call(by.cond[[n[i]]], temp)
             ref <- if(is.null(out2$x)) out2[[1]] else out2$x

             out2$conds <- cond.ids[i]

#this is the bit that gets renamed 
#in each of the functions
if("reset.conds" %in% names(extra.args)){



    if(length(extra.args$reset.conds) > 1) {
        for(k in 1:length(extra.args$reset.conds)){
            out2[[extra.args$reset.conds[k]]] <- unlist(lapply(strsplit(out2$conds, "#by#"), function(x) x[k]))
        }
    } else {
        names(out2)[names(out2) == "conds"] <- extra.args$reset.conds
    }
    out2 <- out2[names(out2)[!names(out2) %in% "conds"]]
}




             temp <- temp[names(temp) %in% panel.elements]
             temp <- temp[!names(temp) %in% names(out2)]
             temp <- lapply(temp, function(x) rep(x, length.out=length(ref)))

             listUpdate(out2, temp)

         })

####merge the lists         

         out <- as.list(do.call(rbind, lapply(out, function(x) as.data.frame(x, stringsAsFactors=FALSE))))
         
         

     } else {  


         for(i in 1:length(n)){


#this needs to be better
#terminates badly

             #process and plot
             #or plot only
             #or null case


             temp <- extra.args[names(extra.args) %in% panel.elements]
             temp <- lapply(temp, function(x) x[conds == cond.ids[i]])

             temp <- listUpdate(extra.args, temp)

             for(j in temp$cond.args)
                 temp[[j]] <- temp[[j]][i]
             temp$subscripts <- 1:length(temp$x)


#            temp[names(temp) %in% c(temp$cond.args] <- temp[names(temp) %in% c(temp$cond.args][i]

             if(is.list(by.cond[[n[i]]])) do.call(panel, listUpdate(temp, by.cond[[n[i]]]))
             if(is.function(by.cond[[n[i]]])) do.call(by.cond[[n[i]]], temp)
         }
     }    

}







###########################
###########################
#
###########################
###########################

groupsHandler <- function(z = NULL, groups = NULL, ..., group.ids = NULL, 
                          handler = "zzz"){

    #setup
    extra.args <- list(...)

    #do nothing if groups not set
    if(is.null(groups)) return(z)
    
    #if groups and group.ids set
    #group.ids should be robust
    group.ids <- unique(c(group.ids, 
                          if(is.factor(groups)) levels(groups) else groups))
    
    #assign z
    if(is.null(z)){

        #if z is absent set groups as default
        z <- groups

#this could be restructured

        #we also assume this needs handling
        known.handlers <- c("col", "cex", "pch")
        if(handler %in% known.handlers){

              temp <- as.numeric(factor(group.ids))
              extra.args$z <- temp
              extra.args$ref <- NULL
              extra.args <- extra.args[names(extra.args) != "zlim"]
              z <- do.call(paste(handler[1], "Handler", sep=""), extra.args)
              
        } else {
              
              z <- as.numeric(factor(group.ids))

        }


    } else {

        #if z is present set groups robustly
        z <- zHandler(z = z, ref = group.ids)
        #transpose onto groups
##here        ##z <- z[as.numeric(factor(groups, group.ids))]

        #we assume this is fine for 'as.is' use
        #because this is what the user set arg as

    }

    #output
    z


}









############################################

#fix
################################
#to be checked should this
#be zcases.ids or zcase.ids
################################


zcasesHandler <- function(z = NULL, zcases = NULL, ..., zcases.ids = NULL, 
                          handler = "zzz"){

     #restructure 
     #so we can use groupsHandler code
     extra.args <- listUpdate(list(...),
                              list(z = z, handler=handler,
                                   groups = zcases, 
                                   group.ids = zcases.ids))
     do.call(groupsHandler, extra.args)

}









#############################################


groupsAndZcasesHandler <- function(..., loa.settings = NULL){

    #set up
    extra.args <- list(...)
    if(is.null(loa.settings)) 
        loa.settings <- list()
    
#tidy this
#groups and zcases to formals
#groups and zcase overlap test to top

    #check for groups
    if("groups" %in% names(extra.args)){
        if(is.null(extra.args$group.ids))
            extra.args$group.ids <- if(is.factor(extra.args$groups))
                                                   levels(extra.args$groups) else unique(extra.args$groups)

        if(is.null(extra.args$group.args))
            extra.args$group.args <- loa.settings$group.args
        extra.args$ignore <- unique(c(extra.args$ignore, extra.args$group.args))

        for(i in extra.args$group.args)
            extra.args[[i]] <- do.call(groupsHandler, 
                                                  listUpdate(extra.args, 
                                                             list(z = extra.args[[i]], handler=i)))
    }

    #check for zcases
    if("zcases" %in% names(extra.args)){
        if(is.null(extra.args$zcase.ids))
            extra.args$zcase.ids <- if(is.factor(extra.args$zcases))
                                                   levels(extra.args$zcases) else unique(extra.args$zcases)
        if(is.null(extra.args$zcase.args))
            extra.args$zcase.args <- loa.settings$zcase.args

#to think through
#if zcases and groups are set 
#common args need to be compared
        if("groups" %in% names(extra.args)){
 

#print(
#any(extra.args$zcase.args %in% extra.args$group.args)
#)

            if(any(extra.args$zcase.args %in% extra.args$group.args)){
                temp <- extra.args$zcase.args[!extra.args$zcase.args %in% extra.args$group.args]
                if(length(temp)<1){
                    extra.args$zcase.args <- NULL
                    warning("zcase.args/group.args conflict\n", 
                            "no requested args not already applied by groups", call.=FALSE)
                } else {
                    extra.args$zcase.args <- temp
                    warning("zcase.args/group.args conflict\n", 
                            "one or more requested args already applied by groups", call.=FALSE)
                }
           }  
        }
        
        extra.args$ignore <- unique(c(extra.args$ignore, extra.args$zcase.args))

        for(i in extra.args$zcase.args)
            extra.args[[i]] <- do.call(zcasesHandler, 
                                                  listUpdate(extra.args, 
                                                             list(z = extra.args[[i]], handler=i)))
    }

    extra.args

}







###########################################################



stepwiseZcasesGlyphHandler <- function(zcases = NULL, ..., zcase.ids = NULL, 
                              panel.elements = NULL, loaGlyph = NULL){

    #########################
    #stepwise Glyph handler 
    #where glyph elements are 
    #by zcases
    #########################

    ######################
    #set up
    ######################

    extra.args <- list(...)


    #########################
    #checks
    #########################

    if(is.null(zcases)){

        #stop?

    }
    
    if(is.null(loaGlyph)){

        #stop?

    }

    if(is.null(panel.elements)){

        #make a substitute panel.elements

    }

    if(is.null(zcase.ids)){

       #make substitute zcase.ids
       zcase.ids <- if(is.factor(zcases)) 
                        levels(zcases) else unique(zcases)

    }


    #######################
    #main variables
    #######################

    zcase.count <- length(zcase.ids)
    zcase.len <- max(sapply(zcase.ids, function(x) length(extra.args$x[zcases==x])))

    #might want to think about warning if different x cases are different
    #might want to drop x from test?
    #this would simplify if we assumed all zcases were same length
    #which I think we can?

    for(a in 1:zcase.len){
    
        for(b in 1:zcase.count){

            temp <- extra.args

            temp$zrow.ids <- temp$z[seq(a, length(extra.args$z), zcase.len)]

            for(i in panel.elements){

                temp[[i]] <- temp[[i]][zcases == zcase.ids[b]][a]

            }

            if("zcase.args" %in% names(extra.args)){

                for(i in extra.args$zcase.args){

                    temp[[i]] <- extra.args[[i]][zcase.ids==zcase.ids[b]]

                }

            }

    #could move a lot of the following out of formals?
 
            temp$zcase.ref <- c(b, zcase.count)
            temp$zcase.ids <- zcase.ids
            temp$count.ref <- c(a, zcase.len)

            do.call(loaGlyph, temp)

        }
    
    }

}







