#in development code
#[TBC - NUMBER] functions 

#panelPal.old
#panelPal
#loaHandler 


#NOTE: much borrowed from lattice 



#panelPal needs a do not subscript...


##############################
##############################
##panelPal.old
##############################
##############################


#this needs a lot of tidying

panelPal.old <- function(x, y, subscripts, at, col.regions, ..., 
                    panel = panel.xyplot, ignore = NULL,
                    group.fun = NULL){

           #############
           #got to be better way of doing this
           #make all panel-specific


###########
#print("HI")
##########
           extra.args <- list(...)
           temp <- if(is.null(ignore))
                       extra.args else extra.args[!names(extra.args) %in% ignore]

           if(!is.null(subscripts)){
               temp <- lapply(temp, function(x)
                           x <- if(length(x)>1) x[subscripts] else x )
                               subscripts <- 1:length(subscripts)
           }
           temp <- listUpdate(list(x = x, y = y, z = temp$z, at = at, 
                              col.regions=col.regions, subscripts=subscripts),
                              temp)
           
           if(!is.null(ignore))
               temp <- listUpdate(temp, extra.args[names(extra.args) %in% ignore])
   
###########
#print(names(temp))
###########

           if(!"groups" %in% names(temp)) return(do.call(panel, temp))

#note this action 

          groups <- temp$groups
          temp <- temp[names(temp) != "groups"]
 


           grp <- if(is.factor(groups)) levels(groups) else unique(groups)
  ##         temp.fun <- function(temp, i){
  ##                         temp2 <- lapply(temp, function(x) x <- if (length(x) > 1) 
  ##                                    x[groups==i] else x)
  ##                        
  ##                    }



#does the next bit need nicer coloring?
#does the next bit an option from above
           
          if(is.null(group.fun)){
              group.fun <- list()
              for(i in grp){
                  group.fun[[i]] <- list(col= i, pch = i)
              }}

#print(group.fun[1])


           for(i in grp){

##add catcher for missing 
##group.fun
##update col and symbol 
## if mono?

                temp2 <- lapply(temp, function(x) x <- if (length(x) > 1) 
                                x[groups==i] else x)


###remove later
#           print(temp2$x)
####


#does the next bit need error catching
#only do else if function?

#does the next bit need 

           if(is.list(group.fun[[i]]))
                do.call(panel, listUpdate(temp2, group.fun[[i]])) else 
                do.call(group.fun[[i]], temp2)
           }


}


##############################
##############################
##panelPal
##############################
##############################


panelPal <- function(ans, panel = NULL, preprocess = FALSE,
         reset.xylims = FALSE, legend = NULL,
         by.group = NULL, by.zcase = NULL, ...){

#panelPal v2
#kr

#changes 
#######################
#from formals
#removed safe.mode
#changed group.fun to by.group 
#    or group.apply 
#added ... 
#######################
#safe.mode need changing 
#in loaPlot panelPal
#call
#######################



#to do
#######################
#other caption terms 
#need to be ignored
#when extrapolating
#######################
#group handling
#######################
#

#to think about
#######################
#drop ref from common after transfer?
#######################
#col/cex handling when either are 
#safe.mode managed  
#######################
#

    #move all required and expanded information to panel.args[[1]]
    #drop ignores
    #make expanded list

    #check for for set x and y lims
    panel.checks <- c("xlim", "ylim")[c("xlim", "ylim") %in% names(ans$call)]
    if(!"xlim" %in% names(ans$panel.args.common))
        ans$panel.args.common$xlim <- ans$x.limits
    if(!"ylim" %in% names(ans$panel.args.common))
        ans$panel.args.common$ylim <- ans$y.limits

#new bit - testing

#make zcase.sublim
#need to rethink labels so 
#zcase.sublim named to say 
#unique, etc if not range
  
    #note: this is pre-rearrange

    if("z" %in% names(ans$panel.args.common) && "zcases" %in% names(ans$panel.args.common)){
        #have inputs to make zcase.sublim
        if(!"zcase.zlim" %in% names(ans$panel.args.common)){
            #for each panel get the range (or unique) of z
            if(!"zcase.ids" %in% ans$panel.args.common){
                ans$panel.args.common$zcase.ids <- if(is.factor(ans$panel.args.common$zcases))
                                                       levels(ans$panel.args.common$zcases) else 
                                                           sort(unique(ans$panel.args.common$zcases))
                #might want to rethink the sort
                #for some cases

#track down where else this is generated
#conmpare outputs

            }
            temp <- lapply(ans$panel.args.common$zcase.ids, function(x) { 
                               temp <- ans$panel.args.common$z[ans$panel.args.common$zcases==x]
                               if(is.numeric(temp)) range(temp) else 
                                   if(is.factor(temp)) levels(temp) else 
                                       unique(temp)
                               #might need to rethink this
                           })
            ans$panel.args.common$zcase.zlim <- temp
        }
    }

    loa.settings <- list()

    if(is.function(panel)) 
        ans$panel <- panel 
        if("loa.settings" %in% names(formals(panel))){
            loa.settings <- panel(loa.settings = TRUE)           
        }

    ignore <- c("xlim", "ylim", "zlim", "xlab", "ylab", "zlab", "main",
                "at", "col.regions", "alpha.regions", "cex.range", "pch.order", 
                "group.elements", "group.ids", "group.args", "scheme", 
                "allowed.scales", "disallowed.scales", "panel.scales",
                "reset.xylims", "load.lists", "lim.borders", "x.elements", 
                "zcase.ids", "zcase.args", "zcase.zlim", "z.rowsum.lim")
    ignore <- unique(c(ignore, loa.settings$common.args))

    temp <- listUpdate(ans$panel.args.common, list(loa.settings = loa.settings, 
                                        ignore = ignore))

    temp <- listUpdate(temp, do.call(groupsAndZcasesHandler, temp))
    ignore <- temp$ignore
    ans$panel.args.common <- temp[!names(temp) %in% c("ignore", "loa.settings")]


#send ignore and ans$panel.args.common
###################################
#
#    #check for groups
#    if("groups" %in% names(ans$panel.args.common)){
#        if(is.null(ans$panel.args.common$group.ids))
#            ans$panel.args.common$group.ids <- if(is.factor(ans$panel.args.common$groups))
#                                                   levels(ans$panel.args.common$groups) else unique(ans$panel.args.common$groups)
#        if(is.null(ans$panel.args.common$group.args))
#            ans$panel.args.common$group.args <- loa.settings$group.args
#        ignore <- unique(c(ignore, ans$panel.args.common$group.args))
#
#        for(i in ans$panel.args.common$group.args)
#            ans$panel.args.common[[i]] <- do.call(groupsHandler, 
#                                                  listUpdate(ans$panel.args.common, 
#                                                             list(z = ans$panel.args.common[[i]], handler=i)))
#    }
#
#    #check for zcases
#    if("zcases" %in% names(ans$panel.args.common)){
#        if(is.null(ans$panel.args.common$zcase.ids))
#            ans$panel.args.common$zcase.ids <- if(is.factor(ans$panel.args.common$zcases))
#                                                   levels(ans$panel.args.common$zcases) else unique(ans$panel.args.common$zcases)
#        if(is.null(ans$panel.args.common$zcase.args))
#            ans$panel.args.common$zcase.args <- loa.settings$zcase.args
#
##to think through
##if zcases and groups are set 
##common args need to be compared
#        if("groups" %in% names(ans$panel.args.common)){
# 
#
##print(
##any(ans$panel.args.common$zcase.args %in% ans$panel.args.common$group.args)
##)
#
#            if(any(ans$panel.args.common$zcase.args %in% ans$panel.args.common$group.args)){
#                temp <- ans$panel.args.common$zcase.args[!ans$panel.args.common$zcase.args %in% ans$panel.args.common$group.args]
#                if(length(temp)<1){
#                    ans$panel.args.common$zcase.args <- NULL
#                    warning("zcase.args/group.args conflict\n", 
#                            "no requested args not already applied by groups", call.=FALSE)
#                } else {
#                    ans$panel.args.common$zcase.args <- temp
#                    warning("zcase.args/group.args conflict\n", 
#                            "one or more requested args already applied by groups", call.=FALSE)
#                }
#           }  
#        }
#        
#        ignore <- unique(c(ignore, ans$panel.args.common$zcase.args))
#
#        for(i in ans$panel.args.common$zcase.args)
#            ans$panel.args.common[[i]] <- do.call(zcasesHandler, 
#                                                  listUpdate(ans$panel.args.common, 
#                                                             list(z = ans$panel.args.common[[i]], handler=i)))
#    }
#
###################################################


#updates to panels

    #what to move to panels

    transfers <- sapply(names(ans$panel.args.common), 
                            function(x) length(ans$panel.args.common[[x]])>1)

    transfers <- names(ans$panel.args.common)[transfers]
    transfers <- ans$panel.args.common[transfers[!transfers %in% ignore]]
    transfers <- listExpand(transfers, ref = ans$panel.args.common$ref)

    #if something to work with
    #transfer subscripted versions to panels
    #wipe these from common
    if(length(transfers)>0){
        temp <- lapply(ans$panel.args, function(y)
                       listUpdate(y, lapply(transfers, function(x) x[y$subscripts])))
        ans$panel.args <- temp
        ans$panel.args.common <- ans$panel.args.common[!names(ans$panel.args.common) %in% names(ans$panel.args[[1]])]
    }

    if(is.null(ans$panel.args.common$panel.elements))
        ans$panel.args.common$panel.elements <- names(ans$panel.args[[1]])

#preprocess

    if(preprocess){

        if(all(c("plot", "process") %in% names(formals(ans$panel)))){


            formals(ans$panel)$plot <- FALSE
            formals(ans$panel)$process <- TRUE 

            ina <- lapply(ans$panel.args, function(x) listUpdate(ans$panel.args.common, x))


#do groups

          out <- if(!"groups" %in% names(ans$panel.args[[1]]) & !"zcases" %in% names(ans$panel.args[[1]]))
                     lapply(ina, function(x) do.call(ans$panel, x)) else {
#                         panel <- if("groups" %in% names(ans$panel.args[[1]]))
#                                      groupsPanelHandler else zcasesPanelHandler
                         panel <- groupsAndZcasesPanelHandler
                         formals(panel)$panel <- ans$panel
                         formals(panel)$plot <- FALSE
                         formals(panel)$process <- TRUE
                         lapply(ina, function(x) do.call(panel, x))
                     }

            ans$panel.args.common$panel.elements <- if(is.null(ans$panel.args.common$panel.elements))
                                                        names(out[[1]]) else 
                                                        unique(c(ans$panel.args.common$panel.elements, names(out[[1]])))

            ans$panel.args <- lapply(1:length(out), function(x)
                                         listUpdate(ans$panel.args[[x]], out[[x]]))

            formals(ans$panel)$plot <- TRUE
            formals(ans$panel)$process <- FALSE 

#move the process set
        }

    }    


#reset x/y lims if wanted
   
##############################################
##############################################
##updated bit
##############################################
##############################################

#could some of this be done at top?
#in checks?

##new bit to allow functions
    if(is.function(reset.xylims)){
        ans <- reset.xylims(ans)
        reset.xylims <- "no.action"
    }
    

##new
    if(is.logical(reset.xylims)){
        reset.xylims <- if(reset.xylims) "refit.xylims" else "no.action" 
    }

##note addition of lim.borders to limsHandler
##this could be tidier


    if("refit.xylims" %in% reset.xylims){
        if(!"xlim" %in% panel.checks){

####################
#testing POSIXct handling
#            temp <- range(lapply(ans$panel.args, 
#                              function(x) range(x$x, na.rm=TRUE, finite=TRUE))
#                    , na.rm=TRUE, finite=TRUE)
            temp <- lapply(ans$panel.args, 
                              function(x) range(x$x, na.rm=TRUE, finite=TRUE))
            temp2 <- do.call(c, temp)
            if("tzone" %in% names(attributes(temp[[1]])))
                attributes(temp2)$tzone <- attributes(temp[[1]])$tzone
            temp <- range(temp2, na.rm=TRUE, finite=TRUE)
#######################

            temp <- limsHandler(x=temp, lim.borders=if(is.null(ans$panel.args.common$lim.borders)) 0.05 else ans$panel.args.common$lim.borders)$xlim
            ans$panel.args.common$xlim <- temp
            ans$x.limits <- temp
        }
        if(!"ylim" %in% panel.checks){

########################
#see above
#            temp <- range(lapply(ans$panel.args, 
#                              function(x) range(x$y, na.rm=TRUE, finite=TRUE))
#                    , na.rm=TRUE, finite=TRUE)
            temp <- lapply(ans$panel.args, 
                              function(x) range(x$y, na.rm=TRUE, finite=TRUE))
            temp2 <- do.call(c, temp)
            if("tzone" %in% names(attributes(temp[[1]])))
                attributes(temp2)$tzone <- attributes(temp[[1]])$tzone
            temp <- range(temp2, na.rm=TRUE, finite=TRUE)
#########################

            temp <- limsHandler(y=temp, lim.borders=if(is.null(ans$panel.args.common$lim.borders)) 0.05 else ans$panel.args.common$lim.borders)$ylim
            ans$panel.args.common$ylim <- temp
            ans$y.limits <- temp
        }
   }


##new
   if("max.xylims" %in% reset.xylims){
        temp<- sqrt(c(ans$x.limits, ans$y.limits)^2)

##if either side of 0? Do the sqrt(x^2)?


        temp <- c(-max(temp), max(temp))
##print(temp)
        ans$panel.args.common$xlim <- temp
        ans$x.limits <- temp
        ans$panel.args.common$ylim <- temp
        ans$y.limits <- temp
   }

#################################################
#################################################

    ###########################
    #update other lims/uniques
    ###########################

#could make this tidier
#could this bit be done better at start?

#compare this and new bit and see if 
#complexity is needed


    #get names of elements in panel.args[[1]] that are numeric
    ranges <- sapply(names(ans$panel.args[[1]]), 
                        function(x) is.numeric(ans$panel.args[[1]][[x]]))
    ranges <- names(ans$panel.args[[1]])[ranges]

    #remove common names 

#need to pull out allowed.scales
#may need to check x/y handling here

    ranges <- ranges[!ranges %in% c("x", "y", "subscripts")]

    #remove elements for which lims are set in panel.args.common
    temp <- gsub("lim$", "", grep("lim$", names(ans$panel.args.common), value=T))
    ranges <- ranges[!ranges %in% temp]

    #recalculate lims for these
    temp.fun <- function(x){
                    range(lapply(ans$panel.args, function(y){
                              range(y[[x]], na.rm=TRUE, finite=TRUE)
                          }), na.rm=TRUE, finite=TRUE)
                }
    if(length(ranges)>0){
        temp <- lapply(ranges, temp.fun)
        names(temp) <- paste(ranges, "lim", sep="")
        ans$panel.args.common <- listUpdate(ans$panel.args.common, temp)
    }

    #repeat for all not characterised
    ranges <- ans$panel.args.common$group.elements
    temp <- gsub("lim$", "", grep("lim$", names(ans$panel.args.common), value=T))
    ranges <- ranges[!ranges %in% temp]

#need to pull out allowed.scales
#may need to check x/y handling here

    ranges <- ranges[!ranges %in% c("x", "y", "subscripts")]

    if(length(ranges)>0){
        temp <- lapply(ranges, function(x)
                     unlist(lapply(ans$panel.args, function(y) 
                         unique(y[[x]]))))
        names(temp) <- paste(ranges, "unique", sep="")
        ans$panel.args.common <- listUpdate(ans$panel.args.common, temp)    
    }


#current allowed.scales fix
    if("allowed.scales" %in% names(ans$panel.args.common))
        ans$panel.args.common <- ans$panel.args.common[!names(ans$panel.args.common) %in% ans$panel.args.common$allowed.scales]
    

##############
#new
##############
    
   if("zlim.in.ylim" %in% reset.xylims & !"ylim" %in% panel.checks){
        if("zlim" %in% names(ans$panel.args.common)){
            temp <- range(c(ans$panel.args.common$zlim, ans$y.limits))
            temp <- limsHandler(y=temp, lim.borders=if(is.null(ans$panel.args.common$lim.borders)) 0.05 else ans$panel.args.common$lim.borders)$ylim
            ans$panel.args.common$ylim <- temp
            ans$y.limits <- temp
        }
   }


    #################
    #add in key
    #################

    if(is.list(legend)){
        legend[[1]]$args$key <- listUpdate(legend[[1]]$args$key, ans$panel.args.common)                 
           
        ans <- update(ans, legend = legend)
    }

    return(ans)

}









##############################
##############################
##loaHandler
##############################
##############################


#old version

#loaHandler <- function(panel = NULL,...){

#    if(is.function(panel)){
#        if("loa.settings" %in% names(formals(panel)))
#            return(panel(loa.settings=TRUE)) else return(FALSE)
#    }

#    return(NULL)

#}

#to test
#new more complex loaHandler

loaHandler <- function(panel = NULL,...){

    if(is.function(panel)){

        if(!"loa.settings" %in% names(formals(panel))) return(FALSE)
        out <- panel(loa.settings=TRUE)
        
        if("data.panel" %in% names(out$default.settings)){

            #if data.panel present
            #update data.panel then panel
            temp <- out$default.settings$data.panel(loa.settings=TRUE)
            out$group.args <- unique(c(out$group.args, temp$group.args))
            out$zcase.args <- unique(c(out$zcase.args, temp$zcase.args))
            out$common.args <- unique(c(out$common.args, temp$common.args))
            if("load.lists" %in% names(out$default.settings) & 
               "load.lists" %in% names(temp$default.settings))            
                   out$default.settings$load.lists <- 
                      unique(c(out$default.settings$load.lists, 
                               temp$default.settings$load.lists))
            out$default.settings <- listUpdate(temp$default.settings, 
                                               out$default.settings)

        }

        return(out) 
    }

    return(NULL)

}





