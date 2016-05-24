#in development code
#[1 -TBC] functions 



#######################
##might want to make 
##own space for conds
#######################


#formulaHandler - handles the x formula
#formulaHandler.old - previous version 

#matrixHandler - converts matrx to x formula

#####not exported
#stripHandler   - handles strips
#getZcaseDimensions

#urgent
##########################
#stripHandler 
##########################
#temp fix for cond 
#need a better fix
#

#fixes 
#to dos
#suggestions


###########################
###########################
#formulaHandler
###########################
###########################

formulaHandler <- function (x, data = NULL, groups = NULL, ..., 
    expand.plot.args = TRUE, formula.type = "z~x*y|cond", panel.zcases = FALSE, 
    coord.conversion = NULL, lattice.like = NULL, check.xy.dimensions = TRUE, 
    check.coord.dimensions = TRUE, get.zcase.dimensions = TRUE, output = "extra.args") 
{

#test new#

    extra.args <- list(...)

###################
##new  1
###################

    if(is.null(lattice.like)){

#think about this
#might be better other way around

#NOTE: all temp.fun are given unique names 
#      for CRAN submission/rcmd checks


        if (!"loa.err.message" %in% names(extra.args)) 
            extra.args$loa.err.message <- "problem with x formula/data combination"

    #get formula info

#could simplify this?

        temp <- strsplit(as.character(deparse(formula.type)), "~")
        temp <- temp[[1]][length(temp[[1]])]
        temp <- strsplit(temp, "[|]")[[1]][1]
        coord.sep <- if(length(grep("[+]", temp))>0) "[+]" else "[*]"
        coord.cases <- gsub(" ", "", strsplit(temp, coord.sep)[[1]])
        coord.cases <- gsub("\"", "", coord.cases)

     
#length of coord.cases should be length of coords
#as check for check...

#####################################
#from get_all_vars, model.frame, etc.
#####################################

#err messages need tidying

        if (missing(x)) {
            if (!missing(data) && inherits(data, "data.frame") && 
                length(attr(data, "terms"))) 
                return(data)
            x <- as.formula(data)
        }
        else if (missing(data) && inherits(x, "data.frame")) {
            if (length(attr(x, "terms"))) 
                return(x)
            data <- x
            x <- as.formula(data)
        }
        x <- as.formula(x)
        if (missing(data)) 
            data <- environment(x)
        else if (!is.data.frame(data) && !is.environment(data) && 
            !is.null(attr(data, "class"))) 
            data <- as.data.frame(data)
        else if (is.array(data)) 
            stop("'data' must be a data.frame, not a matrix or an array")

        #########################################
        #knocked out to allow 1:10~1:10*1:10, etc
        #if (!inherits(x, "terms")) 
        #    x <- terms(x, data = data)
        #########################################

        env <- environment(x)
        rownames <- .row_names_info(data, 0L)
        varnames <- all.vars(x)

#######################################


        ####################
        #get variables
        ####################

#maybe rebuild the next bit
#case condition
        temp <- as.character(x)
        if(length(temp)==3){
               test1 <- temp[2]
               test2 <- temp[3]
        } else {
               test1 <- NULL
               test2 <- temp[2]
        }
        temp <- strsplit(test2, "[|]")[[1]]
        if(length(temp)>1){
            test2 <- temp[1]
            test3 <- temp[2]
        } else test3 <- NULL

        temp.fun <- function(x, str="[+]", unique=TRUE) {
                        if(is.null(x)) return(x)
                        x = strsplit(gsub(" ", "", x), str)[[1]]
                        ans <- lapply(x, function(x){
                                    x <- parse(text=x)
                                    try(eval(x, data, env), silent=TRUE)})
#make.names unique
#BUT allow 1 not x1
                        names(ans) <- if(unique) make.unique(x) else x
                        ans
                        } 
        zcases <- temp.fun(test1)
        coords <- temp.fun(test2, str=coord.sep, unique=FALSE)
        conds <- temp.fun(test3)

        #any missing
        #if so warn about all
        #not just the first encountered

        temp.fun2 <- function(x){
                        temp <- sapply(x, function(x) class(x)[[1]]=="try-error")
                        names(x)[temp]
                        }
        test1 <- temp.fun2(zcases)
        test2 <- temp.fun2(coords)
        test3 <- temp.fun2(conds)

#could put this in temp.fun2
#further reduce this
        if(length(c(test1, test2, test3))>0){
              reply <- extra.args$loa.err.message
              if(length(test1)>0) reply <- paste(reply, paste("[unknown z terms: ", paste(test1, collapse=", "), "]", sep=""), sep="\n       ")
              if(length(test2)>0) reply <- paste(reply, paste("[unknown coord terms: ", paste(test2, collapse=", "), "]", sep=""), sep="\n       ")
              if(length(test3)>0) reply <- paste(reply, paste("[unknown cond terms: ", paste(test3, collapse=", "), "]", sep=""), sep="\n       ")
              stop(reply, call.=FALSE)
        }

###########################

        ############################
        #check.dimensions
        ############################

#might need to move this
#think about add index
#if coord.cases is 2 and only one coord    
#could simply this when check.xy.dimensions droppped

        if(length(coords)!=length(coord.cases)){
            reply <- "coordinate dimensions do not match assigned formula.type"
            if(length(coord.cases)==2 && all(check.xy.dimensions, check.coord.dimensions))
                stop(reply, call.=FALSE)
            if(length(coord.cases)!=2 && check.coord.dimensions)
                stop(reply, call.=FALSE)
        }
    
        ###########################
        #pad zcases, coords, conds
        ###########################
        if(expand.plot.args){
            temp.fun3 <- function(x) if(is.null(x)) x else 
                                        max(sapply(x, length))
            temp <- max(c(temp.fun3(zcases), temp.fun3(coords), temp.fun3(conds)))
        
            temp.fun4 <- function(x){
                            if(is.list(x))
                                for(i in 1:length(x)){
                                    if(!is.null(x[[i]]))
                                         x[[i]] <- rep(x[[i]], length.out=temp)
                                }
                             x} 


#this is the bit that gives you the current wrapping numbers behaviour
#this is where the alternatives would go in
#1:10 ~ 1 -> 1:10~ rep(1,10)
#         -> 1:10~ c(1, rep(NA, 9))
#         -> 1 ~ 1

#I think this only reworks x, y, etc for z1 ~ y * x
#NOT z1+z2~y (y1 +y2, wrapped) + x (x1 + x2, wrapped)
#I think that is done below when zcases are made

            zcases <- temp.fun4(zcases)
            coords <- temp.fun4(coords)
            conds <- temp.fun4(conds)

            if(!is.null(groups))
                 groups <- temp.fun4(list(groups))[[1]]

        }
    

        ############################
        #make zcases
        ############################

        t1 <- length(zcases)
        t2 <- sapply(zcases, length)

#make this case conditioning 
#z <- as.vector(unlist(zcases)) works for all
#currently not using zcase.args
##think the equivalent is actuall zcase.ids anyway!
#could simply this

        zlab <- names(zcases)
        if(!is.null(zlab))
             zlab <- paste(zlab, collapse=" + ")

        if(t1==0) {
            z <- zcases
            zcases <- NULL
            zcase.args <- NULL
        }
        if(t1==1){ 
            z <- zcases[[1]]
            zcases <- NULL
            zcase.args <- NULL
        }
        if(t1>1){
            z <- as.vector(unlist(zcases))
            zcase.args <- names(zcases)
            temp <- lapply(1:length(zcases), function(x){
                                 rep(names(zcases[x]), t2[x])
                           })


#############################
#new bit 0.2.26
#############################
#self ordering factor in previous version
#level added to order based on supplied order
#############################
            zcases <- factor(as.vector(unlist(temp)), levels = unique(as.vector(unlist(temp))))

            temp.fun5 <- function(x){
                            if(is.list(x)) {
                                for(i in 1:length(x)){
                                    temp <- lapply(1:length(t2), function(y)
                                                        rep(x[[i]], length.out=t2[y]))
###############################
#new bit 0.20.28
###############################
#handle posixct better
#                                    x[[i]] <- as.vector(unlist(temp))
                                     x[[i]] <- do.call(c, temp)
                                     if("tzone" %in% names(attributes(temp[[1]])))
                                          attributes(x[[i]])$tzone <- attributes(temp[[1]])$tzone
###############################
                                }
                                x
                            } else return(x)
            }

            coords <- temp.fun5(coords)
#shingle in previous version
            conds <- temp.fun5(conds)

            if(!is.null(groups))
                groups <- temp.fun5(list(groups))[[1]]
        }

        ###############################
        #panel.zcases
        ###############################
    
        if(!is.null(zcases) && panel.zcases){
            conds <- listUpdate(conds, list(zcases = zcases))
            zcases <- NULL
        }

        ###############################
        #make lattice.like output
        ###############################

#could simplify this
#bundle listUpdates

        lattice.like <- coords
        names(lattice.like) <- coord.cases
        temp <- as.list(names(coords))
        names(temp) <- paste(coord.cases, "lab", sep="")
        lattice.like <- listUpdate(lattice.like, temp)
        lattice.like <- listUpdate(lattice.like, list(z=z, zlab=zlab, zcases=zcases))
        lattice.like <- listUpdate(lattice.like, list(groups=groups))
        lattice.like <- listUpdate(lattice.like, list(panel.condition=conds))

        if(get.zcase.dimensions)
            lattice.like <- do.call(getZcaseDimensions, lattice.like)

    
        ##################################
        #return if output lattice.like
        ##################################

        if(output=="lattice.like") return(lattice.like)

##################
##new 1
##just }
##################

    }

    ##################################
    #extra.args updates
    ##################################


#################
###new 2
#################

#allow for coord.conversion
#others >> x,y
#before plot call

#    if(!is.null(coord.conversion)){
#        lattice.like <- listUpdate(lattice.like, do.call(coord.conversion, lattice.like))
#    }

#new fix BUT this needs revisiting

#    if (!is.null(coord.conversion)) {
#        temp <- extra.args[names(extra.args) %in% names(formals(coord.conversion))]
#        lattice.like <- listUpdate(listUpdate(lattice.like, temp), do.call(coord.conversion, 
#            lattice.like))
#    }


#new fix from version 0.2.28
#(stackPlot introduction)
#this may still need revisiting

#    if (!is.null(coord.conversion)) {
#        temp <- extra.args[names(extra.args) %in% names(formals(coord.conversion))]
#        lattice.like <- do.call(coord.conversion, listUpdate(lattice.like, temp))
#    }

#new fix from version 0.2.29
#(z/zcase realignment after above 'fix')
#this may still need revisiting

    if (!is.null(coord.conversion)) {
        temp <- extra.args[names(extra.args) %in% names(formals(coord.conversion))]
        temp <- do.call(coord.conversion, listUpdate(lattice.like, 
            temp))
        lattice.like <- listUpdate(lattice.like, temp)
    }


#################
#this makes some of 
#below redundant
#(if it works)
################

#############


###############
#testing return labels
###############

#this might now be redundant

     temp <- listUpdate(lattice.like, extra.args, use=c("xlab", "ylab", "zlab"))
     lattice.like <- listUpdate(lattice.like, temp)


#this might now be redundant

    extra.args <- listUpdate(extra.args, lattice.like, 
                          ignore.b = c("panel.condition", "x", "y"))

    extra.args <- do.call(stripHandler, listUpdate(list(striplab = names(lattice.like$panel.condition)), 
        extra.args))
    ..loa.x <- lattice.like$x
    ..loa.y <- lattice.like$y
    extra.args$z <- lattice.like$z
    extra.args$ref <- lattice.like$x
    extra.args$groups <- lattice.like$groups

    extra.args <- listUpdate(lattice.like, extra.args, use.a=c("xlab","ylab","zlab"))
#    if ("zcases" %in% names(lattice.like)) 
#        extra.args$zcases <- lattice.like$zcases
    extra.args <- listUpdate(extra.args, lattice.like, use.b=c("zcases", "zcase.ids", "zcase.zlim", "z.rowsum.lim"))

    x <- "..loa.y~..loa.x"
    if (!is.null(lattice.like$panel.condition) && length(lattice.like$panel.condition) > 
        0) {
        ..loa.cond <- lattice.like$panel.condition
        temp <- paste("..loa.cond[[", 1:length(..loa.cond), sep = "")
        temp <- paste(temp, "]]", sep = "", collapse = "+")
        x <- paste(x, temp, sep = "|")
    }
    extra.args$x <- as.formula(x)

    return(extra.args)

}







#############################
#############################
##matrixHandler
#############################
#############################


#this is based on levelplot.matrix in lattice

#started 
#kr 26/04/2015

matrixHandler <- function (x, data = NULL, row.values=NULL, 
                           column.values=NULL, ...){

    extra.args <- list(...)

    if(is.null(row.values)) row.values <- seq_len(nrow(x))
    if(is.null(column.values)) column.values <- seq_len(ncol(x))

###tidy
    stopifnot(length(row.values) == nrow(x), length(column.values) == 
        ncol(x))
    
##tidy 
    if (!is.null(data)) 
        warning("supplied 'data' ignored; x matrix content used")


    form <- z ~ row * column
    data <- expand.grid(row = row.values, column = column.values)
    data$z <- as.vector(as.numeric(x))


#    if (!"xlim" %in% names(extra.args)) 
#        extra.args$xlim <- if (!is.null(rownames(x))) 
#            rownames(x)
#        else range(row.values, finite = TRUE) + c(-0.5, 0.5)
#    if (!"ylim" %in% names(extra.args)) 
#        extra.args$ylim <- if (!is.null(colnames(x))) 
#            colnames(x)
#        else range(column.values, finite = TRUE) + c(-0.5, 0.5)

    if (!"xlim" %in% names(extra.args)) 
        extra.args$xlim <- range(row.values, finite = TRUE) + c(-0.5, 0.5)
    if (!"ylim" %in% names(extra.args)) 
        extra.args$ylim <- range(column.values, finite = TRUE) + c(-0.5, 0.5)

     listUpdate(list(x=form, data=data), extra.args)
}






















##########################
##########################
##triFormulaHandler
##########################
##########################

#triABCFormulaHandler <- function(x, data = NULL, ..., formula.type="z~a+b+c|cond"){

#    extra.args <- list(...)
#    extra.args$output <- "lattice.like"

#    temp <- formulaHandler(x=x, data=data, formula.type=formula.type, output="lattice.like")
    
#    temp


#}




#triABCFormulaHandler(1~1:2+3:4+5:7|8:9)


####see notes in formulaHandler for what makes this wrap? and should it???

###see expand.plot.args























###########################
###########################
#formulaHandler.old
###########################
###########################

formulaHandler.old <- function(x, data = NULL, ..., formula.type="z~x*y|cond", 
                           panel.zcases = FALSE,
                           check.xy.dimensions=TRUE, output = "extra.args"){

    #extra.args
    extra.args <- list(...)

    if(!"loa.err.message" %in% names(extra.args))
        extra.args$loa.err.message <- "problem with x/data combination"

    #get base equation and z terms

    allowed.args <- names(formals(latticeParseFormula))
    temp <- listUpdate(list(dimension=3, multiple=TRUE, outer=TRUE), 
                             extra.args, use.b=allowed.args)
    temp <- temp[names(temp) != "model"]
    temp <- listUpdate(list(model=x, data=data), temp)

    d1 <- try(do.call(latticeParseFormula, temp), silent = TRUE)
    if(is(d1)[1] == "try-error")
        stop(extra.args$loa.err.message, call. = FALSE)

##we get warning here if a duplicated term 
##is used in formula z component
##e.g. z1 + z1 ~ a * b in 
##suppressWarnings would hide this


##next bit is needed for end tests
##if not using anywhere else anymore
##could do better/combine at top

    #the end conditioning
    temp <- as.character(x)
    temp <- temp[length(temp)]
    cond <- strsplit(temp, "[|]")[[1]]
#    if(length(cond)>1)
#        d1$panel.condition <- cond[2]    



####################
#tidy these
####################

    #the z, z condition
    names(d1)[names(d1)=="left"] <- "z"
    names(d1)[names(d1)=="left.name"] <- "z.name"
    z <- d1$left


##this fails if cond a>0 + b>0
##but not g + g even when g = a>0
##but so lattice xyplot!
##see what panel.condition is then

##also z1+z2 does not set zcases


    if(!is.null(d1$condition)){
        if(is.null(names(d1$condition)))
            d1$zcases <- d1$condition[[1]] else 
                if("" %in% names(d1$condition))
                    d1$zcases <- d1$condition[[which(names(d1$condition) == "")]]

        if(is.factor(d1$zcases)){
            levels(d1$zcases) <- make.unique(levels(d1$zcases))
            d1$z.name <- paste(levels(d1$zcases), collapse = " + ")
        }
        
        if(any(names(d1$condition) != ""))
            d1$panel.condition <- d1$condition[names(d1$condition) != ""]

    }
    d1$condition <- NULL


##    if(length(d1$condition)>0){
##        if(is.null(names(d1$condition)))
##            names(d1)[names(d1)=="condition"] <- "zcases" else
##                if(length(d1$condition) > length(names(d1$condition)[names(d1$condition)!=""]))
##                    d1$zcases <- d1$condition[names(d1$condition)==""]
##        for(i in 1:length(d1$zcases)){
#testing make these unique
##            levels(d1$zcases[[i]]) <- make.unique(levels(d1$zcases[[i]]))
##        }
#testing remove list, so it vectorises
##        d1$zcases <- d1$zcases[[1]]

##    }

    #x,y
    if(formula.type=="z~y*x|cond"){
        names(d1)[names(d1)=="right.y"] <- "x"
        names(d1)[names(d1)=="right.y.name"] <- "x.name"
        names(d1)[names(d1)=="right.x"] <- "y"
        names(d1)[names(d1)=="right.x.name"] <- "y.name"
    } else {
        names(d1)[names(d1)=="right.x"] <- "x"
        names(d1)[names(d1)=="right.x.name"] <- "x.name"
        names(d1)[names(d1)=="right.y"] <- "y"
        names(d1)[names(d1)=="right.y.name"] <- "y.name"
    }


#######################
#check for extra dimensions in x, y
#currently x or y not both!
#if allowed?
#######################
#may rethink this
#move to top?


    if(check.xy.dimensions==TRUE){
        temp <- cond[1]
        temp <- strsplit(temp, "[+]|[*]")[[1]]
        if(length(temp) > 2)
            stop("multiple 'x' and/or 'y' dimensions currently not allowed", call. = FALSE)
    }

#######################
#maybe work into extra.args
#in future version
#######################

    #structure zcases
    if("zcases" %in% names(d1) && panel.zcases){
        d1$panel.condition <- listUpdate(list(zcases = d1$zcases), d1$panel.condition)
        d1$zcases <- NULL
    }

    #export results if lattice.like
    if(output=="lattice.like") return(d1)

    #temp fix for conditioning labels
    extra.args <- do.call(stripHandler,
                          listUpdate(list(striplab = names(d1$panel.condition)), extra.args)
                         )


    ..loa.x <- d1$x
    ..loa.y <- d1$y

     extra.args$z <- d1$z
     extra.args$ref <- d1$x
     extra.args$groups <- d1$groups
     extra.args <- listUpdate(list(xlab = d1$x.name, ylab = d1$y.name, zlab = if(is.null(extra.args$z)) NULL else d1$z.name),
                              extra.args)

    

    if("zcases" %in% names(d1))
        extra.args$zcases <- d1$zcases


    x <- "..loa.y~..loa.x"
    if(!is.null(d1$panel.condition) && length(d1$panel.condition)>0){
        ..loa.cond <- d1$panel.condition
        temp <- paste("..loa.cond[[" , 1:length(..loa.cond), sep="")
        temp <- paste(temp, "]]", sep="", collapse="+")
        x <- paste(x, temp, sep="|")

    }

    
    extra.args$x <- as.formula(x)

    extra.args
 
}





####################################
####################################
##stripHandler
####################################
####################################


stripHandler <- function(..., striplab=NULL){

##########################
#messy 
#needs a rethink
##########################

##########################
#pass list 
#with listLoad handling
##########################

    extra.args <- list(...)

    #if not striplab
    #nothing to do
    if(is.null(striplab)) return(extra.args)

    if(!"strip" %in% names(extra.args)) extra.args$strip <- TRUE 

    if("strip" %in% names(extra.args)){

        if(is.logical(extra.args$strip) && !extra.args$strip) return(extra.args)

        temp <- if(is.function(extra.args$strip)) extra.args$strip else strip.default

        extra.args$strip <- function(var.name, ...) temp(..., var.name=striplab)

    }

    return(extra.args)
}




##############################
##############################
##getZcaseDimensions
##############################
##############################


getZcaseDimensions <- function(...){

    #######################
    #calculate zcase dimensions if missing
    #######################
    
    extra.args <- list(...)

    if ("z" %in% names(extra.args) && "zcases" %in% names(extra.args)) {

        if (!"zcase.ids" %in% extra.args) {
            extra.args$zcase.ids <- if (is.factor(extra.args$zcases)) 
              levels(extra.args$zcases)
            else sort(unique(extra.args$zcases))
        }


        z <- extra.args$z

#return(z)
#return(extra.args$zcases)

        if(!is.list(z))
             z <- sapply(extra.args$zcase.ids, function(x) z[extra.args$zcases==x], 
                         simplify = FALSE)

#return(z)

        
        if(!"zcase.zlim" %in% names(extra.args)) {
             extra.args$zcase.zlim <- sapply(z, function(x) {
                                                if (is.numeric(x)) range(x) else
                                                if (is.factor(x)) levels(x) else 
                                                unique(x)}, simplify = FALSE)
        }
        if(!"z.rowsum.lim" %in% names(extra.args)) {
            extra.args$z.rowsum.lim <- range(
                  sapply(1:length(z[[1]]), function(i)
                       sum(sapply(z, function(x) x[i]), na.rm=TRUE)), 
                       na.rm=TRUE)
        }
    }

    extra.args
    
}

