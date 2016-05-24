#in development code
#[TBC - NUMBER] functions 


#getArgs
#getPlotArgs
#isGood4LOA
#parHandler


#NOTE: much borrowed from lattice 

#to do
##############################
#fix parHandler error
#colHandler return col vector by default
#set output="col" 
#might be mutliples
##############################
#callWithThis document or drop
#very minor






##############################
##############################
##getArgs
##############################
##############################



getArgs <- function (source = TRUE, local.resets = TRUE, user.resets = TRUE, 
    is.scales.lines = FALSE, elements = NULL, ..., defaults = list(), 
    defaults.only = FALSE) 
{


    targets <- if (defaults.only) 
        c(names(defaults), "isGood4LOA")
    else NULL
    fix.source <- function(n, a = 1, b = 1, sc.el, is.sc) {
        if (is.null(n)) 
            return(list(isGood4LOA = FALSE))
        if (is.logical(n)) 
            if (n) 
                return(list())
            else return(list(isGood4LOA = FALSE))
        if (is.list(n)) {
            if (!is.null(sc.el) && sc.el %in% names(n)) {
                if (is.vector(n[[sc.el]]) & !is.list(n[[sc.el]])) {
                  if (length(n[[sc.el]]) < b) 
                    n[[sc.el]] <- zHandler(n[[sc.el]], TRUE, 
                      1:b)
                  n[[sc.el]] <- list(col = n[[sc.el]][seq(a, 
                    length(n[[sc.el]]), by = b)])
                }
                n <- listUpdate(n, n[[sc.el]])
                n <- n[names(n) != sc.el]
                if (is.sc) {
                  n <- n[!names(n) %in% c("col", "alpha")]
                  if ("col.line" %in% names(n)) 
                    names(n)[names(n) == "col.line"] <- "col"
                  if ("alpha.line" %in% names(n)) 
                    names(n)[names(n) == "alpha.line"] <- "alpha"
                }
            }
            return(n)
        }
        if (is.data.frame(n)) 
            return(as.data.frame(n))
        if (length(n) < b) 
            n <- zHandler(n, TRUE, 1:b)
        list(col = n[seq(a, length(n), by = b)])
    }
    el.handler <- function(def, sc.el, a = 1, b = 1, so, lr, 
        ur, is.sc = FALSE) {
        output <- listUpdate(def, fix.source(so, sc.el = sc.el, 
            a = a, b = b, is.sc = is.sc), use.b = targets)
        output <- listUpdate(output, fix.source(lr, sc.el = sc.el, 
            a = a, b = b, is.sc = FALSE), use.b = targets)


        output <- listUpdate(output, fix.source(ur, sc.el = sc.el, a = a, 
            b = b, is.sc = FALSE), use.b = targets)

        if(!is.null(output$col)){
            if(is.character(output$col) && output$col[1] == "transparent")
                if(!is.null(output$local.col)) output$col <- output$local.col
        }
        output
    }

    if (is.null(elements)) 
        return(el.handler(defaults, NULL, 1, 1, source, local.resets, 
            user.resets, is.scales.lines))
    if (length(elements) < 2){ 
        output <- el.handler(defaults, elements, 1, 1, source, local.resets, 
            user.resets, is.scales.lines)

return(output)

}
    output <- lapply(1:length(elements), function(x) el.handler(defaults, 
        elements[x], x, length(elements), source, local.resets, 
        user.resets, is.scales.lines))
    names(output) <- elements

    output
}







##############################
##############################
##getPlotArgs
##############################
##############################



getPlotArgs <- function(defaults.as = "axis.line", source = TRUE, local.resets = TRUE, user.resets = TRUE, 
                        elements = NULL, ..., is.scales.lines = NULL, defaults.only = TRUE){



    #check defaults.as
    if(!defaults.as[1] %in% names(trellis.par.get())){
         warning("could not recover 'defaults.as' [guessing axis.text]", call.=FALSE)
         defaults.as <- "axis.text"
    }

    if(is.null(is.scales.lines))
        is.scales.lines <- if(is.character(defaults.as)  && length(grep("line", defaults.as)) > 0)
                               TRUE else FALSE

    #get defaults
    defaults <- trellis.par.get(defaults.as[1])

    ans <- getArgs(source = source, local.resets = local.resets, user.resets = user.resets, 
                   is.scales.lines = is.scales.lines, defaults.only = defaults.only, defaults = defaults, 
                   elements = elements)

    ans
}




##############################
##############################
##isGood2Plot
##############################
##############################

#temp function

#arg could be replaced with ...
#function could be recursive
#and would be needed

isGood4LOA <- function(arg){

    if(is.null(arg)) return(FALSE)
    if(is.logical(arg) && !all(arg)) return(FALSE)
    if(is.list(arg) && !is.null(arg$isGood4LOA) && !arg$isGood4LOA) return(FALSE)
    TRUE

}





########################
########################
##parHandler
########################
########################

parHandler <- function (scheme = NULL, ...) {

#######################
#users
#######################
#loaPlot
# 


#######################
#urgent
#######################
#rethink the inputs
#

#think about order
#theme first means unnamed get assigned as this

#could hide this in ...
#then
#if theme named
#if output set par.settings
#currently only allowed output is par.settings
#could rethink using output option
#

#could use theme and set.theme
#then if garbage sent it will ignore

    extra.args <- list(...)

#in loaPlot
#    extra.args <- do.call(listLoad, listUpdate(extra.args, list(load="par.settings")))
    if(is.null(scheme) && is.null(extra.args$par.settings))
         return(NULL)
    if(!is.null(extra.args$par.settings) && "scheme" %in% names(extra.args$par.settings))
           scheme <- extra.args$par.settings$scheme


    ###################
    #pre-defined themes
    ###################

    myschm <- NULL

    if(is.character(scheme)){

        if (scheme == "greyscale") {
            symbol <- gray(1:8/8)
            fill <- "grey"
            region <- gray(11:1/11)
            reference <- "black"
            bg <- "transparent"
            fg <- "black"
            myschm <- list(plot.polygon = list(col = fill[1], border = fg[1]), 
                         box.rectangle = list(col = symbol[1]), box.umbrella = list(col = symbol[1]), 
                         dot.line = list(col = reference), dot.symbol = list(col = symbol[1]), 
                         plot.line = list(col = symbol[1]), plot.symbol = list(col = symbol[1]), 
                         regions = list(col = colorRampPalette(region)(100)), 
                         reference.line = list(col = reference), superpose.line = list(col = symbol), 
                         superpose.symbol = list(col = symbol), superpose.polygon = list(col = fill, 
                         border = fg), background = list(col = bg), add.line = list(col = fg), 
                         add.text = list(col = fg), box.dot = list(col = fg), 
                         axis.line = list(col = fg), axis.text = list(col = fg), 
                         strip.border = list(col = fg), strip.background = list(col = "white"), 
                         strip.shingle = list(col = "grey"), box.3d = list(col = fg), 
                         par.xlab.text = list(col = fg), par.ylab.text = list(col = fg), 
                         par.zlab.text = list(col = fg), par.main.text = list(col = fg), 
                         par.sub.text = list(col = fg))
        }

        if (scheme == "kr.web") {
            symbol <- colHandler(1:10, col.regions="Blues", output = "col")[3:10]
            fill <- "white"
            region <- colHandler(1:14, col.regions="Blues", output = "col")[3:14]
            reference <- "white"
            bg <- "black"
            fg <- "white"
            myschm <- list(plot.polygon = list(col = fill[1], border = fg[1]), 
                         box.rectangle = list(col = symbol[1]), box.umbrella = list(col = symbol[1]), 
                         dot.line = list(col = reference), dot.symbol = list(col = symbol[1]), 
                         plot.line = list(col = symbol[1]), plot.symbol = list(col = symbol[1]), 
                         regions = list(col = colorRampPalette(region)(100)), 
                         reference.line = list(col = reference), superpose.line = list(col = symbol), 
                         superpose.symbol = list(col = symbol), superpose.polygon = list(col = fill, 
                         border = fg), background = list(col = bg), add.line = list(col = fg), 
                         add.text = list(col = fg), box.dot = list(col = fg), 
                         axis.line = list(col = fg), axis.text = list(col = fg), 
                         strip.border = list(col = fg), strip.background = list(col = "black"), 
                         strip.shingle = list(col = "lightgrey"), box.3d = list(col = fg), 
                         par.xlab.text = list(col = fg), par.ylab.text = list(col = fg), 
                         par.zlab.text = list(col = fg), par.main.text = list(col = fg), 
                         par.sub.text = list(col = fg))
        }

        if (scheme == "kr.blues") {
            symbol <- colHandler(1:10, col.regions="Blues", output = "col")[3:10]
            fill <- "darkblue"
            region <- colHandler(1:14, col.regions="Blues", output = "col")[3:14]
            reference <- "darkblue"
            bg <- "white"
            fg <- "darkblue"
            myschm <- list(plot.polygon = list(col = fill[1], border = fg[1]), 
                         box.rectangle = list(col = symbol[1]), box.umbrella = list(col = symbol[1]), 
                         dot.line = list(col = reference), dot.symbol = list(col = symbol[1]), 
                         plot.line = list(col = symbol[1]), plot.symbol = list(col = symbol[1]), 
                         regions = list(col = colorRampPalette(region)(100)), 
                         reference.line = list(col = reference), superpose.line = list(col = symbol), 
                         superpose.symbol = list(col = symbol), superpose.polygon = list(col = fill, 
                         border = fg), background = list(col = bg), add.line = list(col = fg), 
                         add.text = list(col = fg), box.dot = list(col = fg), 
                         axis.line = list(col = fg), axis.text = list(col = fg), 
                         strip.border = list(col = fg), strip.background = list(col = "azure"), 
                         strip.shingle = list(col = "lightblue"), box.3d = list(col = fg), 
                         par.xlab.text = list(col = fg), par.ylab.text = list(col = fg), 
                         par.zlab.text = list(col = fg), par.main.text = list(col = fg), 
                         par.sub.text = list(col = fg))
        }



    }

    if(is.list(scheme))
        myschm <- scheme

####################
#temp fix
####################
    
     return(listUpdate(myschm, extra.args$par.settings))

####################
    
    if(is.list(myschm)){
        if (is.null(extra.args$par.settings)) 
            extra.args$par.settings <- myschm
        else extra.args$par.settings[!names(myschm) %in% names(extra.args$par.settings)] <- myschm[!names(myschm) %in% 
            names(extra.args$par.settings)]
    }

#warning if myschm not list?

#return par.setting or extra.args?

    extra.args$par.settings
}











