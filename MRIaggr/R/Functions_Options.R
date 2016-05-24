MRIaggr.env <- new.env()

assign(".ls_optionsMRIaggr", 
       list(asp = 1, 
            axes = TRUE, 
            bg = "lightblue", 
            breaks = 50, 
            cex = 1, 
            cex.index = c(1, 1, 1), 
            cex.legend = 1.5, 
            cex.main = 1.5, 
            checkArguments = TRUE,
            col.index = c("red","purple","green"), 
            col.midplane = "red",
            col.NA = "lightyellow", 
            digit.legend = 2,
            digit.result = 2,
            digit.epsilon = 5,
            digit.percentage = 3,
            filter.index = "2D_N4",
            height = 500, 
            hemisphere = "both", 
            legend = TRUE, 
            main = NULL, 
            main.legend = NULL, 
            mar = rep(1.5, 4), 
            mar.legend = c(2, 7, 2, 2), 
            mfrow = NULL, 
            mgp = c(2, 0.5, 0), 
            norm_mu = FALSE, 
            norm_sigma = FALSE, 
            num.main = TRUE, 
            numeric2logical = FALSE, 
            outline.index = FALSE, 
            palette = "terrain.colors", 
            path = NULL, 
            pch.index = 20:22, 
            pch.NA = 8, 
            pty = NULL, 
            quantiles.legend = TRUE, 
            res = NA,
            slice_var = "k", 
            type.breaks = "range", 
            unit = "px", 
            verbose = TRUE,
            xlab = "", 
            ylab = "", 
            width = 1000, 
            window = FALSE
       ), 
       envir = MRIaggr.env)


optionsMRIaggr <- function(..., reinit.options = FALSE){
  
  if(reinit.options == TRUE){
    
    assign(".ls_optionsMRIaggr", 
           list(asp = 1, 
                axes = TRUE, 
                bg = "lightblue", 
                breaks = 50,
                cex = 1, 
                cex.index = c(1,1,1), 
                cex.legend = 1.5, 
                cex.main = 1.5,
                checkArguments = TRUE,
                col.index = c("red", "purple", "green"), 
                col.midplane = "red", 
                col.NA = "lightyellow",
                digit.legend = 2,
                digit.result = 2,
                digit.epsilon = 5,
                digit.percentage = 3,
                filter.index = "2D_N4", 
                height = 500, 
                hemisphere = "both",
                legend = TRUE, 
                main = NULL, 
                main.legend = NULL, 
                mar = rep(1.5, 4), 
                mar.legend = c(2,7,2,2), 
                mfrow = NULL, 
                mgp = c(2,0.5,0), 
                norm_mu = FALSE,
                norm_sigma = FALSE, 
                num.main = TRUE, 
                numeric2logical = FALSE, 
                outline.index = FALSE,
                palette = "terrain.colors", 
                path = NULL, 
                pch.index = 20:22, 
                pch.NA = 8, 
                pty = NULL, 
                quantiles.legend = TRUE, 
                res = NA, 
                slice_var = "k",
                type.breaks = "range", 
                unit = "px", 
                verbose = TRUE,
                xlab = "", 
                ylab = "", 
                width = 1000, 
                window = FALSE),               
           envir = MRIaggr.env     
    )
    
    return(invisible(get(".ls_optionsMRIaggr", envir = MRIaggr.env)))
  }
  
  args <- list(...)
  n.args <- length(args)
  
  #### si lecture 
  if(n.args == 0){ # retourne tout
    
    value <- selectOptionsMRIaggr()
    single <- FALSE
    
  }else if(all(unlist(lapply(args, is.character)))){ # retourne uniquement les arguments demandes
    
    args <- unlist(args)
    if (length(args) == 1) {single <- TRUE} else {single <- FALSE}
    value <- selectOptionsMRIaggr(args)
    
  } else { #### si ecriture 
    if (length(args) == 1) {single <- TRUE} else {single <- FALSE}
    value <- allocOptionsMRIaggr(args, n.args)
  }
  
  #### export
  if (single){ 
    value <- value[[1L]]
  }
  
  if(!is.null(names(args))){
    invisible(value)
  }else{
    value
  } 
}

selectOptionsMRIaggr <- function(field = NULL){
  
  .ls_optionsMRIaggr <- get(".ls_optionsMRIaggr", envir = MRIaggr.env)
  
  if(is.null(field)){
    
    return(.ls_optionsMRIaggr)
    
  }else{
    
    validCharacter(value = field, validLength = NULL, validValues = names(.ls_optionsMRIaggr), method = "selectOptionsMRIaggr")
    
    return(.ls_optionsMRIaggr[field])
  }
  
}

allocOptionsMRIaggr <- function(field, n.args){
  
  names.field <- names(field)
  .ls_optionsMRIaggr <- get(".ls_optionsMRIaggr", envir = MRIaggr.env)
  
  validCharacter(value = names.field, name = "field", validLength = NULL, validValues = names(.ls_optionsMRIaggr), refuse.NULL = TRUE, method = "allocOptionsMRIaggr")
  
  for(iter_field in 1:n.args){
    .ls_optionsMRIaggr[[names.field[iter_field]]] <- field[[iter_field]]
  }
  
  validNumeric(value = .ls_optionsMRIaggr$asp, name = "asp", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validLogical(value = .ls_optionsMRIaggr$axes, name = "axes", validLength = 1, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$bg, name = "bg", validLength = 1, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$breaks, name = "breaks", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$cex, name = "cex", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$cex.index, name = "cex.index", validLength = 3, min = 0, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$cex.legend, name = "cex.legend", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$cex.main, name = "cex.main", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validLogical(value = .ls_optionsMRIaggr$checkArguments, name = "checkArguments", validLength = 1, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$col.index, name = "col.index", validLength = 3, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$col.NA, name = "col.NA", validLength = 1, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$col.midplane, name = "col.midplane", validLength = 1, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$digit.legend, name = "digit.legend", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$digit.result, name = "digit.result", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$digit.epsilon, name = "digit.epsilon", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$digit.percentage, name = "digit.percentage", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$filter.index, name = "filter.index", validLength = 1, validValues = c("2D_N4", "2D_N8", "3D_N4", "3D_N6", "3D_N8", "3D_N10", "3D_N18", "3D_N26"), refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$height, name = "height", validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$hemisphere, name = "hemisphere", validLength = 1, validValues = c("both", "left", "right", "lesion", "contralateral"), method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$legend, name = "legend", validLength = 1, validValues = c(TRUE, FALSE, "only"), refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$main, name = "main", validLength = NULL, refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$main.legend, name = "main.legend", validLength = NULL, refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$mar, name = "mar", validLength = 4, min = 0, refuse.NULL = TRUE, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$mar.legend, name = "mar.legend", validLength = 4, min = 0, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$mfrow, name = "mfrow", validLength = 2, min = 0, refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$mgp, name = "mgp", validLength = 3, min = 0, refuse.NULL = TRUE, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$norm_mu, name = "norm_mu", validLength = 1, 
                 validValues = c(FALSE, "global", "global_1slice", "global_3slices", "contralateral", "contralateral_1slice", "contralateral_3slices", "default_value"),
                 method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$norm_sigma, name = "norm_sigma", validLength = 1, 
                 validValues = c(FALSE, "global", "global_1slice", "global_3slices", "contralateral", "contralateral_1slice", "contralateral_3slices", "default_value"),
                 method = "allocOptionsMRIaggr")
  validLogical(value = .ls_optionsMRIaggr$num.main, name = "num.main", validLength = 1, method = "allocOptionsMRIaggr")
  validLogical(value = .ls_optionsMRIaggr$numeric2logical, name = "numeric2logical",  validLength = 1, method = "allocOptionsMRIaggr")
  validLogical(value = .ls_optionsMRIaggr$outline.index, name = "outline.index", validLength = 1, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$palette, name = "palette", validLength = 1, validValues = c("rgb", "hsv", "rainbow", "grey.colors", "heat.colors", "terrain.colors", "topo.colors", "cm.colors"), method = "allocOptionsMRIaggr")
  validPath(value = .ls_optionsMRIaggr$path, name = "path", method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$pch.index, name = "pch.index", validLength = 3, min = 0, method = "allocOptionsMRIaggr")
  validInteger(value = .ls_optionsMRIaggr$pch.NA, name = "pch.NA", validLength = 1, min = 0, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$pty, name = "pty", validLength = 1, validValues = c("m", "s"), refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  validLogical(value = .ls_optionsMRIaggr$quantiles.legend, name = "quantiles.legend", validLength = 1, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$res, name = "res", validLength = 1, refuse.NA = FALSE, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$slice_var, name = "slice_var", validLength = 1, validValues = c("i", "j", "k"), method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$type.breaks, name = "type.breaks", validLength = 1, validValues = c("range", "range_center", "quantile"), method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$unit, name = "unit", validLength = 1, validValues = c("px", "in", "cm", "mm"), refuse.NULL = FALSE, method = "allocOptionsMRIaggr")  
  validLogical(value = .ls_optionsMRIaggr$verbose, name = "verbose",  validLength = 1, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$xlab, name = "xlab", validLength = 1, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$ylab, name = "ylab", validLength = 1, method = "allocOptionsMRIaggr")
  validNumeric(value = .ls_optionsMRIaggr$width, name = "width", validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = "allocOptionsMRIaggr")
  validCharacter(value = .ls_optionsMRIaggr$window, name = "window", validLength = 1, validValues = c(TRUE, FALSE, "eps", "pdf", "png", "svg"), refuse.NULL = FALSE, method = "allocOptionsMRIaggr")
  
  assign(".ls_optionsMRIaggr", .ls_optionsMRIaggr, envir = MRIaggr.env)
  
  return(field)
}
