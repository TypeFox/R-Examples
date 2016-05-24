activeDataSetDesignccd <- function (){
  aus <- FALSE
  if (activeDataSetDesignP()){
     di <- design.info(eval(parse(text=ActiveDataSet())))
     if (length(grep("ccd",di$type))>0)
        aus <- TRUE
     }
  aus
}

activeDataSetDesignlhs <- function (){
  aus <- FALSE
  if (activeDataSetDesignP()){
     di <- design.info(eval(parse(text=ActiveDataSet())))
     if (length(grep("lhs",di$type))>0)
        aus <- TRUE
     }
  aus
}


activeDataSetDesignparam <- function (){
  aus <- FALSE
  if (activeDataSetDesignP()){
     di <- design.info(eval(parse(text=ActiveDataSet())))
     aus <- length(grep("param",di$type))>0
     }
  aus
}

activeDataSetDesignparamlong <- function ()
  activeDataSetDesignparam() & activeDataSetDesignLongToWide()

activeDataSetDesignLongToWide <- function (){
  aus <- FALSE
  if (activeDataSetDesignP()){
     di <- design.info(eval(parse(text=ActiveDataSet())))
     aus <- di$repeat.only
     if (length(grep("param",di$type))>0 & is.null(di$responselist))
        aus <- TRUE
     }
  aus
}

activeDataSetDesignWide <- function (){
  aus <- FALSE
  if (activeDataSetDesignPResp()){
     di <- design.info(eval(parse(text=ActiveDataSet())))
     aus <- !is.null(di$responselist)
     }
  aus
}

activeDataSetDesignP <- function (){
  aus <- FALSE
  if (activeDataSetP())
     aus <- "design" %in% class(eval(parse(text=ActiveDataSet())))
  aus
}

activeDataSetDesignPResp <- function (){
  aus <- FALSE
  if (activeDataSetDesignP())
     aus <- !is.null(response.names(eval(parse(text=ActiveDataSet()))))
  aus
}


activeDataSetDesign2P <- function (){
  aus <- FALSE
  if (activeDataSetDesignP())
     aus <- isDesign2pb(eval(parse(text=ActiveDataSet()))) | isDesign2FrF(eval(parse(text=ActiveDataSet())))
  aus
}

activeDataSetDesign2Pwoc <- function (){
  aus <- FALSE
  if (activeDataSetDesignP())
     aus <- is.null(design.info(eval(parse(text=ActiveDataSet())))$ncenter)
  aus
}

activeDataSetDesign2Pnoccd <- function (){
  aus <- FALSE
  if (activeDataSetDesign2P())
     aus <- !activeDataSetDesignccd()
  aus
}

activeDataSetDesign2PResp <- function (){
  aus <- FALSE
  if (activeDataSetDesign2P())
     aus <- !is.null(response.names(eval(parse(text=ActiveDataSet()))))
  aus
}

activeModelRSM <- function(){
  aus <- FALSE
  if (activeModelP())
     if ("rsm" %in% class(get(.activeModel))) aus <- TRUE
  aus
}

activeModelLM <- function(){
  aus <- FALSE
  if (activeModelP())
     if ("lm" %in% class(get(.activeModel)) & 
         !(any(c("glm","mlm") %in% class(get(.activeModel))))) aus <- TRUE
  aus
}

exist2Designs <- function(){
   length(listDesigns()) > 1
}

existDesigns <- function(){
   length(listDesigns()) > 0
}

existDesign2 <- function(){
   length(listDesigns2()) > 0
}

existDesign2pb <- function(){
   length(listDesigns2(type=="pb")) > 0
}
existDesign2FrF <- function(){
   length(listDesigns2(type=="FrF2")) > 0
}

existDesignsWithResp <- function(){
    length(listDesignsWithResp()) > 0
}

existcatlgPkgs <- function(){
    "FrF2.catlg128" %in% .packages()
    ## das funktioniert nicht
    }

existRSMs <- function()
    length(listRSMs())>=1


isDesign2pb <- function(design){
        aus <- FALSE
        if (design.info(design)$type=="pb"){ 
           aus <- TRUE 
           return(aus)}
        if (design.info(design)$type=="oa"){ 
           nlevels <- design.info(design)$nlevels
           if (all(nlevels==2)) aus <- TRUE}
        aus
}
isDesign2FrF <- function(design){
        aus <- FALSE
        if (substr(design.info(design)$type,1,4) == "FrF2"){ 
           aus <- TRUE 
           return(aus)}
        if (substr(design.info(design)$type,1,14) == "full factorial"){ 
           nlevels <- design.info(design)$nlevels
           if (all(nlevels==2)) aus <- TRUE}
        aus
}

activeDataSetDesignQualP <- function (){
  aus <- FALSE
  if (activeDataSetDesignP()){
     hilf <- design.info(eval(parse(text=ActiveDataSet())))$quantitative
     aus <- !all(hilf) | is.null(hilf)
     }
  aus
}

activeDataSetDesignRemovableP <- function (){
  aus <- FALSE
  if (activeDataSetDesignP()){
     di <- design.info(eval(parse(text=ActiveDataSet())))
     if (length(setdiff(colnames(eval(parse(text=ActiveDataSet()))), c(names(di$factor.names), di$block.name)))>0)
     aus <- TRUE
     }
  aus
}

