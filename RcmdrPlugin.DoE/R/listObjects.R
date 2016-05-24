listObjects <- function (envir = .GlobalEnv, ...)
{
    ls(envir = envir, all.names = TRUE)
}

listLists <- function (envir = .GlobalEnv, ...)
{
   Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    names(which(sapply(Vars, function(.x) is.list(get(.x,
        envir = envir)))))
}

listCatlgs <- function (envir = .GlobalEnv, ...)
{
   Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    c("catlg", names(which(sapply(Vars, function(.x) "catlg" %in% class(get(.x,
        envir = envir))))))
}

listDesigns <- function (envir = .GlobalEnv, ...)
{
   Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    names(which(sapply(Vars, function(.x) "design" %in% class(get(.x,
        envir = envir)))))
}

listDesignsWithResp <- function (envir = .GlobalEnv, ...)
{
   designs <- listDesigns()
   if (length(designs)>0){
     mitresp <- sapply(designs, function(obj) !is.null(response.names(eval(parse(text=obj)))))
     designs <- designs[mitresp]}
   designs
}


listDesigns2 <- function (envir = .GlobalEnv, type = NULL, ...)
{
   Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    hilf <- names(which(sapply(Vars, function(.x) "design" %in% class(get(.x,
        envir = envir)))))
    if (length(hilf) == 0)
        return(hilf)
    if (is.null(type)) 
      aus <- hilf[which(sapply(hilf, function(.x) (isDesign2pb(get(.x,
        envir = envir)) | isDesign2FrF(get(.x,
        envir = envir)))))]
    else if (type=="FrF2") 
      aus <- hilf[which(sapply(hilf, function(.x) isDesign2FrF(get(.x,
        envir = envir))))]
    else if (type=="pb") 
      aus <- hilf[which(sapply(hilf, function(.x) isDesign2pb(get(.x,
        envir = envir))))]
    aus
}

listRSMs <- function (envir = .GlobalEnv, ...)
{
   Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
   names(which(sapply(Vars, function(.x) "rsm" %in% class(get(.x,
        envir = envir)))))
}

