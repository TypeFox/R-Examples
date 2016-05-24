.has_method <-
function(o, name)
{
    ## Essentially the same code as in .jmethods().
    
    ## Determine whether Java object o has method 'name'.
    cl <- .jcall(o, "Ljava/lang/Class;", "getClass")
    ms <- .jcall(cl, "[Ljava/lang/reflect/Method;", "getMethods")
    ss <- unlist(lapply(ms, function(x) .jcall(x, "S", "toString")))
    length(grep(paste("\\.", name, "\\(", sep = ''), ss)) > 0L
}

.has_Java_method <- 
function(object, name) 
{
    object <- .jcall(object, "Ljava/lang/Class;", "getClass")
    object <- .jcall(object, "[Ljava/lang/reflect/Method;", "getMethods")
    object <- sapply(object, .jcall, "S", "getName")
    match(name, object, nomatch = 0L) > 0L
}

make_R_Weka_interface <-
function(f, meta)
    .structure(f,
               class = unique(c(meta$kind, "R_Weka_interface")),
               meta = meta)

make_R_Weka_interface_metadata <-
function(name, kind, class = NULL, init = NULL)
    list(name = name, kind = kind, class = class, init = init)

as_JNI_name <-
function(x)
    gsub("\\.", "/", x)

as_qualified_name <-
function(x)
    gsub("/", ".", x)

Java_class_base_name <-
function(x)
    sub(".*[/.]", "", x)

get_Java_class <-
function(x, packages = NULL)
{
    ## For consistency (and simplicity), return qualified names.
    
    .find_Java_class_in_packages <- function(x, packages) {
        classes <- paste(as_JNI_name(packages), x, sep = "/")
        for(cl in classes)
            if(!is.null(.jfindClass(cl, silent = TRUE))) return(cl)
        NULL
    }

    cls <- if(is.character(x)) {
        if(regexpr("[/.]", x) > -1L) {
            ## If possibly a full Java class name, leave alone.
            x
        }
        else {
            ## Otherwise, try treating as the base class name of a Weka
            ## class interfaced and registered.
            cls <- Weka_interfaces[[x]]$name
            ## And finally, try to find within the given packages ...
            if(is.null(cls) && !is.null(packages))
                cls <- .find_Java_class_in_packages(x, packages)
            ## (Shouldn't we do something if we only "find" NULL?
            cls
        }
    }
    else if(inherits(x, "R_Weka_interface"))
        attr(x, "meta")$name
    else
        NULL

    ## Canonicalize.
    if(!is.null(cls)) cls <- as_qualified_name(cls)

    cls
}

get_R_classes_returned <-
function(x)
{
    if(is.character(x))
        Weka_interfaces[[x]]$class
    else if(inherits(x, "R_Weka_interface"))
        attr(x, "meta")$class
    else
        NULL
}

.compose_and_funcall <-
function(flist, x)
{
    if(is.function(flist))
        flist(x)
    else {
        ## flist should be a list of functions.
        for(i in seq_along(flist))
            if(is.function(f <- flist[[i]]))
                x <- f(x)
        x
    }
}

.structure <-
function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))
