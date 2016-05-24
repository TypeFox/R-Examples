### Weka filters
##
## Note that a lot of the subclasses of Filter implement basic data 
## manipulation, such as removing an attribute.  As these may not be of
## direct interest to R users, they have not been registered.
##
## ceeboo 2006

make_Weka_filter <-
function(name, class = NULL, init = NULL)
{
    classes <- c(class, "data.frame")
    kind <- "R_Weka_filter_interface"
    name <- as_JNI_name(name)
    meta <- make_R_Weka_interface_metadata(name, kind, classes, init)
    Weka_interfaces[[Java_class_base_name(name)]] <- meta    

    out <- function(formula, data, subset, na.action, control = NULL) {
        mc <- match.call()
        mf <- mc[c(1L, match(c("formula", "data", "subset", "na.action"),
                             names(mc), 0L))]
        ## Need 'stats::' for non-standard evaluation:
        mf[[1L]] <- quote(stats::model.frame)
        mf <- eval(mf, parent.frame())
        
        RWeka_use_filter(mf, control, name, init)
    }

    make_R_Weka_interface(out, meta)
}

RWeka_use_filter <-
function(mf, control, name, init)
{
    if(is.function(init)) init()
    
    ## We do not always need a response variable, e.g., for the class 
    ## of unsupervised filters.

    ## <NOTE>
    ## We do not check the Weka model class.  Thus, the formula may not
    ## fit the model class.  In this case Weka throws, but the rJava
    ## calls do not stop :-(
    ## </NOTE>
    
    if (attr(attr(mf, "terms"), "response") == 0)
       instances <- read_data_into_Weka(mf)
    else
       instances <- read_model_frame_into_Weka(mf)

    ## Build filter.
    filter <- .jnew(name)
    control <- as.character(control)
    if (length(control))
       .jcall(filter, "V", "setOptions", .jarray(control))
    ## This is strange ...
    .jcall(filter, "Z", "setInputFormat", instances)
                       
    ## Unlike Classifier Filter provides a class method for a data set,
    ## i.e., for class Instances :-) 

    instances <- .jcall("weka/filters/Filter", "Lweka/core/Instances;",
                        "useFilter",
                        .jcast(instances,"weka/core/Instances"), 
                        .jcast(filter,"weka/filters/Filter"))
   
    read_instances_from_Weka(instances)
}

## Has no toString method.

print.Weka_filter <-
function(x, ...)
{
    writeLines(.jcall(x$filter, "S", "globalInfo"))
    invisible(x)
}

