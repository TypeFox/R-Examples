make_Weka_file_saver <-
function(name, handlers = list())
{
    ## Create an interface to a Weka class for saving data.

    ## This is somewhat awful: if we want to allow the specification of
    ## control arguments (and we should so that e.g. we can given the
    ## index of the class attribute where needed), it seems that we
    ## *must* give non-empty input ('-i') and output ('-o') file
    ## arguments as well (this is "inherited" from the setOptions()
    ## method for class AbstractFileSaver).  But then we really have to
    ## get the instances from an ARFF file ...

    ## Of course, the upside is that handling control arguments and in
    ## particular finding out about them using WOW() works "usual way".

    kind <- "R_Weka_file_saver_interface"
    name <- as_JNI_name(name)
    meta <- make_R_Weka_interface_metadata(name, kind)
    Weka_interfaces[[Java_class_base_name(name)]] <- meta

    out <- function(x, file, control = NULL) {
        ## For the time being, all we provide is an interface for saving
        ## to a *file*, given by its path.  More general connection
        ## stuff could be handled by saving to a temporary file first
        ## and then put its contents to the desired connection.
        if(!is.character(file) || (length(file) != 1L))
            stop("Argument 'file' must be a character string.")
        ## Write the data to a temporary ARFF file.
        arfff <- tempfile()
        on.exit(unlink(arfff))
        write.arff(x, arfff)
        ## Call the saver with the ARFF file as input, the given output
        ## file, and the given control arguments.
        saver <- .jnew(name)
        control <- as.character(.compose_and_funcall(handlers$control,
                                                     control))
        .jcall(saver, "V", "setOptions",
               .jarray(c(control, "-i", arfff, "-o", path.expand(file))))
        .jcall(saver, "V", "writeBatch")
        invisible()
    }

    make_R_Weka_interface(out, meta)
}

## <NOTE>
## The following variant does not work in general, for the reasons
## indicated above.  Also, we have not found a way to set the output
## file stem for the C4.5 savers which "works" and does not use '-o'
## (implying '-i', see above).
##
## make_Weka_file_saver <-
## function(name)
## {
##     function(x, file, control = NULL) {
##         if(!is.character(file) || (length(file) != 1L))
##             stop("Argument 'file' must be a character string.")
##         saver <- .jnew(name)
##         .jcall(saver, "V", "setInstances", read_data_into_Weka(x))
##         .jcall(saver, "V", "setFile",
##                .jnew("java/io/File", path.expand(file)))
##         .jcall(saver, "V", "writeBatch")
##         invisible()
##     }
## }
## </NOTE>

make_Weka_file_loader <-
function(name)
{
    ## Create an interface to a Weka class for saving data.

    kind <- "R_Weka_file_loader_interface"
    name <- as_JNI_name(name)    
    meta <- make_R_Weka_interface_metadata(name, kind, "data.frame")
    Weka_interfaces[[Java_class_base_name(name)]] <- meta

    out <- function(file) {
        ## As usual, only files for the time being ...
        if(!is.character(file) || (length(file) != 1L))
            stop("Argument 'file' must be a character string.")
        
        loader <- .jnew(name)
        .jcall(loader, "V", "setSource",
               .jnew("java/io/File", path.expand(file)))
        instances <-
            .jcall(loader, "Lweka/core/Instances;", "getDataSet")
        
        out <- read_instances_from_Weka(instances)
        ## For some formats (e.g., XRFF), we may get additional metadata
        ## such as the class index and attribute or instance weights.
        ## For the time being, return at least the name of the class
        ## attribute, and non-trivial instance and attribute weights.
        attr(out, "Weka_class_attribute_name") <-
            .jcall(.jcall(instances, "Lweka/core/Attribute;",
                          "classAttribute"),
                   "S", "name")
        weights <- .jcall("RWekaInterfaces", "[D",
                          "getAttributeWeights", instances)
        if(any(weights != 1))
            attr(out, "attribute_weights") <- weights
        weights <- .jcall("RWekaInterfaces", "[D",
                          "getInstanceWeights", instances)
        if(any(weights != 1))
            attr(out, "instance_weights") <- weights
        out
    }

    make_R_Weka_interface(out, meta)
}
