make_Weka_stemmer <-
function(name)
{
    ## Create an interface to a Weka class for stemming.

    kind <- "R_Weka_stemmer_interface"
    name <- as_JNI_name(name)
    meta <- make_R_Weka_interface_metadata(name, kind, "character")
    Weka_interfaces[[Java_class_base_name(name)]] <- meta
    
    out <- function(x, control = NULL) {
        stemmer <- .jnew(name)
        control <- as.character(control)
        if(length(control)) {
            if(.has_method(stemmer, "setOptions"))
                .jcall(stemmer, "V", "setOptions", .jarray(control))
            else
                warning("Stemmer cannot set control options.")
        }
        .jcall("RWekaInterfaces", "[S", "stem",
               .jcast(stemmer, "weka/core/stemmers/Stemmer"),
               .jarray(as.character(x)))
    }

    make_R_Weka_interface(out, meta)
}

