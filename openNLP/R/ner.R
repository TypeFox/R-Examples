Maxent_Entity_Annotator <-
function(language = "en", kind = "person", probs = FALSE, model = NULL)
{
    f <- Maxent_Simple_Entity_Detector(language, kind, probs, model)
    description <-
        sprintf("Computes entity annotations using the Apache OpenNLP Maxent name finder employing %s.",
                environment(f)$info)

    Simple_Entity_Annotator(f, list(description = description))
}
    
Maxent_Simple_Entity_Detector <-
function(language = "en", kind = "person", probs = FALSE, model = NULL)    
{
    force(language)
    force(probs)
    force(kind)
    
    info <- if(is.null(model)) {
        ## <NOTE>
        ## Alternatively could have used 'type' instead of 'kind' ...
        kind <- match.arg(kind,
                          ## Should be good enough for now ...
                          c("date", "location", "money", "organization",
                            "percentage", "person", "misc"))
        ## </NOTE>
        package <- sprintf("openNLPmodels.%s", language)
        model <- system.file("models",
                             sprintf("%s-ner-%s.bin", language, kind),
                             package = package)
        if(model == "") {
            msg <-
                paste(gettextf("Could not find model file for language '%s' and kind '%s'.",
                               language, kind),
                      if(system.file(package = package) == "") {
                          gettextf("Please make sure package '%s' is installed,\navailable from <http://datacube.wu.ac.at/>.",
                                   package)
                      } else {
                          gettextf("Apparently, package '%s' is installed\nbut does not provide this model.",
                                   package)
                      },
                      sep = "\n")
            stop(msg)            
        }
        sprintf("the default model for language '%s' and kind '%s'",
                language, kind)
    }
    else
        "a user-defined model"

    ## See
    ## <http://opennlp.apache.org/documentation/1.5.3/manual/opennlp.html#tools.namefind.recognition.api>.

    model <- .jnew("opennlp.tools.namefind.TokenNameFinderModel",
                   .jcast(.jnew("java.io.FileInputStream", model),
                          "java.io.InputStream"))

    ref <- .jnew("opennlp.tools.namefind.NameFinderME", model)

    function(x) {
        y <- .jcall(ref, "[Lopennlp/tools/util/Span;", "find",
                    .jarray(x))
        y <- if(!length(y))
            Annotation()
        else {
            start <- as.integer(sapply(y, .jcall, "I", "getStart")) + 1L
            end <- as.integer(sapply(y, .jcall, "I", "getEnd"))
            kind <- as.character(sapply(y, .jcall, "S", "getType"))
            type <- rep.int("entity", length(start))
            features <- lapply(kind, single_feature, "kind")
            if(probs) {
                ## Apparently need the probabilities for the obtained
                ## spans, see
                ## <http://opennlp.apache.org/documentation/1.5.3/apidocs/opennlp-tools/index.html>.
                probs <- .jcall(ref, "[D", "probs",
                                .jcast(.jarray(y),
                                       "[Lopennlp/tools/util/Span;"))
                features <- Map(c,
                                features,
                                lapply(probs, single_feature, "prob"))
            }
            Annotation(NULL, type, start, end, features)                
        }
        .jcall(ref, "V", "clearAdaptiveData")
        y
    }

}
