Maxent_Sent_Token_Annotator <-
function(language = "en", probs = FALSE, model = NULL)
{
    f <- Maxent_Simple_Sent_Tokenizer(language, probs, model)
    description <-
        sprintf("Computes sentence annotations using the Apache OpenNLP Maxent sentence detector employing %s.",
                environment(f)$info)
    
    Simple_Sent_Token_Annotator(f, list(description = description))
}

Maxent_Simple_Sent_Tokenizer <-
function(language = "en", probs = FALSE, model = NULL)
{
    force(language)
    force(probs)
    
    info <- if(is.null(model)) {
        package <- if(language == "en")
            "openNLPdata"
        else
            sprintf("openNLPmodels.%s", language)
        model <- system.file("models",
                             sprintf("%s-sent.bin", language),
                             package = package)
        if(model == "") {
            msg <-
                paste(gettextf("Could not find model file for language '%s'.",
                               language),
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
        sprintf("the default model for language '%s'", language)
    }
    else
        "a user-defined model"

    ## See
    ## <http://opennlp.apache.org/documentation/1.5.3/manual/opennlp.html#tools.sentdetect.detection.api>.

    model <- .jnew("opennlp.tools.sentdetect.SentenceModel",
                   .jcast(.jnew("java.io.FileInputStream", model),
                          "java.io.InputStream"))

    ref <- .jnew("opennlp.tools.sentdetect.SentenceDetectorME", model)

    function(x) {
        y <- .jcall(ref,
                    "[Lopennlp/tools/util/Span;",
                    "sentPosDetect",
                    x)
        start <- as.integer(sapply(y, .jcall, "I", "getStart")) + 1L
        end <- as.integer(sapply(y, .jcall, "I", "getEnd"))
        if(probs) {
            probs <- .jcall(ref, "[D", "getSentenceProbabilities")
            Annotation(NULL,
                       rep.int("sentence", length(start)),
                       start,
                       end,
                       lapply(probs, single_feature, "prob"))
        } else 
            Span(start, end)
    }
    
}
