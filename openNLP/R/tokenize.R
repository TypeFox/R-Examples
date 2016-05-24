Maxent_Word_Token_Annotator <-
function(language = "en", probs = FALSE, model = NULL)
{
    f <- Maxent_Simple_Word_Tokenizer(language, probs, model)
    description <-
        sprintf("Computes word token annotations using the Apache OpenNLP Maxent tokenizer employing %s.",
                environment(f)$info)

    Simple_Word_Token_Annotator(f, list(description = description))
}

Maxent_Simple_Word_Tokenizer <-    
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
                             sprintf("%s-token.bin", language),
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
    ## <http://opennlp.apache.org/documentation/1.5.3/manual/opennlp.html#tools.tokenizer.api>.

    model <- .jnew("opennlp.tools.tokenize.TokenizerModel",
                   .jcast(.jnew("java.io.FileInputStream", model),
                          "java.io.InputStream"))

    ref <- .jnew("opennlp.tools.tokenize.TokenizerME", model)

    function(x) {
        y <- .jcall(ref, "[Lopennlp/tools/util/Span;", "tokenizePos", x)
        start <- as.integer(sapply(y, .jcall, "I", "getStart")) + 1L
        end <- as.integer(sapply(y, .jcall, "I", "getEnd"))
        if(probs) {
            probs <- .jcall(ref, "[D", "getTokenProbabilities")
            Annotation(NULL,
                       rep.int("word", length(start)),
                       start,
                       end,
                       lapply(probs, single_feature, "prob"))
        } else 
            Span(start, end)
    }

}
