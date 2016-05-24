Parse_Annotator <-
function()
{
    ## Currently, only provide the chunking parser with the default
    ## model.
    ## Add type = "chunking" and model = NULL arguments evetually.

    package <- "openNLPmodels.en"
    model <- system.file("models/en-parser-chunking.bin",
                         package = package)
    if(model == "")
        stop(gettextf("Could not find model file.\nPlease make sure package '%s' is installed,\navailable from <http://datacube.wu.ac.at/>.",
                      package),
             domain = NA)
    model <- .jnew("opennlp.tools.parser.ParserModel",
                   .jcast(.jnew("java.io.FileInputStream", model),
                          "java.io.InputStream"))
    ## Need to instantiate via ParserFactory, see
    ## <http://opennlp.apache.org/documentation/1.5.3/manual/opennlp.html#tools.parser>.
    ref <- .jcall(.jnew("opennlp.tools.parser.ParserFactory"),
                  "Lopennlp/tools/parser/Parser;",
                  "create",
                  model)

    f <- function(s, a) {

        parse_sentence <- function(as, at) {
            ## Function to parse a single sentence based on the
            ## annotations for the sentence and its word tokens.
            ## The "tricky" part is generating a parse tree from the
            ## sentence and word token spans first, see e.g.
            ## See
            ## <http://blog.dpdearing.com/2011/12/how-to-use-the-opennlp-1-5-0-parser/>.
            ## Note that it should be more efficient to do these
            ## computations directly in Java.
            t_inc <- .jfield("opennlp.tools.parser.AbstractBottomUpParser",
                             name = "INC_NODE")
            t_tok <- .jfield("opennlp.tools.parser.AbstractBottomUpParser",
                             name = "TOK_NODE")
            p <- .jnew("opennlp.tools.parser.Parse",
                       s,
                       .jnew("opennlp.tools.util.Span",
                             as$start - 1L,
                             as$end),
                       t_inc,
                       1,
                       0L)
            for(i in seq_along(at)) {
                p$insert(.jnew("opennlp.tools.parser.Parse",
                               s,
                               .jnew("opennlp.tools.util.Span",
                                     at$start[i] - 1L,
                                     at$end[i]),
                               t_tok,
                               0,
                               as.integer(i) - 1L))
            }
            p <- ref$parse(p)
            sb <- .jnew("java.lang.StringBuffer")
            p$show(sb)
            sb$toString()
        }

        as <- a[a$type == "sentence"]
        at <- annotations_in_spans(a[a$type == "word"], as)

        parses <- sapply(seq_along(at),
                         function(i) parse_sentence(as[i], at[[i]]))

        as$features <- lapply(parses, single_feature, "parse")

        as
    }

    description <-
        "Computes Treebank parse annotations using the Apache OpenNLP chunking parser for English."

    Annotator(f, list(description = description))
}
