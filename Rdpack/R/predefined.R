.rd_sections <-      # note: user defined sections are typeset before note, seealso, examples.
    c("\\name", "\\Rdversion", "\\docType", "\\alias", "\\encoding", "\\concept",
      "\\title", "\\description", "\\usage", "\\format", "\\source", "\\arguments",
      "\\details", "\\value", "\\references", "\\section",
      "\\note", "\\author", "\\seealso", "\\examples", "\\keyword")

    ##  todo: "\\Sexpr", "\\Rdopts", "\\newcommand", "\\renewcommand"

rdo_top_tags <- unique(c(.rd_sections,
                         "#ifdef", "#ifndef",
                         "\\newcommand", "\\renewcommand",
                         "COMMENT", "TEXT"
                         ))


Rdo_piece_types <- c(  name = "VERB"             # sections
                     , alias = "VERB"
                     , concept =  "TEXT"
                     , docType =  "TEXT"
                     , title =  "TEXT"
                     , description =  "TEXT"
                     , examples = "RCODE"
                     , usage = "RCODE"
                     , Rdversion = "VERB"
                     , synopsis = "VERB"

                     , Sexpr = "RCODE"
                     , RdOpts = "VERB"

                     , code = "RCODE"            # macros
                     , dontshow = "RCODE"
                     , donttest = "RCODE"
                     , testonly = "RCODE"

                     , dontrun = "VERB"
                     , env  = "VERB"
                     , kbd  = "VERB"
                     , option  = "VERB"
                     , out  = "VERB"
                     , preformatted = "VERB"
                     , samp  = "VERB"
                     , special  = "VERB"
                     , url  = "VERB"
                     , verb = "VERB"
                     , deqn = "VERB"
                     , eqn = "VERB"
                     , renewcommand  = "VERB"
                     , newcommand  = "VERB"
                     )

                                        # todo: not complete. Skaniray instalirana
                                        #                     dokumentatsiya za tezi raboti!
Rdo_predefined_sections <- c(  name = "VERB"
                             , alias = "VERB"
                             , concept =  "TEXT"
                             , docType =  "TEXT"
                             , title =  "TEXT"
                             , description =  "TEXT"
                             , examples = "RCODE"
                             , usage = "RCODE"
                             , Rdversion = "VERB"
                             , synopsis = "VERB"
                             , section = "TEXT"   # not clear what to put for this element
                             # (and the following) 2012-09-23 new entries below.  (maybe
                             # their absence was the reason for some strange behaviour of
                             # char2Rdpiece). Slagam vsichkite "TEXT", no (todo:) tryabva da
                             # se proveri.
                             , arguments = "TEXT"
                             , keyword = "TEXT"
                             , note = "TEXT"

                             , format = "TEXT"
                             , source = "TEXT"
                             , details = "TEXT"
                             , value = "TEXT"
                             , references = "TEXT"
                             , author = "TEXT"
                             , seealso = "TEXT"
                             )
