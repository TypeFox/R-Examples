MissingNS = c(gating = "http://www.crap.org",
              'data-type' = "http://www.morecrap.org")

fixXMLNamespaces =
  #
  #  call as
  #    dd = fixXMLNamespaces("~/v75_step6.wsp", .namespaces = MissingNS)
  #  or
  #   dd = fixXMLNamespaces("~/v75_step6.wsp", gating = "http://www.crap.org", 'data-type' = "http://www.morecrap.org")
  #
function(doc = "~/v75_step6.wsp", ..., .namespaces = list(...)) 
{
    # collect the error messages
  e = xmlErrorCumulator(, FALSE)
  doc = xmlParse(doc, error = e)

  if(length(e) == 0)
     return(doc)

     # find the ones that refer to prefixes that are not defined
  ns = grep("^Namespace prefix .* not defined", unique(environment(e)$messages), value = TRUE)
  ns = unique(gsub("Namespace prefix ([^ ]+) .*", "\\1", ns))

    # now set those name spaces on the root of the document
  if(is(.namespaces, "list"))
    .namespaces = structure(as.character(unlist(.namespaces)), names = names(.namespaces))

  uris = .namespaces[ns]
  if(length(uris)) {
     mapply(function(id, uri)
              newXMLNamespace(xmlRoot(doc), uri, id),
            names(uris), uris)
     xmlParse(saveXML(doc), asText = TRUE)
  } else
     doc
}
