setMethodS3("listToXml", "list", function(tree, indentStep=" ", collapse="\t", ...) {
##   tree <- list(
##     chipType = "GenomeWideSNP_6",
##     createdOn = "2008-02-13",
##     source = list(
##       url = "http://www.affymetrix.com/",
##       files = ""
##     ),
##     creator = list(
##       name  = "Henrik Bengtsson",
##       email = "hb@stat.berkeley.edu"
##     )
##   )

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asXml <- function(tree, indentStep=" ", .indent="", collapse="\t") {
    xml <- c();

    names <- names(tree);
    for (name in names) {
      child <- tree[[name]];
      if (is.list(child)) {
        xmlChild <- asXml(child, indentStep=indentStep, 
                               .indent=paste(.indent, indentStep, sep=""));
        xmlChild <- sprintf("%s<%s>\n%s%s</%s>\n", 
                     .indent, name, as.character(xmlChild), .indent, name);
      } else {
        child <- paste(child, collapse=collapse);
        xmlChild <- sprintf("%s<%s>%s</%s>\n",
                                 .indent, name, child, name);
      }
      xml <- c(xml, xmlChild);
    }

    paste(xml, collapse="");
  } # asXml()  


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asXml(tree, indentStep=indentStep, collapse=collapse, ...);
}) # listToXml()



setMethodS3("xmlToList", "character", function(xml, ..., drop=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parseXml <- function(xml, beginTag=NA, ...) {
#    cat("############################################\n");
#    cat("<", beginTag, ">\n", sep="");
#    cat("############################################\n");

    xmlTree <- list();
    if (length(xml) == 0)
      return(xmlTree);

    state <- "inBody";
    body <- tag <- "";
    while (length(xml) > 0) {
#      cat(sprintf("%s (%d): %s\n", state, length(xml), paste(xml, collapse="")));
#      cat(sprintf("tag: (%d): %s\n", nchar(tag), tag));
#      cat(sprintf("body: (%d): %s\n", nchar(body), body));
#      cat("\n");

      xmlNext <- xml[1];
      xml <- xml[-1];

      if (state == "inBeginTag") {
        if (xmlNext == ">") {
          res <- parseXml(xml, beginTag=trim(tag), ...);
          xmlNode <- res[[1]];

          # Drop dimension of values
          if (drop && is.null(names(xmlNode))) {
            xmlNode <- xmlNode[[1]];
          }

          xmlNode <- list(xmlNode);
          names(xmlNode) <- trim(tag);

          xmlTree <- c(xmlTree, xmlNode);
          xml <- res$xml;
##          str(xmlNode);
##          cat("############################################\n");
##          cat("</", trim(tag), ">\n", sep="");
##          cat("############################################\n");
          tag <- "";
          state <- "inBody";
        } else {
          tag <- paste(tag, xmlNext, sep="");
        }
      } else if (state == "inEndTag") {
        if (xmlNext == ">") {
          endTag <- trim(tag);
          tag <- "";
          if (endTag != beginTag) {
            throw("End tag does not match expected begin tag: ", 
                                               endTag, " != ", beginTag);
          }
          if (length(xmlTree) == 1) {
            xmlNode <- xmlTree;
          } else {
            xmlNode <- list(xmlTree);
          }
          return(list(xmlNode=xmlNode, xml=xml));
        } else {
          tag <- paste(tag, xmlNext, sep="");
        }
      } else if (state == "inBody") {
        if (xmlNext == "<") {
          if (xml[1] == "/") {
            xml <- xml[-1];
            body <- trim(body);
            if (nchar(body) > 0) {
              xmlNode <- list(body);
              xmlTree <- c(xmlTree, xmlNode);
            }
            body <- "";
            state <- "inEndTag";
          } else {
#            xmlNode <- list(body);
#            xmlTree <- c(xmlTree, list(xmlNode));
            state <- "inBeginTag";
          }
        } else {
          body <- paste(body, xmlNext, sep="");
        }
      }
    } # while(...)

    xmlTree;
  } # parseXml()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  xml <- paste(xml, collapse="");
  xml <- trim(xml);
  xml <- unlist(strsplit(xml, split=""));
  parseXml(xml);
}) # xmlToList()



############################################################################
# HISTORY:
# 2009-01-12
# o Added argument 'drop' to xmlToList().
# 2008-05-23
# o UPDATED: Rewrote internal parseXml() to be a "real" parser.  Hopefully,
#   it is a bit more generic now.
# 2008-03-05
# o BUG FIX: Regular expression pattern 'a-Z' is illegal on (at least) some
#   locale (where 'A-z' works). The only way to get specify the ASCII
#   alphabet is to specify all characters explicitly.
# 2008-02-13
# o Validation test: identical(tree, xmlToList(listToXml(tree)))
# o Created to support read/write of footer in AromaTabularBinaryFile.R.
############################################################################
