##
## NOTES: To Pandoc
##
## Pandoc - API - Changes
##
## pandoc 1.16 (2016-01-02): 
##     Added Attr field to Link and Image (Mauro Bieg, #261, API change).
## pandoc 1.12 (2013-09-15): 
##     various changes
##
pandoc <- new.env()
pandoc$version <- 1.16

## Imports
#' @importFrom "stats" "setNames"

#  -----------------------------------------------------------
#  Filter
#  ======
#' @title Filter JSON-formatted AST
#' @description Apply a filter on the JSON-formatted abstract syntax tree (AST).
#' @param FUN the function to be applied on the AST.
#' @param ... optional arguments to \code{FUN}.
#' @param input a connection object or a character string from which the JSON-formatted AST is read
#' @param output a connection object or a character string to which the JSON-formatted AST is written
#' @export
#  -----------------------------------------------------------
filter <- function(FUN, ..., input = stdin(), output = stdout()) {
    ## read ast (in json format)
    
    json <- character(0)
    while(length(line <- readLines(input, n=1)) > 0) {
        json <- c(json, as.character(line))
    }
    
    if ( length(json) == 1 ) {
        if ( nchar(json) == 0 ){
            stop("InputError: The JSON-formatted AST for read in is empty!")
        }
    }

    ## convert json to native r
    x <- jsonlite::fromJSON(json, simplifyVector = FALSE, flatten=TRUE)

    ## modify the ast tree
    x <- astrapply(x, FUN, ...)

    out <- as.character( jsonlite::toJSON(x, auto_unbox=TRUE) )
    writeLines(out, con=output)
}

## Walk the ast tree and apply the function FUN to every object.
## Returns a new modified AST.
#  -----------------------------------------------------------
#  astrapply
#  =========
#' @title Apply a Function on a AST
#' @description Apply the function \code{FUN} on the abstract syntax tree (AST) obtained from pandoc.
#' @param x a list representing the AST obtained from pandoc.
#' @param FUN the function to be applied to the AST.
#' @param ... optional arguments to \code{FUN}.
#' @return A list containing the modified AST.
#' @export
#  -----------------------------------------------------------
astrapply <- function(x, FUN, ...) {
    if ( is.list(x) ) {
        if ( is.null(names(x)) ) {
            obj <- list()
            for (item in x) {
                if ( is.list(item) & ("t" %in% names(item)) ) {
                    res <- FUN(item[['t']], item[['c']], ...)
                    if ( is.null(res) ) {
                        obj[[length(obj) + 1]] <- astrapply(item, FUN, ...)
                    } else if ( is.list(res) & is.null(names(res)) ) {
                        for (z in res) {
                            obj[[length(obj) + 1]] <- astrapply(z, FUN, ...)
                        }
                    } else {
                        obj[[length(obj) + 1]] <- astrapply(res, FUN, ...)
                    }
                } else {
                    obj[[length(obj) + 1]] <- astrapply(item, FUN, ...)
                }
            }
            return( obj )
        } else {
            obj <- nlist()
            for (k in names(x)) {
                obj[[k]] <- astrapply(x[[k]], FUN, ...)
            }
            return( obj )
        }
    }
    return( x )
}

#  -----------------------------------------------------------
#  document
#  ========
#' @title Create a new Document
#' @description A constructor of an object of type \code{"document"}.
#' @details Each document has the following methods:
#' \itemize{
#'   \item[] \code{to_json()}
#'      \itemize{ 
#'        \item[] \emph{Description}
#'        \item[] \itemize{ Returns the \code{JSON} representation of the document. } }
#'   \item[] \code{write(con, format="markdown", writer=pandocfilters_writer)}
#'      \itemize{ 
#'        \item[] \emph{Description}
#'        \item[] \itemize{\item[] Write the JSON-formatted AST to a connection.}
#'        \item[] \emph{Arguments}
#'        \item[] \itemize{\item[] \code{con}    \sspace a connection object or a character string to which the document is written  }
#'        \item[] \itemize{\item[] \code{format} \sspace a character string giving the format (e.g. \code{"latex"}, \code{"html"})  }
#'        \item[] \itemize{\item[] \code{writer} \sspace an optional writer function, see \link{pandocfilters_writer} }    
#'        \item[] \emph{Note}
#'        \item[] \itemize{\item[] Any function with the three arguments \code{x}, \code{con} and \code{format} can be used as writer function.}
#'      }
#'   \item[] \code{append(x)}
#'      \itemize{ 
#'        \item[] \emph{Description}
#'        \item[] \itemize{\item[] Append a new block to the document.}
#'        \item[] \emph{Arguments}
#'        \item[] \itemize{\item[] \code{x} \sspace a block object or list of block objects} }
#'   \item[] \code{append_plain(x)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{Plain}. }
#'   \item[] \code{append_para(x)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{Para}.}
#'   \item[] \code{append_code_block(attr, code)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{CodeBlock}.}
#'   \item[] \code{append_block_quote(blocks)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{BlockQuote}.}
#'   \item[] \code{append_ordered_list(lattr, lblocks)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{OrderedList}.}
#'   \item[] \code{append_bullet_list(lblocks)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{BulletList}.}
#'   \item[] \code{append_definition_list(x)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{DefinitionList}.}
#'   \item[] \code{append_header(x, level=1L, attr=Attr())} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{Header}.}
#'   \item[] \code{append_horizontal_rule()} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{HorizontalRule}.}
#'   \item[] \code{append_table(rows, col_names=NULL, aligns=NULL, col_width=NULL, caption=list())} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{Table}.}
#'   \item[] \code{append_div(blocks, attr)} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{Div}.}
#'   \item[] \code{append_null()} \sspace
#'      \itemize{ \item[] For more information about the arguments see \link{Null}.}
#' }
#' @export
#  -----------------------------------------------------------
document <- function() {
    env <- new.env()
    env$doc <- list()
    env$meta <- nlist()
    env$append <- function(x) {
        self <- parent.env(environment())$env
        self$doc <- c(self$doc, as.lobo(x))
        invisible(NULL)
    }
    ##  1. Plain
    env$append_plain <- function(x) {
        self <- parent.env(environment())$env
        self$append(Plain(x))
    }
    ##  2. Para
    env$append_para <- function(x) {
        self <- parent.env(environment())$env
        self$append(Para(x))
    }
    ##  3. CodeBlock
    env$append_code_block <- function(attr, code) {
        self <- parent.env(environment())$env
        self$append(CodeBlock(attr, code))
    }
    ##  4. RawBlock
    ##  5. BlockQuote
    env$append_block_quote <- function(blocks) {
        self <- parent.env(environment())$env
        self$append(BlockQuote(blocks))
    }
    ##  6. OrderedList
    env$append_ordered_list <- function(lattr, lblocks) {
        self <- parent.env(environment())$env
        self$append(OrderedList(lattr, lblocks))
    }
    ##  7. BulletList
    env$append_bullet_list <- function(lblocks) {
        self <- parent.env(environment())$env
        self$append(BulletList(lblocks))
    }
    ##  8. DefinitionList
    env$append_definition_list <- function(x) {
        self <- parent.env(environment())$env
        self$append(DefinitionList(x))
    }
    ##  9. Header
    env$append_header <- function(x, level=1L, attr=Attr()) {
        self <- parent.env(environment())$env
        self$append(Header(x, level, attr))
    }
    ## 10. HorizontalRule 
    env$append_horizontal_rule <- function() {
        self <- parent.env(environment())$env
        self$append(HorizontalRule())
    }
    ## 11. Table
    env$append_table <- function(rows, col_names=NULL, aligns=NULL, col_width=NULL, caption=list()) {
        self <- parent.env(environment())$env
        self$append(Table(rows, col_names, aligns, col_width, caption))
    }
    ## 12. Div
    env$append_div <- function(blocks, attr) {
        self <- parent.env(environment())$env
        self$append(Div(blocks, attr))
    }
    ## 13. Null
    env$append_null <- function() {
        self <- parent.env(environment())$env
        self$append(Null())
    }
    ## to_json
    env$to_json <- function() {
        self <- parent.env(environment())$env
        d <- list(list(unMeta=self$meta), self$doc)
        return( jsonlite::toJSON(d, auto_unbox=TRUE) )
    }
    ## Write
    env$write <- function(con, format="markdown", writer=pandocfilters_writer) {
        self <- parent.env(environment())$env
        writer(self$to_json(), con, format)
    }
    structure(env, class="document")
}

##' @noRd
##' @export
print.document <- function(x, ...) print("A pandoc document.")

## ----------------------------------------------------------------------------- 
##
##   Additional Constructors (only used as function arguments)
##
## -----------------------------------------------------------------------------

## -----------------------------------------------------------
##  Inline Objects
##  ==============
#' @title Inline Objects
#' @description 
#'    Objects of the classes \code{"NULL"} and \code{"character"} 
#'    can be coerced to \code{"inline"}.
#' @param x an object of type \code{"NULL"}, \code{"character"} or \code{"inline"}.
#' @return an object of class \code{"inline"}.
#' @examples
#' as.inline("some text")
#' as.inline(NULL)
#' @export
as.inline <- function( x ) {
    UseMethod( "as.inline" )
}

##' @noRd
##' @export
as.inline.inline <- identity

##' @noRd
##' @export
as.inline.character <- function( x ) {
    Str(paste(x, collapse=" "))
}

##' @noRd
##' @export
as.inline.NULL <- function( x ) structure(list(), class=c("inline", "list"))

#' @title Inline Objects
#' @description 
#'   Tests if an object has the class attribute \code{"inline"}.
#' @param x an object to be tested.
#' @return a logical indicating if the provided object is of type \code{"inline"}.
#' @examples
#' is.inline(as.inline(NULL))
#' @export
is.inline <- function(x) class(x)[1] == "inline"

combine_two <- function(x, y) {
    if ( is.null(x) ) return(y)
    if ( is.null(y) ) return(x)
    if ( is.inline(x) & is.inline(y)) {
        return(list(x, y))
    }
    if ( is.loio(x) & is.inline(y) ) {
        return(c(x, list(y)))
    }
    if ( is.inline(x) & is.loio(y) ) {
        return(c(list(x), y))
    }
    return( c(as.loio(x), as.loio(y)) )
}

#' @title Combine Inline Objects
#' @description 
#'    Objects of class \code{"inline"} can be combined by using the generic 
#'    default method \code{"c"} (combine).
#' @param ... objects to be concatenated.
#' @return an list of \code{"inline"} objects.
#' @examples
#' c(Str("some"), Strong("text"))
#' @export
c.inline <- function(...) {
    x <- lapply(list(...), as.inline)
    rval <- Reduce(combine_two, x)
    if ( length(rval) == 1 ) return(structure(rval, class=c("inline", "list")))
    return(structure(rval, class=c("loio", "list")))
}

## -----------------------------------------------------------
##  loio (List of Inline Objects)
##  =============================

as.loio <- function( x ) UseMethod( "as.loio" )

##' @noRd
##' @S3method as.loio loio
as.loio.loio <- identity

##' @noRd
##' @S3method as.loio NULL
as.loio.NULL <- function( x ) structure(list(), class=c("loio", "list"))

##' @noRd
##' @S3method as.loio inline
as.loio.inline <- function( x ) structure(list(x), class=c("loio", "list"))

##' @noRd
##' @S3method as.loio character
as.loio.character <- function( x ) structure(list(as.inline(x)), class=c("loio", "list"))

##' @noRd
##' @S3method as.loio list
as.loio.list <- function( x ) {
    x <- lapply(x, as.inline)
    structure(x, class=c("loio", "list"))
}

is.loio <- function(x) class(x)[1] == "loio"

## -----------------------------------------------------------
##  Block Objects
##  =============
#' @title Block Objects
#' @description 
#'   In pandoc \code{"block"} objects are used as container for 
#'   \code{"inline"} objects and to give them specific roles. 
#'   Objects of the classes \code{"NULL"} and \code{"character"} 
#'   can be coerced to \code{"block"}.
#' @param x an object of type \code{"NULL"} or \code{"character"} or \code{"block"}.
#' @return an object of class \code{"block"}.
#' @examples
#' as.block("some text")
#' as.block(NULL)
#' @export
as.block <- function( x ) {
    UseMethod( "as.block" )
}

##' @noRd
##' @export
as.block.NULL <- function(x) structure(list(), class=c("block", "list"))

##' @noRd
##' @export
as.block.character <- function(x) Plain(x)

#' @title Block Objects
#' @description 
#'   Tests if an object has the class attribute \code{"block"}.
#' @param x an object to be tested.
#' @return a logical indicating if the provided object is of type \code{"block"}.
#' @examples
#' is.block(as.block(NULL))
#' @export
is.block <- function(x) class(x)[1] == "block"

combine_two_blocks <- function(x, y) {
    if ( is.null(x) ) return(y)
    if ( is.null(y) ) return(x)
    if ( is.block(x) & is.block(y)) {
        return(list(x, y))
    }
    if ( is.lobo(x) & is.block(y) ) {
        return(c(x, list(y)))
    }
    if ( is.block(x) & is.lobo(y) ) {
        return(c(list(x), y))
    }
    return( c(as.lobo(x), as.lobo(y)) )
}

#' @title Combine Block Objects
#' @description 
#'    Objects of class \code{"block"} can be combined by using the generic 
#'    default method \code{"c"} (combine).
#' @param ... objects to be concatenated.
#' @return an list of \code{"block"} objects.
#' @examples
#' c(Header( "R Basics" ), Header("What is R?", level=2),
#' Plain(c(Emph("R"), Space(), "is a system for ", Strong("statistical computation"))))
#' @export
c.block <- function(...) {
    x <- list(...)
    rval <- Reduce(combine_two_blocks, x)
    if ( length(rval) == 1 ) return(structure(rval, class=c("block", "list")))
    return( structure(rval, class=c("lobo", "list")) )
}

## -----------------------------------------------------------
##  lobo (List of Block Objects)
##  ============================
as.lobo <- function( x ) UseMethod( "as.lobo" )

##' @noRd
##' @S3method as.lobo lobo
as.lobo.lobo <- identity

##' @noRd
##' @S3method as.lobo NULL
as.lobo.NULL <- function( x ) structure(list(), class=c("lobo", "list"))

##' @noRd
##' @S3method as.lobo block
as.lobo.block <- function( x ) structure(list(x), class=c("lobo", "list"))

##' @noRd
##' @S3method as.lobo list
as.lobo.list <- function( x ) {
    b <- sapply(x, is.block)
    if ( !all(b) ) {
        stop(sprintf("TypeError: elements %s are not of type block ", 
             paste(which(!b), collapse=", ")), "All elements must be of type block!")
    }
    class(x) <- c("lobo", class(x))
    x
}

is.lobo <- function(x) class(x)[1] == "lobo"

## -----------------------------------------------------------
## lolobo (List of List of Block Objects)
## ======================================
as.lolobo <- function( x ) UseMethod( "as.lolobo" )

##' @noRd
##' @S3method as.lolobo lolobo
as.lolobo.lolobo <- identity

##' @noRd
##' @S3method as.lolobo NULL
as.lolobo.NULL  <- function( x ) structure(list(as.lobo(x)), class=c("lolobo", "lobo", "list"))

##' @noRd
##' @S3method as.lolobo block
as.lolobo.block <- function( x ) structure(list(as.lobo(x)), class=c("lolobo", "lobo", "list"))

##' @noRd
##' @S3method as.lolobo lobo
as.lolobo.lobo  <- function( x ) structure(list(x), class=c("lolobo", "lobo", "list"))

##' @noRd
##' @S3method as.lolobo list
as.lolobo.list  <- function( x ) {
    structure(lapply(x, as.lobo), class=c("lolobo", "lobo", "list"))
}

#  -----------------------------------------------------------
#  Citation
#  ========
#' @title Citation
#' @description A constructor of an object of type \code{"Citation"}.
#' @param suffix a inline object or list of inline objects
#' @param id a character string (not visible in the text)
#' @param note_num an integer 
#' @param mode a character string giving the citation mode, possible values are 
#'             \code{"AuthorInText"}, \code{"SuppressAuthor"} and \code{"NormalCitation"}.
#' @param prefix a inline object or list of inline objects
#' @param hash an integer
#' @export
#  -----------------------------------------------------------
Citation <- function(suffix, id, note_num=0L, mode="AuthorInText", prefix=list(), hash=0L) {
    suffix <- as.loio(suffix)
    prefix <- as.loio(prefix)
    x <- list(citationSuffix = suffix,
              citationNoteNum = note_num,
              citationMode = list(t=mode, c=list()),
              citationPrefix = list(),
              citationId = id,
              citationHash = hash)
    structure(x, class=c("Citation", "list"))
}

is.citation <- function(x) class(x)[1] == "Citation"

## -----------------------------------------------------------
## Type
## ====
## title pandoc Types
## description A constructor for pandoc types.
## param x a character string giving the type
## details A convenience function to create the following data structure
##          \code{list(t=x, c=list())} by only providing x.
## examples
## Type("SmallCaps")
## -----------------------------------------------------------
Type <- function(x) structure(list(t=x, c=list()), class=c("Type", "list"))

#  -----------------------------------------------------------
#  Attr
#  ====
#' @title Attributes
#' @description A constructor for pandoc attributes.
#' @param identifier a character string
#' @param classes a character giving the classes
#' @param key_val_pairs a list of tuple of type \code{"character"}
#' @examples
#' Attr("A", c("B", "C"), list(c("D", "E")))
#' @export
## type Attr = (String, [String], [(String, String)]) 
Attr <- function(identifier="", classes=character(), key_val_pairs=list()) {
    if ( !is.character(classes) ) stop("'classes' has to be of type character!")
    structure(list(identifier, as.list(classes), key_val_pairs), class=c("Attr", "list"))
}

#  -----------------------------------------------------------
#  ListAttributes
#  ==============
#' @title ListAttributes
#' @description A constructor for pandoc list attributes.
#' @param first_number an integer giving the first number of the list
#' @param style a character string giving the style, possible values are \code{"DefaultStyle"}, 
#'              \code{"Example"}, \code{"Decimal"}, \code{"LowerRoman"}, 
#'              \code{"UpperRoman"}, \code{"LowerAlpha"} and \code{"UpperAlpha"}.
#' @param delim a character string giving the delimiter, possible values are \code{"DefaultDelim"},
#'              \code{"Period"}, \code{"OneParen"} and \code{"TwoParens"}.
#' @export
## ListAttributes = (Int, ListNumberStyle, ListNumberDelim) 
ListAttributes <- function(first_number=1L, style="DefaultStyle", delim="DefaultDelim") {
    structure(list(first_number, list(t=style, c=list()), list(t=delim, c=list())), class=c("ListAttributes", "list"))
}

#  -----------------------------------------------------------
#  TableCell
#  =========
#' @title Table Cell
#' @description Table cells is a constructor for plain table cells.
#' @param x a character string giving the content of the table cell
#' @details In general table cells are a list of block elements, the 
#'          constructor \code{TableCell} creates a plain table cell.
#' @examples
#' TableCell("Cell 1")
#' @export
#  -----------------------------------------------------------
TableCell <- function(x) structure(list(Plain(list(Str(x)))), class=c("TableCell", "list"))

## ----------------------------------------------------------------------------- 
##
##   Inline Element Constructors
##
## -----------------------------------------------------------------------------
##  1. Str
##  2. Emph
##  3. Strong
##  4. Strikeout
##  5. Superscript
##  6. Subscript
##  7. SmallCaps
##  8. Quoted
##  9. Cite
## 10. Code
## 11. Space
## 12. SoftBreak
## 13. LineBreak
## 14. Math
## 15. RawInline
## 16. Link
## 17. Image
## 18. Note
## 19. Span

#  -----------------------------------------------------------
#   1. Str
#   ======
#' @title Text (String)
#' @description A constructor of an inline object of type \code{"Str"}.
#' @param x a character string
#' @details 
#'   To minimize the amount of unnecessary typing, pandoc filters automatically 
#'   converts character strings to pandoc objects of type \code{"Str"} if needed.   
#'   Furthermore, if a single inline object is provided where a list of inline 
#'   objects is needed \pkg{pandocfilters} automatically converts this inline 
#'   object into a list of inline objects. For example, the canonical way to emphasize 
#'   the character string \code{"some text"} would be \code{Emph(list(Str("some text")))} 
#'   since single inline objects are automatically transformed to lists of inline objects, 
#'   this is equivalent to \code{Emph(Str("some text"))}. Since a character 
#'   string is automatically transformed to an inline object, this is is equivalent 
#'   to \code{Emph("some string")}. In short, whenever a list of inline objects 
#'   is needed one can also use a single inline object or a character string.
#' @examples
#' Str("SomeString")
#' @export
#  -----------------------------------------------------------
Str <- function(x) {
    structure(list(t="Str", c=x), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  2. Emph
#  =======
#' @title Emphasized Text
#' @description A constructor of an inline object of type \code{"Emph"}.
#' @param x a inline object or a list of inline objects
#' @examples
#' Emph("emphasize")
#' @export
#  -----------------------------------------------------------
Emph <- function(x) {
    structure(list(t="Emph", c=as.loio(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  3. Strong
#  =========
#' @title Strongly Emphasized Text
#' @description A constructor of an inline object of type \code{"Strong"}.
#' @param x a inline object or a list of inline objects
#' @examples
#' Strong("strong")
#' @export
#  -----------------------------------------------------------
Strong <- function(x) {
    structure(list(t="Strong", c=as.loio(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  4. Strikeout
#  ============
#' @title Strikeout Text
#' @description A constructor of an inline object of type \code{"Strikeout"}.
#' @param x a inline object or a list of inline objects
#' @examples
#' Strikeout("strikeout")
#' @export
#  -----------------------------------------------------------
Strikeout <- function(x) {
    structure(list(t="Strikeout", c=as.loio(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  5. Superscript
#  ==============
#' @title Superscripted Text
#' @description A constructor of an inline object of type \code{"Superscript"}.
#' @param x a inline object or a list of inline objects
#' @examples
#' Superscript("some text written in superscript")
#' @export
#  -----------------------------------------------------------
Superscript <- function(x) {
    structure(list(t="Superscript", c=as.loio(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  6. Subscript
#  ============
#' @title Subscripted Text
#' @description A constructor of an inline object of type \code{"Subscript"}.
#' @param x a inline object or a list of inline objects
#' @examples
#' Subscript("some text written in superscript")
#' @export
#  -----------------------------------------------------------
Subscript <- function(x) {
    structure(list(t="Subscript", c=as.loio(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  7. SmallCaps
#  ============
#' @title Small Caps Text
#' @description A constructor of an inline object of type \code{"SmallCaps"}.
#' @param x a inline object or a list of inline objects
#' @examples 
#' SmallCaps("The latex command for 'small caps' is 'textsc'!")
#' @export
#  -----------------------------------------------------------
SmallCaps <- function(x) {
    structure(list(t="SmallCaps", c=as.loio(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  8. Quoted
#  =========
#' @title Quoted Text
#' @description A constructor of an inline object of type \code{"Quoted"}.
#' @param x a inline object or a list of inline objects
#' @param quote_type a character giving the quote type,
#'                   valid types are \code{"SingleQuote"} and \code{"DoubleQuote"}
#' @examples
#' Quoted("some text", quote_type="SingleQuote")
#' Quoted("some text", quote_type="DoubleQuote")
#' @export
#  -----------------------------------------------------------
## Quoted QuoteType [Inline]
## Quoted text (list of inlines)
Quoted <- function(x, quote_type="DoubleQuote") {
    structure(list(t="Quoted", c=list(list(t=quote_type, c=list()), as.loio(x))), 
        class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  9. Cite
#  =======
#' @title Citation
#' @description A constructor of an inline object of type \code{"Cite"}.
#' @param citation an object of type \code{"Citation"}
#' @param x a inline object or a list of inline objects
#' @examples
#' ci <- Citation(suffix=list(Str("Suffix_1")),
#'                id="Citation_ID_1", prefix=list(Str("Prefix_1")))
#' Cite(ci, Str("some text"))
#' @export
## Cite [Citation] [Inline] 
## Citation (list of inlines)
#  -----------------------------------------------------------
Cite <- function(citation, x) {
    if ( is.citation(citation) ) citation <- list(citation)
    structure(list(t="Cite", c=list(citation, as.loio(x))), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  10. Code
#  ========
#' @title Inline Code 
#' @description A constructor of an inline object of type \code{"Code"}.
#' @param code a character string giving the inline code
#' @param name an optional character string giving the name of the inline code chunk
#' @param language an optional character string giving the programming language
#' @param line_numbers a logical which controls if line numbers should be used
#' @param start_from an integer giving the first line number
#' @examples
#' Code("lm(hello ~ world)", "my_r_inline_code", "R", TRUE, 0)
#' Code("lm(hello ~ world)")
#' @export
## Additional material (from the pandoc homepage)
## ===================
## type Attr = (String, [String], [(String, String)])
## Attributes: identifier, classes, key-value pairs
## one example on the pandoc homepage shows the following for
## inline code: `<$>`{.haskell}
## block code:
##     ~~~~ {#mycode .haskell .numberLines startFrom="100"}
##     qsort []     = []
##     qsort (x:xs) = qsort (filter (< x) xs) ++ [x] ++
##                    qsort (filter (>= x) xs)
##     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     Here mycode is an identifier, haskell and numberLines are classes,
##     and startFrom is an attribute with value 100.
## or in html:
##    <pre id="mycode" class="haskell numberLines" startFrom="100"><code>... </code></pre>
#  -----------------------------------------------------------
Code <- function(code, name="", language=NULL, line_numbers=FALSE, start_from=1) {
    if (line_numbers) {
        linum <- list(list("startFrom", sprintf("%i", start_from)))
    } else {
        linum <- list()
    }
    if ( !is.null(language) ) {
        lang <- list(language)
    } else {
        lang <- list()
    }
    meta <- list(name, lang, linum)
    x <- list(meta, code)
    structure(list(t="Code", c=x), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  11. Space
#  =========
#' @title Inter-word space
#' @description A constructor of an inline object of type \code{"Space"}.
#' @examples
#' Space()
#' @export
#  -----------------------------------------------------------
Space <- function() structure(Type("Space"), class=c("inline", "list"))

#  -----------------------------------------------------------
#  12. SoftBreak
#  =============
#
#  NOTE: SoftBreak, the created data structure should be correct but
#        I couldn't figure out what it actually does.
#
#' @title Soft Line Break
#' @description A constructor of an inline object of type \code{"SoftBreak"}.
#' @examples
#' SoftBreak()
#' @export
#  -----------------------------------------------------------
SoftBreak <- function() structure(Type("SoftBreak"), class=c("inline", "list"))

#  -----------------------------------------------------------
#  13. LineBreak
#  =============
#' @title Hard Line Break
#' @description A constructor of an inline object of type \code{"LineBreak"}.
#' @examples
#' LineBreak()
#' @export
#  -----------------------------------------------------------
LineBreak <- function() structure(Type("LineBreak"), class=c("inline", "list"))

#  -----------------------------------------------------------
#  14. Math
#  ========
#' @title TeX Math
#' @description A constructor of an inline object of type \code{"Math"}.
#' @param x a character string 
#' @examples
#' Math("3*x^2")
#' @export
#  -----------------------------------------------------------
Math <- function(x) structure(list(t="Math", c=list(Type("InlineMath"), x)), class=c("inline", "list"))

#  -----------------------------------------------------------
#  15. RawInline
#  =============
#' @title Raw Inline
#' @description A constructor of an inline object of type \code{"RawInline"}.
#' @param format a character string giving the format (e.g. \code{"latex"}, \code{"html"})
#' @param x a character string giving the inline
#' @examples
#' RawInline("latex", "some RawInline")
#' @export
#  -----------------------------------------------------------
RawInline <- function(format, x) {
    structure(list(t="RawInline", c=list(format, x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  16. Link
#  ========
#' @title Hyperlink
#' @description A constructor of an inline object of type \code{"Link"}.
#' @param target a character string giving the target (hyper reference)
#' @param text a inline object or a list of inline objects giving the visible part
#' @param title an optional character string giving the title
#' @param attr an optional object of type \code{"Attr"}
#' @details Further Usage examples can be found in the README.
#' @examples
#' Link("https://cran.r-project.org/", "Text_Shown", "some title")
#' @export
## <A HREF="url.html" TITLE="some_alterinative_text">the_text_shown</A>
## Link Attr [Inline] Target | Hyperlink: "alt text" (list of inlines), target
#  -----------------------------------------------------------
Link <- function(target, text, title="", attr=Attr()) {
    if ( get_pandoc_version() < 1.16 ) {
        return( structure(list(t="Link", c=list(as.loio(text), 
                                                list(target, title))), 
                          class=c("inline", "list")) )
    }
    structure(list(t="Link", c=list(attr,
                                    as.loio(text), 
                                    list(target, title))), 
              class=c("inline", "list"))
}
## Link <- function(target, text, alt_text="") {
##     structure(list(t="Link", c=list(as.loio(text), list(target, alt_text))), 
##         class=c("inline", "list"))
## }

#  -----------------------------------------------------------
#  17. Image
#  =========
#' @title Image
#' @description A constructor of an inline object of type \code{"Image"}.
#' @param target a character string giving the target (hyper reference)
#' @param text a inline object or a list of inline objects giving the visible part
#' @param caption a character string describing the picture
#' @param attr an optional object of type \code{"Attr"}
#' @details Further Usage examples can be found in the README.
#' @examples
#' Image("https:://Rlogo.jpg", "some_text", "fig:some_caption")
#' @export
#  -----------------------------------------------------------
## Image Attr [Inline] Target | Image: alt text (list of inlines), target
Image <- function(target, text, caption="", attr=Attr()) {
    if ( get_pandoc_version() < 1.16 ) {
        return( structure(list(t="Image", c=list(as.loio(text), 
                                                 list(target, caption))), 
                          class=c("inline", "list")) )
    }
    structure(list(t="Image", c=list(attr,
                                     as.loio(text), 
                                     list(target, caption))), 
              class=c("inline", "list"))
}
## Image <- function(target, text, caption="") {
##     structure(list(t="Image", c=list(as.loio(text), list(target, caption))), 
##         class=c("inline", "list"))
## }


#  -----------------------------------------------------------
#  18. Note
#  ========
#' @title Note
#' @description A constructor of an inline object of type \code{"Note"}.
#' @param x a pandoc block object or a list of pandoc block objects
#' @examples
#' block <- Plain("x")
#' Note(block)
#' @export
## Note [Block]
## note <- Note(block)
## pandocfilters:::test(list(Plain(note)))
Note <- function(x) {
    structure(list(t="Note", c=as.lobo(x)), class=c("inline", "list"))
}

#  -----------------------------------------------------------
#  19. Span
#  ========
#' @title Generic Inline Container with Attributes
#' @description A constructor of an inline object of type \code{"Span"}.
#' @param attr an object of type \code{"Attr"}
#' @param inline a inline object or a list of inline objects which will be shown
#' @examples
#' attr <- Attr("A", "B", list(c("C", "D")))
#' Span(attr, "some inline string")
#' @export
#  -----------------------------------------------------------
Span <- function(attr, inline) {
    structure(list(t="Span", c=list(attr, as.loio(inline))), class=c("inline", "list"))
}

## ----------------------------------------------------------------------------- 
##
##   Block Element Constructors
##
## -----------------------------------------------------------------------------

##  1. Plain
##  2. Para
##  3. CodeBlock
##  4. RawBlock
##  5. BlockQuote
##  6. OrderedList
##  7. BulletList
##  8. DefinitionList
##  9. Header
## 10. HorizontalRule
## 11. Table
## 12. Div
## 13. Null

## NOTE:


#  -----------------------------------------------------------
#  1. Plain
#  ========
#' @title Plain Text
#' @description A constructor of a block object of type \code{"Plain"}.
#' @param x a inline object or list of inline objects
#' @examples
#' Plain("x")
#' @export
#  -----------------------------------------------------------
Plain <- function(x) {
    structure(list(t="Plain", c=as.loio(x)), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  2. Para
#  =======
#' @title Paragraph
#' @description A constructor of a block object of type \code{"Para"}.
#' @param x a inline object or list of inline objects
#' @examples
#' Para("x")
#' @export
#  -----------------------------------------------------------
Para <- function(x) {
    structure(list(t="Para", c=as.loio(x)), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  3. CodeBlock
#  ============
#' @title Code Block
#' @description A constructor of a block object of type \code{"CodeBlock"}.
#' @param attr an object of type \code{"Attr"}
#' @param code a character string containing the source code.
#' @examples
#' attr <- Attr("id", "Programming Language", list(c("key", "value")))
#' code <- "x <- 3\nprint('Hello R!')"
#' CodeBlock(attr, code)
#' @export
## Attr String
#  -----------------------------------------------------------
CodeBlock <- function(attr, code) {
    structure(list(t="CodeBlock", c=list(attr, code)), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  4. RawBlock
#  ===========
# #TODO: RawBlock <- elt('RawBlock', 2)
# #NOTE: Currently not implemented since

#  -----------------------------------------------------------
#  5. BlockQuote
#  =============
#' @title Block Quote
#' @description A constructor of a block object of type \code{"BlockQuote"}.
#' @param blocks a block object or list of block objects
#' @examples
#' BlockQuote(Plain("Hello R!"))
#' @export
## Attr String
#  -----------------------------------------------------------
BlockQuote <- function(blocks) {
    structure(list(t="BlockQuote", c=as.lobo(blocks)), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  6. OrderedList
#  ==============
#' @title Ordered List
#' @description A constructor of a block object of type \code{"OrderedList"}.
#' @param lattr a list of attributes
#' @param llblocks a list of lists of blocks
#' @examples
#' ordered_1 <- Plain("A")
#' ordered_2 <- list(Plain(Str("B")))
#' ordered_3 <- list(Plain(list(Str("C"))))
#' OrderedList(ListAttributes(), ordered_1)
#' OrderedList(ListAttributes(), list(ordered_1, ordered_2, ordered_3))
#' @export
## Attr String
#  -----------------------------------------------------------
OrderedList <- function(lattr, llblocks) {
    structure(list(t="OrderedList", c=list(lattr, as.lolobo(llblocks))), class=c("block", "list"))
}

## attach(getNamespace("pandocfilters"))

#  -----------------------------------------------------------
#  7. BulletList
#  =============
#' @title Bullet List
#' @description A constructor of a block object of type \code{"BulletList"}.
#' @param llblocks a list of lists of blocks
#' @examples
#' bullet_1 <- Plain("A")
#' bullet_2 <- Plain(Str("B"))
#' bullet_3 <- list(Plain(list(Str("C"))))
#' BulletList(list(bullet_1, bullet_2, bullet_3))
#' @export
## Attr String
#  -----------------------------------------------------------
BulletList <- function(llblocks) {
    structure(list(t="BulletList", c=as.lolobo(llblocks)), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  8.0 Definition
#  ==============
#' @title Definition
#' @description A constructor of a \code{Definition} which can be used as an element of a \code{"DefinitionList"}.
#' @param key a inline object or list of inline objects 
#' @param value a block object or list of block objects
#' @examples
#' Definition("some key", Plain("some value"))
#' @export
#  -----------------------------------------------------------
Definition <- function(key, value) {
    list(as.loio(key), as.lolobo(value))
}

#  -----------------------------------------------------------
#  8. DefinitionList
#  =================
#' @title Definition List
#' @description A constructor of a block object of type \code{"DefinitionList"}.
#' @param x a list of key value pairs, the key is a list of \code{"inline"} objects and
#'          the values are a list of lists of objects of type \code{"block"}.
#' @details In the pandoc API \url{http://johnmacfarlane.net/BayHac2014/doc/pandoc-types/Text-Pandoc-Definition.html} 
#'          the \code{DefinitionList} is described as follows, each list item is a pair consisting of a term 
#'          (a list of \code{"inline"} objects) and one or more definitions (each a list of blocks).
#' @examples
#' key <- list(Str("key"))
#' value <- list(list(Plain(list(Str("value")))))
#' DefinitionList(list(list(key, value), Definition("some key", Plain("some value"))))
#' @export
## Attr String
#  -----------------------------------------------------------
DefinitionList <- function(x) {
    structure(list(t="DefinitionList", c=x), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  9. Header
#  =========
#' @title Header
#' @description A constructor of a block object of type \code{"Header"}.
#' @param x a inline object or a list of inline objects
#' @param level an integer giving the level
#' @param attr an object of type \code{"Attr"}
#' @examples
#' Header("My Header")
#' @export
#  -----------------------------------------------------------
Header <- function(x, level=1L, attr=Attr()) {
    structure(list(t="Header", c=list(level, attr, as.loio(x))), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  10. HorizontalRule
#  ==================
#' @title Horizontal Rule
#' @description A constructor of a block object of type \code{"HorizontalRule"}.
#' @examples
#' HorizontalRule()
#' @export
## Attr String
#  -----------------------------------------------------------
HorizontalRule <- function() {
    structure(list(t="HorizontalRule", c=list()), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  11. Table
#  =========
#' @title Table
#' @description A constructor of a block object of type \code{"Table"}.
#' @param rows an object of class \code{"matrix"}, \code{"data.frame"}, \code{"table"} 
#'   or a list of lists of pandoc objects of type \code{"TableCell"}
#' @param col_names a list of objects of type \code{"TableCell"}
#' @param aligns a character vector of alignments, possible values are \dQuote{l} for left,
#'               \dQuote{r} for right, \dQuote{c} for center and \dQuote{d} for default.
#' @param col_width a numeric vector
#' @param caption a inline object or a list of inline objects giving the caption
#' @details Table, with caption, column alignments (required), relative column widths 
#'          (0 = default), column headers (each a list of blocks), 
#'          and rows (each a list of lists of blocks)
#' @export
## Table [Inline] [Alignment] [Double] [TableCell] [[TableCell]]
#  -----------------------------------------------------------
Table <- function(rows, col_names=NULL, aligns=NULL, col_width=NULL, caption=list() ) {

    if ( is.null(col_names) & (! is.null(colnames(rows))) ) {
        col_names <- colnames(rows)
    }
    
    if ( is.matrix(rows) | is.data.frame(rows) | is.table(rows) ) {

        if ( is.table(rows) ) {
            rows <- as.matrix(rows)
            if ( min(dim(rows)) == 1 ) {
                rows <- t(rows)
            }
        }

        col_fun <- function(m, n) TableCell(as.character(rows[[m, n]]))
        row_fun <- function(m) lapply(seq_len(ncol(rows)), function(n) col_fun(m, n))
        rows <- lapply(seq_len(nrow(rows)), row_fun)
    }

    number_of_columns <- length( rows[[1]] )
    if ( is.null(col_names) ) col_names <- rep("", number_of_columns)

    if ( length(col_names) ==  number_of_columns ) {
        col_names <- lapply(col_names, function(x) TableCell(paste(x)))
    } else {
        msg <- sprintf("argument 'col_names' has length %i but the Table has %i columns.", 
                       length(col_names), number_of_columns)
        stop(msg, "The number of columns have to match the number of 'col_names'.")
    }
    
    if ( is.null(aligns) ) {
        aligns <- rep("d", number_of_columns)
    } else {
        if ( length(aligns) != number_of_columns ) {
            msg <- sprintf("argument 'aligns' has length %i but the Table has %i columns.", 
                           length(aligns), number_of_columns)
            stop(msg, "The number of columns have to match the number of 'aligns'.")
        }
    }
    if ( is.null(col_width) ) {
        col_width <- integer(number_of_columns)
    } else {
        if ( length(col_width) != number_of_columns ) {
            msg <- sprintf("argument 'col_width' has length %s but the Table has %i columns.", 
                           length(col_width), number_of_columns)
            stop(msg, "The number of columns have to match the number of 'col_width'.")
        }
    }

    alignments <- setNames(c("AlignLeft", "AlignRight", "AlignCenter", "AlignDefault"), 
                           c("l", "r", "c", "d") )
    if ( !all(aligns %in% names(alignments)) ) {
        stop("wrong alignment, possible values are 'l', 'r', 'c' or 'd'")
    }
    aligns <- unname(lapply(alignments[aligns], FUN=function(x) list(t=unname(x), c=list())))
    if ( is.character(caption) ) {
        caption <- Str(caption)
    }
    structure(list(t="Table", c=list(as.loio(caption), aligns, 
                   as.list(col_width), col_names, rows)), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  12. Div
#  =======
#' @title Generic Block Container with Attributes
#' @description A constructor of a block object of type \code{"Div"}.
#' @param blocks a block object or list of block objects
#' @param attr an object of type \code{"Attr"}
#' @examples
#' blocks <- Plain("Hello R!")
#' Div(blocks)
#' @export
## Attr String
#  -----------------------------------------------------------
Div <- function(blocks, attr=Attr()) {
    structure(list(t="Div", c=list(attr, as.lobo(blocks))), class=c("block", "list"))
}

#  -----------------------------------------------------------
#  13. Null
#  ========
#' @title Nothing
#' @description A constructor of a block object of type \code{"Null"}.
#' @examples
#' Null()
#' @export
#  -----------------------------------------------------------
Null <- function() structure(list(t="Null", c=list()), class=c("block", "list"))

##' @noRd
##' @export
print.inline <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.loio <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.block <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.TableCell <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.ListAttributes <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.Attr <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.Type <- function(x, ...) print(unclass(x))

##' @noRd
##' @export
print.Citation <- function(x, ...) print(unclass(x))
