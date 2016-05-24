
trim =
function(x)
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)


textNodesOnly =
  # Only process the top-level text nodes, not recursively.
  # Could be done as simply as
  #     xmlValue(x, recursive = FALSE)
function(x)
    paste(xmlSApply(x, function(n) if(is(n, "XMLInternalTextNode")) xmlValue(n) else ""), collapse = "")


toNumber =
function(x)
{
  as.numeric(gsub("[%,]", "", x))
}


if(FALSE) {
  doc = htmlParse("http://elections.nytimes.com/2008/results/states/president/california.html")
  tbls = getNodeSet(doc, "//table[not(./tbody)]|//table/tbody")
  o = readHTMLTable(tbls[[1]], skip.rows = c(1, Inf), header = FALSE, colClasses = c("character", replicate(5, toNumber)), elFun = textOnly)


  o = readHTMLTable("http://elections.nytimes.com/2008/results/states/president/california.html")  

  x = readHTMLTable("http://www.usatoday.com/news/politicselections/vote2004/CA.htm", as.data.frame = FALSE)  
}

setGeneric("readHTMLTable",
          function(doc, header = NA,
                    colClasses = NULL, skip.rows = integer(), trim = TRUE, elFun = xmlValue,
                     as.data.frame = TRUE, which = integer(), ...)           
             standardGeneric("readHTMLTable"))

setMethod("readHTMLTable", "character",
          function(doc, header = NA,
                    colClasses = NULL, skip.rows = integer(), trim = TRUE, elFun = xmlValue,
                     as.data.frame = TRUE, which = integer(), encoding = character(), ...) {
              pdoc = htmlParse(doc, encoding = encoding)
              readHTMLTable(pdoc, header, colClasses, skip.rows, trim, elFun, as.data.frame, which, ...)
          })


 # XXX Should vectorize in header, colClasses, i.e. allow different values for different tables.
setMethod("readHTMLTable", "HTMLInternalDocument",
          function(doc, header = NA,
                    colClasses = NULL, skip.rows = integer(), trim = TRUE, elFun = xmlValue,
                     as.data.frame = TRUE, which = integer(), ...)           
{
     #  tbls = getNodeSet(doc, "//table[not(./tbody)]|//table/tbody")
   tbls = getNodeSet(doc, "//table")  # XXX probably want something related to nested tables
                                      # "//table[not(ancestor::table)]" -> outer ones
      # if header is missing, compute it each time.
   if(length(which))
       tbls = tbls[which]

#   ans = lapply(tbls, readHTMLTable,  header, colClasses, skip.rows, trim, elFun, as.data.frame, ...)
   header = rep(header, length = length(tbls))
   ans = mapply(readHTMLTable,
                 tbls, header,
                 MoreArgs = list(colClasses = colClasses, skip.rows = skip.rows, trim = trim, elFun = elFun, as.data.frame = as.data.frame, ...),
                 SIMPLIFY = FALSE)
   names(ans) = sapply(tbls, getHTMLTableName)

   if(length(which) && length(tbls) == 1)
      ans[[1]]
   else
      ans
})

getHTMLTableName =
function(node)
{
  id = xmlGetAttr(node, "id")
  if(!is.null(id))
    return(id)

  cap = getNodeSet(node, "./caption")
  if(length(cap))
    return(xmlValue(cap[[1]]))
}

setClass("FormattedNumber", contains = "numeric")
setClass("FormattedInteger", contains = "integer")

setAs('character', 'FormattedNumber', function(from) as.numeric(gsub(",", "", from)))
setAs('character', 'FormattedInteger', function(from) as.integer(gsub(",", "", from)))

setClass("Currency", contains = "numeric")
setAs("character", "Currency",
       function(from)
          as.numeric(gsub("[$,]", "", from)))

setClass("Percent", contains = "numeric")
setAs('character', 'Percent', function(from) as.numeric(gsub("%", "", from)))

setMethod("readHTMLTable", "XMLInternalElementNode",
#readHTMLTable.XMLInternalElementNode  =
  #
  #
  # header is computed based on whether we have a table node and it has a thead.
  #  (We don't currently bother with the col spans.)
  #
  #  colClasses can be a character vector giving the name of the type for a column,
  #  an NULL to drop the corresponding column, or a function in which case it will
  #  be passed the contents of the column and can transform it as it wants.
  #  This allows us to clean text before converting it.
  #
  # skip.rows - an integer vector indicating which rows to ignore.
  #
  #  trim - a logical indicating whether to trim white space from the start and end of text.
  #
  #  elFun - a function which is called to process each th or td node to extract the content.
  #      This is typically xmlValue, but one can supply others (e.g. textNodesOnly)

  #  as.data.frame
  #

function(doc, header = NA ,
          colClasses = NULL, skip.rows = integer(), trim = TRUE, elFun = xmlValue,
            as.data.frame = TRUE, encoding = 0L, ...)
{

  node = doc
  headerFromTable = FALSE
  dropFirstRow = FALSE

  
     # check if we have a header
  if(length(header) == 1 && is.na(header))                                     # this node was doc
      header = (xmlName(doc) %in% c("table", "tbody") &&
                      ("thead" %in% names(doc) || length(getNodeSet(node, "./tr[1]/th | ./tr[1]/td")) > 0))

  if(is.logical(header) && (is.na(header) || header) &&  xmlName(node) == "table") {
    if("thead" %in% names(node))
       header = node[["thead"]]
    else {
       if("tr" %in% names(node))
          tmp = node[["tr"]]
       else
          tmp = node[["tbody"]][["tr"]]

       if(!is.null(tmp) && all(names(tmp) %in% c('text', 'th'))) {
          header = xpathSApply(tmp, "./th | ./td", xmlValue, encoding = encoding)
          dropFirstRow = TRUE
        }
    }
  }

     # Moved this from before the check for header as we set node here and that seems
     # premature. Checked on
     #     readHTMLTable("http://www.google.com/finance?q=NASDAQ:MSFT&fstype=ii", header = TRUE, which = 1)
  tbody = getNodeSet(node, "./tbody")
  if(length(tbody))
     node = tbody[[1]]  

  if(is(header, "XMLInternalElementNode"))   {
      # get the last tr in the thead
     if(xmlName(header) == "thead") {
        i = which(names(header) == "tr")
        header = header[[ i [ length(i) ] ]]
        xpath = "./th | ./td"
     } else
        xpath = "./*/th | ./*/td"

     header = as.character(xpathSApply(header, xpath, elFun, encoding = encoding))
     headerFromTable = TRUE

     if(xmlName(node) == "table" && "tbody" %in% names(node))
        node = node[["tbody"]]
   }
  
     # Process each row, by getting the content of each "cell" (th/td)
  rows = getNodeSet(node, ".//tr")
  if(dropFirstRow)
     rows = rows[-1]
  els =  lapply(rows, function(row) {
                           tmp = xpathSApply(row, "./th|./td", elFun)
                           if(trim)
                              trim(tmp)
                           else
                              tmp
                        })


#  spans = getNodeSet(node, ".//td[@rowspan] | .//th[@rowspan]")

  
  if(length(skip.rows)) {
    infs = (skip.rows == Inf)
    if(any(infs))
          # want Inf - 2, Inf - 1, Inf,  to indicate drop last 3, but that won't work
          # take sequence of Inf to identify Inf - 2, Inf - 1, Inf
       skip.rows[skip.rows == Inf] = length(els)  - seq(0, length = sum(infs))  
    els = els[ - skip.rows ]
  }

  if(length(els) == 0)
    return(NULL)
  
   numEls = sapply(els, length)
                                                           # els[[1]] should be a scalar
   if(is.logical(header) && !is.na(header) && header && any(nchar(els[[1]]) < 999)) {
     header = els[[1]]
     els = els[-1]
     numEls = numEls[ - 1]
   }

  if(length(els) == 0)
    return(NULL)  #XXX we should have a header here so return a data frame with 0 rows.

   ans = lapply(seq(length = max(numEls)),
                  function(col) {
                    sapply(els, `[`, col)
                  })

   if(is.character(header) && length(header) == length(ans))
      names(ans) = header

   if(length(colClasses)) {

      colClasses = rep(colClasses, length = length(ans))
     
      n = sapply(colClasses, is.null)
      if(any(n)) {
         ans = ans[ ! n ]
         colClasses = colClasses[ ! n ]
      }

      ans = lapply(seq(along = ans) ,
                      function(i) 
                         if(is.function(colClasses[[i]]))
                            colClasses[[i]](ans[[i]])
                         else if(colClasses[[i]] == "factor")
                             factor(ans[[i]])
                         else if(colClasses[[i]] == "ordered")
                             ordered(ans[[i]])        
                         else 
                            as(ans[[i]], colClasses[[i]])
                   )

   }

   if(as.data.frame)  {
     ans = as.data.frame(ans, ...)
     if(is.character(header) && length(header) == length(ans))
        names(ans) = header
     else if(nrow(ans) > 0)
       names(ans) = paste("V", seq(along = ans), sep = "")
   }
    
  ans
})



getTableWithRowSpan =
function(node, r = xmlSize(node),
          c = max(xmlSApply(node, function(x) length(getNodeSet(x, "./td | ./th")))),
           encoding = 0L)
{
  ans = matrix(NA_character_, r, c)
  for(i in seq(length = r)) {
     col = 1
     kids = getNodeSet(node[[i]], "./th | ./td")
     for(k in seq(along = kids)) {
       sp = xmlGetAttr(k, "rowspan", 1)
       ans[seq(i, length = sp)] = xmlValue(k, encoding = encoding)
     }
  }
}
