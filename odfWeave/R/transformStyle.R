# This function makes a new style from a specified style,
# with the specified name.  The "type" argument specifies
# whether the speccified style is common or automatic.
# The resulting style will cause a page break to occur
# before it.
makePageBreakStyle <- function(name, family=c('paragraph', 'table'),
                               type=c('common', 'automatic'), prevstyle)
{
   family <- match.arg(family)
   type <- match.arg(type)
   propertiesName <- sprintf('style:%s-properties', family)

   if (type == 'common')
   {
      if (is.null(prevstyle) || ! is.character(prevstyle))
         stop('prevstyle must be a style name if type is common')

      # Create the new style
      xmlNode('style:style',
              xmlNode(propertiesName,
                      attrs=c('fo:break-before'='page')),
              attrs=c('style:name'=name,
                      'style:family'=family,
                      'style:parent-style-name'=prevstyle))
   } else if (type == 'automatic') {
      if (is.null(prevstyle) || ! inherits(prevstyle, 'XMLNode'))
         stop('prevstyle must be a style node if type is automatic')

      # Extract the attributes and children from the previous style object
      attrs <- xmlAttrs(prevstyle)
      children <- xmlChildren(prevstyle)

      # Add/update the style name
      attrs['style:name'] <- name

      # Create a new child list from the old
      childrenNames <- unlist(lapply(children, function(c) xmlName(c, full=TRUE)))
      iprops <- which(childrenNames == propertiesName)
      newchildren <- if (length(iprops) > 0)
      {
         if (length(iprops) > 1)
            warning(sprintf('found style with multiple %s children', propertiesName))
         xnode <- children[[iprops[1]]]
         pattrs <- xmlAttrs(xnode)
         pattrs['fo:break-before'] <- 'page'
         xnode$attributes <- pattrs
         c(children[-iprops], list(xnode))
      } else {
         c(children, list(xmlNode(propertiesName, attrs=c('fo:break-before'='page'))))
      }

      # Set the updated attributes and children, and return the modified object
      prevstyle$attributes <- attrs

      # XXX work-around bug in XML 3.4?
      # xmlChildren(prevstyle) <- newchildren
      prevstyle <- makeNode(prevstyle, newchildren)

      # Return the modified node
      prevstyle
   } else {
      stop('illegal style type specified: ', type)
   }
}

# This function makes a new style from a specified style,
# with the specified name.  The "type" argument specifies
# whether the speccified style is common or automatic.
# The resulting style will use the specified page style.
makeSetPageStyle <- function(name, pagestyle, family=c('paragraph', 'table'),
                             type=c('common', 'automatic'), prevstyle)
{
   family <- match.arg(family)
   type <- match.arg(type)

   if (type == 'common')
   {
      if (is.null(prevstyle) || ! is.character(prevstyle))
         stop('prevstyle must be a style name if type is common')

      # Create the new style
      xmlNode('style:style',
              attrs=c('style:name'=name,
                      'style:family'=family,
                      'style:parent-style-name'=prevstyle,
                      'style:master-page-name'=pagestyle))
   } else if (type == 'automatic') {
      if (is.null(prevstyle) || ! inherits(prevstyle, 'XMLNode'))
         stop('prevstyle must be a style node if type is automatic')

      if (family != xmlGetAttr(prevstyle, 'style:family'))
         stop('family argument must match family of prevstyle')

      # Extract the attributes and children from the previous style object
      attrs <- xmlAttrs(prevstyle)

      # Add/update the style name
      attrs['style:name'] <- name
      attrs['style:master-page-name'] <- pagestyle

      # Set the updated attributes
      prevstyle$attributes <- attrs

      # Return the modified node
      prevstyle
   } else {
      stop('illegal style type specified: ', type)
   }
}
