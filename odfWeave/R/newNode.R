makeNode <- function(node, newChildren) {
   if (inherits(node, 'XMLTextNode')) {
      node
   } else {
      name <- xmlName(node)
      attrs <- xmlAttrs(node)
      namespace <- xmlNamespace(node)
      namespaceDefinitions <- xmlNamespaceDefinitions(node)
      xmlNode(name,
              attrs=attrs,
              namespace=namespace,
              namespaceDefinitions=namespaceDefinitions,
              .children=newChildren)
   }
}
