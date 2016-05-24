#' Initialize operators 
#' @rdname initOps
.initOps <- function() {

  setOperator( '::',  type = 'namespace' )
  setOperator( ':::', type = 'namespace' )
 
  setOperator( '@',   type = 'component' )
  setOperator( '$',   type = 'component' )
 
  setOperator( '[',   type = 'indexing' )
  setOperator( '[[',  type = 'indexing' )
 
  setOperator( ':',  type = 'sequence' )
 
  setOperator( '+',  type = 'arithmetic', inverse = '-' )
  setOperator( '-',  type = 'arithmetic', inverse = '+' )
  setOperator( '*',  type = 'arithmetic', inverse = '/' )
  setOperator( '/',  type = 'arithmetic', inverse = '*' )
  setOperator( '^',  type = 'arithmetic' )  # inverse = as.name('log') 
  setOperator( '%%', type = 'arithmetic' )
  setOperator( '%/%', type = 'arithmetic' )


  setOperator( '<',    type = 'relational', inverse = '>=', rel.type = 'lt' )
  setOperator( '<=',   type = 'relational', inverse = '>',  rel.type = 'lt' )
  setOperator( '>',    type = 'relational', inverse = '<=', rel.type = 'gt' )
  setOperator( '>=',   type = 'relational', inverse = '<',  rel.type = 'gt' )
  setOperator( '==',   type = 'relational', inverse = '!=', rel.type = 'eq' )
  setOperator( '!=',   type = 'relational', inverse = '==', rel.type = 'ne' )
  setOperator( '%in%', type = 'relational', rel.type = 'eq' )  # %!in% 
  setOperator( '%in%' , type = 'relational', inverse='%!in%', rel.type = 'eq' )
  setOperator( '%!in%', type = 'relational', inverse='%in%' , rel.type = 'ne' )


  setOperator( '!',  type = 'logical' ) # inverse = identity, !?    
  setOperator( '&',  type = 'logical' )
  setOperator( '&&', type = 'logical' )
  setOperator( '|',  type = 'logical' )
  setOperator( '||', type = 'logical' )

  setOperator( '~', type = 'tilde' )

  # Cannot work on rightward assingment "syntactic equivalents" ->, ->>
  setOperator( '<-',  type = 'assignment' )  # inverse rm
  setOperator( '<<-', type = 'assignment' )  # 
  setOperator( '=',   type = 'assignment' )

  setOperator( '?', type = 'help' )

  setOperator( '%*%', type = 'matrix' )
  setOperator( '%x%', type = 'matrix' )
  setOperator( '%o%', type = 'matrix' )

}
