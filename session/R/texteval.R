# $Id: texteval.R 198 2003-04-01 22:36:07Z warnes $
#
# $Log$
# Revision 1.5  2003/04/01 22:36:07  warnes
# - Added 'capture' function that grabs the result of evaluating an expression.
# - Modified texteval to use capture.
#
# Revision 1.4  2003/04/01 19:54:09  warnes
#
# - Add 'on.exit' so that the output sink will be terminated if there is
#   an error during the evaluation of the sourceText code.
#
# Revision 1.3  2002/06/20 17:01:57  warnesgr
# Fixed typo in printed()
#
# Revision 1.2  2002/06/12 15:57:28  warnesgr
#
# - Belated checkin of changes to texteval.R to add the function
#   'printed' which does not include the evaluated commands in the output.
# - Updated the texteval.Rd help file.
#
# Revision 1.1  2002/04/30 14:20:01  warneg
#
# - Added texteval function & documentation
# - Added keywords to session.Rd
#

capture <- function( expression, collapse=NULL)
  {
    resultConn <- textConnection("resultText", open="w")
    sink(resultConn)
    on.exit( function() { sink(); close(resultConn) } )

    expression

    on.exit( NULL )
    sink()
    close(resultConn)

    return( paste( c(resultText, ""), collapse=collapse, sep="" ) )
    # the reason for c(result, "") is so that we get the line
    # terminator on the last line of output.  Otherwise, it just shows
    # up between the lines.
  }

texteval <- function( sourceText, collapse=NULL, echo=TRUE )
  {
    sourceConn <- textConnection(sourceText, open="r")
    on.exit( close(sourceConn) )
    
    result <- capture(
                      source(file=sourceConn, local=FALSE,
                             echo=echo, print.eval=TRUE),
                      collapse=collapse)

    on.exit( NULL )
    close(sourceConn)
    return( result )
  }


printed <- function( sourceText, collapse=NULL )
  {
    return( texteval(sourceText, collapse, FALSE) )
  }

