#' setClass for ecdb class
#' 
#' setClass for ecdb class
#'
#' @slot call the match.call slot
#' @slot file character, the full path to an elliptic database.
#' @slot conn an object of SQLiteConnection class.
#' @slot is.internal logical, whether the connected db is internal.
#' @slot conf list of configuration for data generation assigned by the constructor.
#'            Typical user should not have to modify this list 
#'            unless you need to generate more data for advanced research.            
#'
#' @keywords ecdb
#'
#' @author Stephen H-T. Lihn
#'
#' @exportClass ecdb
setClass("ecdb",
         representation(call = "call",
                        file = "character",
                        conn = "SQLiteConnection",
                        is.internal = "logical",
                        conf = "list"),
          prototype(call = call("ecdb"),
                    file = "elliptic.db",
                    conn = NULL,
                    is.internal = NA,
                    conf = list())
)
### <---------------------------------------------------------------------->
