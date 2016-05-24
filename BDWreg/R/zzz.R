.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0(
        '\n ========================================================================',
        '\n If you have any question about this package and corresponding paper use ',
		    '\n      hamedhaseli@gmail.com  or visit www.hamedhaseli.webs.com           ',
		    '\n ========================================================================'
         ),
      domain = NULL,  appendLF = TRUE )
}
