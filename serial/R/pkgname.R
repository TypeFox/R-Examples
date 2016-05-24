#' A serial communication interface for R.
#' 
#' @description
#' This R package provides the functionality to use the serial communication ports
#' "COM" or "tty" to use the RS232/RS422/RS485 functionality of the corresponding 
#' hardware. Also virtual COM-ports via USB do work, as long as they are mapped
#' to COM[n] (win) or tty[n] (Linux) in the operating system.
#' 
#' \describe{
#'  \item{\code{open(con)}}{opens a serial connection}
#'  \item{\code{close(con)}}{closes the serial connection}
#'  \item{\code{read.serialConnection(con)}}{byte wise read from the interface as long as the buffer is empty}
#'  \item{\code{write.serialConnection(con,dat)}}{writes a string to the serial interface}
#'  \item{\code{isOpen(con)}}{test a connection, whether it is open or not}
#'  \item{\code{listPorts()}}{list all available ports on the system}
#' }
#' 
#' @examples
#' # for this example I used the 'null-modem' emulator 'com0com' for Windows
#' # which is available on 'http://com0com.sourceforge.net/'
#' # Here the pair of com-ports is 'CNCA0' <-> 'CNCB0'
#' 
#' # Test the functionality:
#' # ======================
#' #
#' # first: install the virtual null-modem connection like
#' #        com0com (win) or tty0tty (linux)
#' #        Hint: Some unix insist on port names like 'ttyS[n]'.
#' # 
#' # second: setup a terminal program (like HTerm or gtkterm) and listen to 
#' #         com-port 'CNCB0' (or what ever you have installed)
#' #         or (for unix only) 'cat /dev/tnt1' will output tnt1 to console
#' 
#' \dontrun{
#' 
#' # Now configure one of the com-ports with appropriate connection properties
#' con <- serialConnection(name = "testcon",port = "CNCA0"
#'                        ,mode = "115200,n,8,1"
#'                        ,newline = 1
#'                        ,translation = "crlf"
#'                        )
#' 
#' # let's open the serial interface
#' 
#' open(con)
#' 
#' # write some stuff
#' write.serialConnection(con,"Hello World!")
#' 
#' # read, in case something came in
#' read.serialConnection(con)
#' 
#' # close the connection
#' close(con)
#' }
#' @concept RS232 serial
#' @concept RS485 serial
#' @concept RS422 serial
#' @concept USB serial
#' 
#' @docType package
#' @name serial
NULL