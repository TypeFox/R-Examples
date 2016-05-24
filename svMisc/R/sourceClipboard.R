sourceClipboard <- function (primary = TRUE, ...)
{
	## Source data from the clipboard, manage clipboard correctly depending
	## on the OS
	if (isWin()) { # Windows OS
		data <- file("clipboard")
	} else if (isMac()) { # Mac OS
		data <- pipe("pbpaste")
	} else { # Must be Linux/Unix
		if (primary) {
			data <- file("X11_clipboard")
		} else {
			data <- file("X11_secondary")
		}
	}
	on.exit(close(data))
	## Invoke source() with the data from the clipboard
	invisible(source(data, ...))
}
