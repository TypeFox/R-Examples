function(){}                      # empty
function(){cat("Hello World!\\n")}  # with body
function()cat("Hello World!\\n")    # single line function
function()
  cat("Hello World!\\n")            # single line on other line
a <- function(){}                  # assignment

# Arguments
a <- function(file){cat("Hello World!\\n", file=file)}
a <- function(file, b=2){cat("Hello World!\\n", file=file)}

# Multiline arguments
a <- function(file,
              b=2,
              ...) {
                cat("Hello World!\\n", file=file)
              }

# embeded
llply(1:9, function(x){x ** 0.5})
