## CLI torture test, copyright (c) Ph. Grosjean (phgrosjean@sciviews.org)
## GNU GPL => 2 license
## A series of commands to check for R CLI (or console widget)
## Version 1.0-0

## Simple instructions
1+1            # Simple command with one line of output
1:100          # Simple command with several lines of output
search()       # Evaluation of a function
ls()           # Idem... look if this function is evaluated in .GlobalEnv!

## Multiple instructions on one line
1+1; 2+2    # Two successive prints
1+1; cat("test\n"); 2+2   # Print, cat, print
1+1; cat("test"); 2+2   # Idem, but cat without \n

## Visible/invisible output
invisible(1)   # Command with invisible result
a <- 1:10      # Idem
(a <- 1:10)    # Idem, with visible output
for (i in 1:3) print(i)  # Output within a loop
for (i in 1:3) i  # Do not explicit use print(), so output nothing

cat("test\n")  # Simple text print with carriage return
cat("test")    # Idem, but without carriage return

## S4 objects
setClass("track", representation(x = "numeric", y = "numeric"))
setMethod("show", "track", function(object) {cat("track:\n"); print(object@x); print(object@y)})
tr <- new("track", x = 1:3, y = 4:6)    # invisible
(tr <- new("track", x = 1:3, y = 4:6))  # visible
tr             # show it
show(tr)       # idem
print(tr)      # This is the same!

## Special characters
cat("\ttabs\n")  # \t (tabulation) should indent the text by 4 characters
cat("tesg\bt\n") # \b (backspace) erases previous character thus it print "test"
alarm()        # idem as cat("\a"), should sound a bip

## Code parsing and multiline commands
log)           # Generate and error
1+1; log)      # Should run first instruction before generating the error
## This is a partial line without comments (should issue the continue prompt)
log(
10) +
1

log(          # This is partial with comments
10)

search(); log(
10)

log(          # Issuing an empty line between instructions

10)

## String on several lines
text <- "some
text"
text         # should print "some\ntext"

## Idem, but with single quote
text <- 'some
text'
text         # should print "some\ntext"

## A horrible code with a variable name on two lines (but correct syntax)!
`var
1` <- 1
`var\n1`

## Truncation of very long output
options(max.print = 1000)   # Put a lower limit
1:1100       # Should be truncated
options(max.print = NULL)   # Restore default value

## Errors messages
nonExistingVar
cos("a")
cos(nonExisting)
stop("Error!"); 1:3                  # Error in .GlobalEnv (no call info)
(function() stop("Error!"))(); 1:3   # Error; second command not evaluated

## Warnings handling (a little tricky!)
options(warn = -1)      # Do not display warnings
{warning("Warn!"); 1:3} # Simple warning
{warning("Warn!", immediate. = TRUE); 1:3}  # Should issue the warning anyway

options(warn = 0)       # Delay warning display
{warning("Warn!"); 1:3} # Simple delayed warning
{warning("Warn!", immediate. = TRUE); 1:3} # Do not delay warning
for (i in 1:3) {print(i); warning("test")}
for (i in 1:4) {print(i); warning("test", immediate. = (i < 3))}

options(warn = 1)       # Display warnings immediatelly
{warning("Warn!"); 1:3}

options(warn = 2)       # Every warning generates an error
{warning("Warn!"); 1:3} # Warning turned into an error
warnings()


## Warnings inside functions
options(warn = -1)
(function() {warning("Warn!"); 1:3})()
options(warn = 0)
(function() {warning("Warn!"); 1:3})()
options(warn = 1)
(function() {warning("Warn!"); 1:3})()
options(warn = 2)
(function() {warning("Warn!"); 1:3})()

## Multiple warnings and/or errors (warn = 0 + 9, 10, 11, 49, 50 & 60 warnings)
options(warn = 0)
for (i in 1:9) warning("Warn ", i)
warnings()      # Redisplay last warnings
for (i in 1:10) warning("Warn ", i)
warnings()
for (i in 1:11) warning("Warn ", i)
warnings()
for (i in 1:49) warning("Warn ", i)
warnings()
for (i in 1:50) warning("Warn ", i)
warnings()
for (i in 1:60) warning("Warn ", i)
warnings()

## warning() and then, error message with warn = 0 ("In addition: ...")
options(warn = 0)
(function() {warning("Warn!"); stop("Error!"); 1:3})()
options(warn = 1)
(function() {warning("Warn!"); stop("Error!"); 1:3})()

## Messages handling
{message("A message"); 1:3} # Should issue the message
simpleMessage("test")
simpleMessage("test", call = "me")

## Multiline or very long message
options(warning.length = 100)
warning("A very long message for my warning that should be truncated or at least flowed on several lines, what does it gives here?")
warning("A multiline warning\nSecond line,\nThird line")
