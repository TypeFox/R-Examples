# Setting up an environment for use within gridSVG
.gridSVGEnv <- new.env()

# Don't show any progress bars or messages by default
assign("showProgress", FALSE, envir = .gridSVGEnv)

# Setting a context level for clipGrobs, clipping paths and masks.
# Allows popping viewports to work correctly if we know many
# SVG groups we need to "pop".
assign("contextLevels", 0, envir = .gridSVGEnv)
assign("contextNames", character(0), envir = .gridSVGEnv)

# The following definitions are required for references to work
# because we cannot assume that everything is determined at
# draw time.
# References require definitions to be created *prior* to drawing
# which means we need naming systems and usage tables present.

# Needed by getID()
assign("uniqueNames", TRUE, envir = .gridSVGEnv)
assign("prefix", "", envir = .gridSVGEnv)

# We are going to be filling this list with lists of information (describing
# reference definitions) keyed by their labels
assign("refDefinitions", list(), envir = .gridSVGEnv)
assign("usageTable", 
       data.frame(name = character(0),
                  suffix = integer(0),
                  type = character(0),
                  selector = character(0),
                  xpath = character(0),
                  stringsAsFactors = FALSE),
       envir = .gridSVGEnv)
assign("pchUsageTable",
       matrix(c(0:127, logical(128)), ncol = 2,
              dimnames = list(NULL, c("pch", "used"))),
       envir = .gridSVGEnv)
assign("refUsageTable",
       data.frame(label = character(0),
                  used = logical(0),
                  stringsAsFactors = FALSE),
       envir = .gridSVGEnv)
