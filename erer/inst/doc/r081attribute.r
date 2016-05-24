# Run the program version for Sun et al. (2007) by source().
# Create a new object for demonstration.
# Show 3 selected observations.
source("C:/aErer/r072sunSJAF.r", echo = FALSE); getwd(); dir()
nam <- daInsNam; nam[c(3, 11, 14), ]

# Show all attributes or a specific one
attributes(x = nam)
attr(x = nam, which = "class")

# Revise or create attributes
attr(x = nam, which = "names") <- c("vari", "defi")  # Revise an old attr
attr(x = nam, which = "years") <- 2005:2008          # New numerical attr
attr(x = nam, which = "color") <- "red"              # New character attr
attr(x = daIns, which = "description") <- daInsNam   # New data.frame attr
attributes(x = nam)  # Show all attributes again
nam[c(3, 11, 14), ]  # Show 3 selected obs again

# Reveal the attribute of an object from regressions
class(ra); mode(ra); typeof(ra); length(ra); tsp(ra); names(ra)[1:4]

# Change the class for an existing object
(class(ra) <- c("glm", "lm", "myClass"))

# Methods for generic function or class
methods(generic.function = mean); methods(class = "maTrend")