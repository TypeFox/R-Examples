pause <- function() invisible(readline())

print("Summary information about emotions")
summary(emotions)
pause()

print("List of attributes in emotions")
emotions$attributes
pause()

print("List of general measures")
emotions$measures
pause()

print("Information about labels in emotions")
emotions$labels
pause()

print("Information about labelsets in emotions")
emotions$labelsets
pause()

print("Dataset content")
emotions$dataset

