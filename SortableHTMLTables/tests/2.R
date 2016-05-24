# Test rendering of factors as character strings.

library('SortableHTMLTables')
library('testthat')

sortable.html.table(iris, 'iris.html', 'iris')

text <- paste(readLines(file.path('iris', 'iris.html'), n = -1), collapse = '\n')

expect_that(text, matches('setosa'))
expect_that(text, matches('versicolor'))
expect_that(text, matches('virginica'))

unlink('iris', recursive = TRUE)
