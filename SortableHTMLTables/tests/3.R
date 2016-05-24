# Test arbitrary page titles.

library('SortableHTMLTables')
library('testthat')

sortable.html.table(iris,
                    output.file = 'iris.html',
                    output.directory = 'iris',
                    page.title = 'Iris Data Set')

text <- paste(readLines(file.path('iris', 'iris.html'), n = -1), collapse = '\n')

expect_that(text, matches('Iris Data Set'))

unlink('iris', recursive = TRUE)
