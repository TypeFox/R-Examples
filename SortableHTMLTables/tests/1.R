# Test basic functionality.

library('SortableHTMLTables')
library('testthat')

df <- data.frame(X = rnorm(10), Y = runif(10), Z = rcauchy(10))

sortable.html.table(df, 'sample.html', 'sandbox')

expect_that(file.exists('sandbox'), is_true())

assets <- c('asc.gif',
            'bg.gif',
            'desc.gif',
            'jquery-1.4.2.js',
            'jquery.tablesorter.js',
            'style.css')

for (asset in assets)
{
  expect_that(file.exists(file.path('sandbox', asset)), is_true())
}

expect_that(file.exists(file.path('sandbox', 'sample.html')), is_true())

unlink('sandbox', recursive = TRUE)
