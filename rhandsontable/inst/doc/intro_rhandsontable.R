## ----echo = FALSE--------------------------------------------------------
library(rhandsontable)
library(knitr)

opts_knit$set(warning = FALSE, error = FALSE, message = FALSE, cache = FALSE,
              fig.width=7, fig.height=3)

## ------------------------------------------------------------------------
DF = data.frame(integer = 1:10,
                   numeric = rnorm(10),
                   logical = rep(TRUE, 10), 
                   character = LETTERS[1:10],
                   factor = factor(letters[1:10], levels = letters[10:1], 
                                   ordered = TRUE),
                   factor_allow = factor(letters[1:10], levels = letters[10:1], 
                                         ordered = TRUE),
                   date = seq(from = Sys.Date(), by = "days", length.out = 10),
                   stringsAsFactors = FALSE)

rhandsontable(DF, width = 600, height = 300) %>%
  hot_col("factor_allow", allowInvalid = TRUE)

## ------------------------------------------------------------------------
DF_na = data.frame(integer = c(NA, 2:10), 
                   logical = c(NA, rep(TRUE, 9)), 
                   character = c(NA, LETTERS[1:9]),
                   factor = c(NA, factor(letters[1:9])),
                   date = c(NA, seq(from = Sys.Date(), by = "days", 
                                    length.out = 9)),
                   stringsAsFactors = FALSE)

DF_na$factor_ch = as.character(DF_na$factor)
DF_na$date_ch = c(NA, as.character(seq(from = Sys.Date(), by = "days", 
                                       length.out = 9)))

rhandsontable(DF_na, width = 550, height = 300)

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

# try updating big to a value not in the dropdown
rhandsontable(DF, rowHeaders = NULL, width = 550, height = 300) %>%
  hot_col(col = "big", type = "dropdown", source = LETTERS) %>%
  hot_col(col = "small", type = "autocomplete", source = letters,
          strict = FALSE)

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, width = 550, height = 300) %>%
  hot_col("small", "password")

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

DF$chart = c(sapply(1:5,
                    function(x) jsonlite::toJSON(list(values=rnorm(10),
                                                      options = list(type = "bar")))),
             sapply(1:5,
                    function(x) jsonlite::toJSON(list(values=rnorm(10),
                                                      options = list(type = "line")))))

rhandsontable(DF, rowHeaders = NULL, width = 550, height = 300) %>%
  hot_col("chart", renderer = htmlwidgets::JS("renderSparkline"))

## ----fig.height = 5, fig.width = 8---------------------------------------
DF = data.frame(
  title = c(
    "<a href='http://www.amazon.com/Professional-JavaScript-Developers-Nicholas-Zakas/dp/1118026691'>Professional JavaScript for Web Developers</a>",
    "<a href='http://shop.oreilly.com/product/9780596517748.do'>JavaScript: The Good Parts</a>",
    "<a href='http://shop.oreilly.com/product/9780596805531.do'>JavaScript: The Definitive Guide</a>"
  ),
  desc = c(
    "This <a href='http://bit.ly/sM1bDf'>book</a> provides a developer-level introduction along with more advanced and useful features of <b>JavaScript</b>.",
    "This book provides a developer-level introduction along with <b>more advanced</b> and useful features of JavaScript.",
    "<em>JavaScript: The Definitive Guide</em> provides a thorough description of the core <b>JavaScript</b> language and both the legacy and standard DOMs implemented in web browsers."
  ),
  comments = c(
    "I would rate it &#x2605;&#x2605;&#x2605;&#x2605;&#x2606;",
    "This is the book about JavaScript",
    "I've never actually read it, but the <a href='http://shop.oreilly.com/product/9780596805531.do'>comments</a> are highly <strong>positive</strong>."
  ), 
  cover = c(
    "http://ecx.images-amazon.com/images/I/51bRhyVTVGL._SL50_.jpg",
    "http://ecx.images-amazon.com/images/I/51gdVAEfPUL._SL50_.jpg",
    "http://ecx.images-amazon.com/images/I/51VFNL4T7kL._SL50_.jpg"
 ),
 stringsAsFactors = FALSE
)

rhandsontable(DF, allowedTags = "<em><b><strong><a><big>", 
              width = 800, height = 450, rowHeaders = FALSE) %>%
  hot_cols(colWidths = c(200, 200, 200, 80)) %>%
  hot_col(1:2, renderer = "html") %>%
  hot_col(1:3, renderer = htmlwidgets::JS("safeHtmlRenderer")) %>%
  hot_col(4, renderer = "
    function(instance, td, row, col, prop, value, cellProperties) {
      var escaped = Handsontable.helper.stringify(value),
        img;
  
      if (escaped.indexOf('http') === 0) {
        img = document.createElement('IMG');
        img.src = value;
  
        Handsontable.Dom.addEvent(img, 'mousedown', function (e){
          e.preventDefault(); // prevent selection quirk
        });
  
        Handsontable.Dom.empty(td);
        td.appendChild(img);
      }
      else {
        // render as text
        Handsontable.renderers.TextRenderer.apply(this, arguments);
      }
  
      return td;
    }")

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, width = 550, height = 300) %>%
  hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)

## ------------------------------------------------------------------------
MAT = matrix(rnorm(50), nrow = 10, dimnames = list(LETTERS[1:10],
                                                   letters[1:5]))

rhandsontable(MAT, width = 550, height = 300) %>%
  hot_context_menu(
    customOpts = list(
      csv = list(name = "Download to CSV",
                    callback = htmlwidgets::JS(
                      "function (key, options) {
                         var csv = csvString(instance);

                         var link = document.createElement('a');
                         link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                           encodeURIComponent(csv));
                         link.setAttribute('download', filename);

                         document.body.appendChild(link);
                         link.click();
                         document.body.removeChild(link);
                       }"))))

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10,
                bool = TRUE,
                big = LETTERS[1:10],
                small = factor(letters[1:10]),
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, search = TRUE, width = 550, height = 300) %>%
  hot_context_menu(
    customOpts = list(
      search = list(name = "Search",
                    callback = htmlwidgets::JS(
                      "function (key, options) {
                         var srch = prompt('Search criteria');

                         this.search.query(srch);
                         this.render();
                       }"))))

## ------------------------------------------------------------------------
DF = data.frame(int = 1:10, float = rnorm(10), cur = rnorm(10) * 1E5,
                lrg = rnorm(10) * 1E8, pct = rnorm(10))

rhandsontable(DF, width = 550, height = 300) %>%
  hot_col("float", format = "0.0") %>%
  hot_col("cur", format = "$0,0.00") %>%
  hot_col("lrg", format = "0a") %>%
  hot_col("pct", format = "0%")

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, readOnly = TRUE, width = 550, height = 300) %>%
  hot_col("val", readOnly = FALSE)

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, width = 550, height = 300) %>%
  hot_cols(columnSorting = TRUE)

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

# click on a cell to see the highlighting
rhandsontable(DF, width = 550, height = 300) %>%
  hot_table(highlightCol = TRUE, highlightRow = TRUE)

## ----fig.height = 6, fig.width = 6---------------------------------------
MAT = matrix(rnorm(50), nrow = 10, dimnames = list(LETTERS[1:10],
                                                   letters[1:5]))

rhandsontable(MAT, width = 600, height = 600) %>%
  hot_cols(colWidths = 100) %>%
  hot_rows(rowHeights = 50)

## ------------------------------------------------------------------------
MAT = matrix(rnorm(26 * 26), nrow = 26, dimnames = list(LETTERS, letters))

# scroll through the table to see the fixed row and column
rhandsontable(MAT, width = 550, height = 300) %>%
  hot_cols(fixedColumnsLeft = 1) %>%
  hot_rows(fixedRowsTop = 1)

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, width = 550, height = 300) %>%
  hot_cell(1, 1, "Test comment")

## ------------------------------------------------------------------------
MAT_comments = matrix(ncol = ncol(DF), nrow = nrow(DF))
MAT_comments[1, 1] = "Test comment"
MAT_comments[2, 2] = "Another test comment"

rhandsontable(DF, comments = MAT_comments, width = 550, height = 300)

## ------------------------------------------------------------------------
MAT = matrix(rnorm(50), nrow = 10, dimnames = list(LETTERS[1:10],
                                                   letters[1:5]))

rhandsontable(MAT, width = 550, height = 300) %>%
  hot_table(customBorders = list(list(
    range = list(from = list(row = 1, col = 1),
                 to = list(row = 2, col = 2)),
    top = list(width = 2, color = "red"),
    left = list(width = 2, color = "red"),
    bottom = list(width = 2, color = "red"),
    right = list(width = 2, color = "red"))))

## ------------------------------------------------------------------------
MAT = matrix(rnorm(50), nrow = 10, dimnames = list(LETTERS[1:10],
                                                   letters[1:5]))

rhandsontable(MAT * 10, width = 550, height = 300) %>%
  hot_validate_numeric(col = 1, min = -50, max = 50, exclude = 40)

rhandsontable(MAT * 10, width = 550, height = 300) %>%
  hot_validate_numeric(col = 1, choices = c(10, 20, 40))

## ------------------------------------------------------------------------
DF = data.frame(val = 1:10, bool = TRUE, big = LETTERS[1:10],
                small = letters[1:10],
                dt = seq(from = Sys.Date(), by = "days", length.out = 10),
                stringsAsFactors = FALSE)

rhandsontable(DF, width = 550, height = 300) %>%
  hot_validate_character(col = "big", choices = LETTERS[1:10])

## ------------------------------------------------------------------------
MAT = matrix(rnorm(50), nrow = 10, dimnames = list(LETTERS[1:10],
                                                   letters[1:5]))

# try to update any cell to 0
rhandsontable(MAT * 10, width = 550, height = 300) %>%
  hot_cols(validator = "
           function (value, callback) {
            setTimeout(function(){
              callback(value != 0);
            }, 1000)
           }",
           allowInvalid = FALSE)

## ---- fig.width = 8------------------------------------------------------
MAT = matrix(runif(100, -1, 1), nrow = 10,
             dimnames = list(LETTERS[1:10], LETTERS[1:10]))
diag(MAT) = 1
MAT[upper.tri(MAT)] = MAT[lower.tri(MAT)]
rhandsontable(MAT, readOnly = TRUE, width = 750, height = 300) %>%
  hot_cols(renderer = "
           function (instance, td, row, col, prop, value, cellProperties) {
             Handsontable.renderers.TextRenderer.apply(this, arguments);
             if (row == col) {
              td.style.background = 'lightgrey';
             } else if (col > row) {
              td.style.background = 'grey';
              td.style.color = 'grey';
             } else if (value < -0.75) {
              td.style.background = 'pink';
             } else if (value > 0.75) {
              td.style.background = 'lightgreen';
             }
           }")

## ------------------------------------------------------------------------
MAT = matrix(rnorm(50), nrow = 10, dimnames = list(LETTERS[1:10],
                                                   letters[1:5]))

rhandsontable(MAT, width = 550, height = 300) %>%
  hot_heatmap()

