#' Browse photo queries in your default browser.
#' 
#' @import whisker
#' @importFrom assertthat assert_that
#' @param input Input, usually output from a call to \code{\link[ecoengine]{ee_photos}}
#' @param output Path and file name for output file. If NULL, a temp file is used.
#' @param browse Browse file in your default browse immediately after file creation.
#'    If \code{FALSE}, the file is written, but not opened.
#' @export
#' @examples 
#' # view_photos(ee_photos())
#' # Pictures of racoons
#' # view_photos(ee_photos(scientific_name = "Procyon lotor", quiet = TRUE))
#' # or the California Condor
#' # view_photos(ee_photos(scientific_name = "Gymnogyps californianus", quiet = TRUE))

view_photos <- function(input = NULL, output = NULL, browse = TRUE)
{
  if(is.null(input))
    stop("Please supply some input")

assert_that(identical(input$type, "photos"))
 photo_list <- apply(input$data, 1, function(x) as.list(x))
 photo_list <- unname(photo_list)
  template <-
    '<!DOCTYPE html>
      <head>
        <meta charset="utf-8">
              <title>ecoengine - view highlighs</title>
              <meta name="viewport" content="width=device-width, initial-scale=1.0">
              <meta name="description" content="View highlights from an ecoengine photo search">
              <meta name="author" content="ecoengine">

              <!-- Le styles -->
              <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.min.css" rel="stylesheet">
              <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.css" rel="stylesheet">
              <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap-responsive.css" rel="stylesheet">
              <script src="http://use.edgefonts.net/quattrocento-sans.js"></script>
              <style>
              body {
                margin: 0;
                font-family: quattrocento-sans, Helvetica, Arial, sans-serif;
                font-size: 16px;
                line-height: 20px;
                color: #333333;
                background-color: #ffffff;
              }
              </style>
      </head>

      <body>

      <div class="container">

      <center><h2>Ecoengine Photo Viewer</h2></center>

      <table class="table table-striped table-hover" align="center">
              <thead>
                      <tr>
                              <th>Photo</th>
                              <th>Authors</th>
                              <th>Locality / County</th>
                              <th>Notes</th>
                              <th>Start Date</th>
                      </tr>
              </thead>
              <tbody>
        {{#photo_list}}
          <tr><td>
          <a href="{{remote_resource}}"><img src="{{media_url}}" height = 250></a></td>
          <td>{{authors}}</td>
          <td>{{locality}}, {{county}}</td>
          <td>{{photog_notes}}</td>
          <td>{{begin_date}}</td>
          </tr>
        {{/photo_list}}
        </tbody>
      </table>
      </div>

      <script src="http://code.jquery.com/jquery-2.0.3.min.js"></script>
      <script src="http://netdna.bootstrapcdn.com/bootstrap/3.0.2/js/bootstrap.min.js"></script>

      </body>
      </html>'
        
  rendered <- whisker.render(template)
  rendered <- gsub("&lt;em&gt;", "<b>", rendered)
  rendered <- gsub("&lt;/em&gt;", "</b>", rendered)
  if(is.null(output))
    output <- tempfile(fileext=".html")
  write(rendered, file = output)
  if(browse) browseURL(output)
}