sortable.html.table <- function(df,
                                output.file = 'output.html',
                                output.directory = getwd(),
                                page.title = 'Untitled Page')
{
  # Create output.directory if it does not exist already.
  if (! file.exists(output.directory))
  {
    dir.create(output.directory, recursive = TRUE)
  }
  
  # Print an HTML file from our template.
  brew(file = system.file('template.brew', package = 'SortableHTMLTables'),
       output = file.path(output.directory, output.file))
  
  # Copy Javascript, CSS and GIF assets.
  assets <- dir(system.file('assets', package = 'SortableHTMLTables'))
  
  for (asset in assets)
  {
    file.copy(file.path(system.file('assets', package = 'SortableHTMLTables'), asset),
              file.path(output.directory, asset))
  }
}
