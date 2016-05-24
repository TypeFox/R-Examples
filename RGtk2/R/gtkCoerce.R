as.GtkStockItem <-
function(x)
{
    x <- as.struct(x, "GtkStockItem", c("stock_id", "label", "modifier", "keyval", "translation_domain"))
    x[[1]] <- as.character(x[[1]])
    x[[2]] <- as.character(x[[2]])
    x[[4]] <- as.numeric(x[[4]])
    x[[5]] <- as.character(x[[5]])

    return(x)
}

as.GtkActionEntry <-
function(x)
{
    x <- as.struct(x, "GtkActionEntry", c("name", "stock_id", "label", "accelerator", "tooltip", "callback"))
    x[[1]] <- as.character(x[[1]])
    x[[2]] <- as.character(x[[2]])
    x[[3]] <- as.character(x[[3]])
    x[[4]] <- as.character(x[[4]])
    x[[5]] <- as.character(x[[5]])
    if (!is.null(x[[6]]))
        x[[6]] <- as.function(x[[6]])

    return(x)
}
as.GtkToggleActionEntry <-
function(x)
{
    x <- as.struct(x, "GtkToggleActionEntry", c("name", "stock_id", "label", "accelerator", "tooltip", "callback", "is_active"))
    x[[1]] <- as.character(x[[1]])
    x[[2]] <- as.character(x[[2]])
    x[[3]] <- as.character(x[[3]])
    x[[4]] <- as.character(x[[4]])
    x[[5]] <- as.character(x[[5]])
    if (!is.null(x[[6]]))
        x[[6]] <- as.function(x[[6]])
    x[[7]] <- as.logical(x[[7]])

    return(x)
}
as.GtkRadioActionEntry <-
function(x)
{
    x <- as.struct(x, "GtkActionEntry", c("name", "stock_id", "label", "accelerator", "tooltip", "value"))
    x[[1]] <- as.character(x[[1]])
    x[[2]] <- as.character(x[[2]])
    x[[3]] <- as.character(x[[3]])
    x[[4]] <- as.character(x[[4]])
    x[[5]] <- as.character(x[[5]])
    x[[6]] <- as.integer(x[[6]])

    return(x)
}
as.GtkFileFilterInfo <-
function(x)
{
	x <- as.struct(x, "GtkFileFilterInfo", c("contains", "filename", "uri", "display.name", "mime.type"))
  
	x[[2]] <- as.character(x[[2]])
	x[[3]] <- as.character(x[[3]])
	x[[4]] <- as.character(x[[4]])
	x[[5]] <- as.character(x[[5]])
	
	return(x)
}
as.GtkSettingsValue <-
function(x)
{
	x <- as.struct(x, "GtkSettingsValue", c("origin", "value"))
	
	x[[1]] <- as.character(x[[1]])
	
	return(x)
}
as.GtkItemFactoryEntry <-
function(x)
{
	x <- as.struct(x, "GtkItemFactoryEntry", c("path", "accelerator", "callback", "callback.action", "item.type", "extra.data"))
	
	x[[1]] <- as.character(x[[1]])
	x[[2]] <- as.character(x[[2]])
	x[[3]] <- as.function(x[[3]])
	x[[4]] <- as.numeric(x[[4]])
	x[[5]] <- as.character(x[[5]])
	
	return(x)
}
as.GtkAllocation <-
function(x)
{
	x <- as.GdkRectangle(x)
	class(x) <- "GtkAllocation"
	return(x)
}
as.GtkTargetEntry <-
function(x)
{
	x <- as.struct(x, "GtkTargetEntry", c("target", "flags", "info"))
	x[[1]] <- as.character(x[[1]])
	x[[3]] <- as.integer(x[[3]])
	return(x)
}
as.GtkRecentFilterInfo <-
function(x)
{
  x <- as.struct("GtkRecentFilterInfo", c("contains", "uri", "display_name", "mime_type", "applications", "groups", "age"))
  
  x[[2]] <- as.character(x[[2]])
  x[[3]] <- as.character(x[[3]])
  x[[4]] <- as.character(x[[4]])
  x[[5]] <- as.list(as.character(x[[5]]))
  x[[6]] <- as.list(as.character(x[[6]]))
  x[[7]] <- as.integer(x[[7]])
  
  return(x)
}
as.GtkRecentData <-
function(x)
{
  x <- as.struct("GtkRecentData", c("display_name", "description", "mime_type", "app_name", "app_exec", "groups", "is_private"))
  x[[1]] <- as.character(x[[1]])
  x[[2]] <- as.character(x[[2]])
  x[[3]] <- as.character(x[[3]])
  x[[4]] <- as.character(x[[4]])
  x[[5]] <- as.character(x[[5]])
  x[[6]] <- as.list(as.character(x[[6]]))
  x[[7]] <- as.logical(x[[7]])
  
  return(x)
}
as.GtkPageRange <-
function(x)
{
  x <- as.struct("GtkPageRange", c("start", "end"))
  x[[1]] <- as.integer(x[[1]])
  x[[2]] <- as.integer(x[[2]])
  
  return(x)
}
