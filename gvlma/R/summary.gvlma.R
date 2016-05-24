"summary.gvlma" <-
function(object, ...)
{
  print(NextMethod("summary", object))
  display.gvlmatests(object)
}

