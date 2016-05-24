searchEntry <-
setRefClass("RzSearchEntry",
  fields = c("entry.search"),
  methods = list(
    initialize    = function(...) {
      initFields(...)
      entry.search  <<- gtkEntryNew()
      entry.search$modifyFont(pangoFontDescriptionFromString(rzSettings$getGlobalFont()))
      entry.search$setWidthChars(20)
      entry.search$setIconFromStock(GtkEntryIconPosition["primary"], "gtk-find")
      entry.search$setTooltipText(gettext("Prev: Up Arrow Key\nNext: Down Arrow Key"))
    },
    searchFunc = function(model, column, key, iter) {
        key <-  gsub("([]\\[\\^?|(){}!$&*+?])", "\\\\\\1", localize(key))
        value <- model$get(iter, column.definition["vars"],  column.definition["var.labs"])
        match <- tryCatch(
          grepl(key, c(localize(value[[1]][1]), localize(value[[2]][1])), perl=TRUE, ignore.case=TRUE),
          error=function(...) FALSE)
        return(!any(match))
    }
))
searchEntry$accessors("entry.search")
