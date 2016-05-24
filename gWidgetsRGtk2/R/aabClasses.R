### these classes need to be defined before their subclasses. Alphabetical doesn't cut
### is so they go here.

### this must come after aaaGenerics, as there gComponentRGtk is defined
setClass("gEditRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )
setClass("gGroupRGtk",
         contains="gContainerRGtk",
         prototype=prototype(new("gContainerRGtk"))
         )

setClass("gNotebookRGtk",
         contains="gComponentRGtk"
         )
