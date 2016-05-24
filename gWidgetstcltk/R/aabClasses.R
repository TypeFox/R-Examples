### these classes need to be defined before their subclasses. Alphabetical doesn't cut
### is so they go here.

## for coerce.with
setClassUnion("NULLorFunction",c("NULL","function"))

### this must come after aaaGenerics, as there gComponenttcltk is defined
setClass("gEdittcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )
setClass("gGrouptcltk",
         representation = representation("gContainertcltk",
           horizontal="logical"),
         contains="gContainertcltk",
         prototype=prototype(new("gContainertcltk"))
         )

setClass("gWindowtcltk",
         representation = representation("gContainertcltk"),
         ##           horizontal="logical"),
         contains="gContainertcltk",
         prototype=prototype(new("gContainertcltk"))
         )

setClass("gNotebooktcltk",
         representation = representation("gContainertcltk",
           closebuttons="logical",
           dontCloseThese="numeric"),
         contains="gContainertcltk"
         )
