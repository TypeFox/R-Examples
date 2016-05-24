summary.Game <-
function(object, ...) {
	#assignInNamespace("summary.Game", summary.Game, ns = asNamespace("base"))
   x<-object
   cat("\n")
   cat("Characteristic form of the game","\n")
   cat("\n")
   cat("Number of agents:",x$Ag,"\n")
   cat("\n")
   cat("Coaliton Value(s)","\n")
   cat("\n")
   print(x$Lex)
   cat("\n")
}
