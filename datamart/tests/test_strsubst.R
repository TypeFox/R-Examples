library(datamart)

test_strsubst <- function() {
  verbose <- FALSE
  
  # one string
  print(strsubst("$(who) likes $(what)!", list(who="tim", what="kung pao"), verbose=verbose))
  print(strsubst("$(who) likes $", list(who="tim", what="kung pao"), verbose=verbose))
  print(strsubst("$(unresolved) likes $(what)", list(who="tim", what="kung pao"), verbose=verbose))
  print(strsubst("$(who) likes $(what)", list(who="tim", what="$10000"), verbose=verbose))
  print(strsubst("$(who) likes $$(what)", list(who="tim", what=10000), verbose=verbose))

  # vector of strings
  print(
    strsubst(
      c("$(who) is human.", "All humans are mortal.", "Therefore, $(who) is mortal.", "$(who) knew $(who) knew nothing."), 
      list(who="Sokrates"), verbose=verbose
    )
  )
  
}

test_strsubst()

