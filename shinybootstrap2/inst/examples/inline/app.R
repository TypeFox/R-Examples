appExpr <- quote({
  shinyApp(
    ui = fluidPage(
      checkboxGroupInput("variable", "Variable:", inline = TRUE,
                         c("Cylinders" = "cyl",
                           "Transmission" = "am",
                           "Gears" = "gear",
                           "Horsepower" = "hp",
                           "Weight" = "wt")),
      radioButtons("radio", "Variable:", inline = TRUE,
                         c("Cylinders" = "cyl",
                           "Transmission" = "am",
                           "Gears" = "gear",
                           "Horsepower" = "hp",
                           "Weight" = "wt"))
    ),
    server = function(input, output) {}
  )
})

# The useBS2 environment var controls whether or not we run the app using the
# wihtBootstrap2 wrapper function.
if (is.null(.GlobalEnv$useBS2))
  .GlobalEnv$useBS2 <- FALSE

if (.GlobalEnv$useBS2) {
  shinybootstrap2::withBootstrap2(appExpr, quoted = TRUE)
} else {
  eval(appExpr)
}
