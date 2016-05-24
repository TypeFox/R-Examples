appExpr <- quote({
  function(input, output, clientData, session) {

    observe({
      # We'll use the input$controller variable multiple times, so save it as x
      # for convenience.
      x <- input$controller

      # TRUE if x is even, FALSE otherwise. Will be used below
      x_even <- x %% 2 == 0


      # Text ============================================================
      # This will change the value of input$inText, based on x
      updateTextInput(session, "inText", value = paste("New text", x))

      # Can also set the label, this time for input$inText2
      updateTextInput(session, "inText2",
        label = paste("New label", x),
        value = paste("New text", x))


      # Number ==========================================================
      # Number inputs have more settable parameters
      updateNumericInput(session, "inNumber", value = x)

      updateNumericInput(session, "inNumber2",
        label = paste("Number label ", x),
        value = x, min = x-10, max = x+10, step = 5)


      # Slider input ====================================================
      # Similar to number and text. only label and value can be set for slider
      updateSliderInput(session, "inSlider",
        label = paste("Slider label", x),
        value = x)


      # Slider range input ==============================================
      # For sliders that pick out a range, pass in a vector of 2 values.
      updateSliderInput(session, "inSlider2", value = c(x-1, x+1))

      # An NA means to not change that value (the low or high one)
      updateSliderInput(session, "inSlider3", value = c(NA, x+2))


      # Date input ====================================================
      # Similar to number and text. only label and value can be set for slider
      updateDateInput(session, "inDate",
        label = paste("Date label", x),
        value = paste("2013-04-", x, sep=""))


      # Date range input ================================================
      # Similar to number and text. only label and value can be set for slider
      updateDateRangeInput(session, "inDateRange",
        label = paste("Date range label", x),
        start = paste("2013-01-", x, sep=""),
        end = paste("2013-12-", x, sep=""),
        min = paste("2001-01-", x, sep=""),
        max = paste("2030-12-", x, sep="")
      )

      # Checkbox ========================================================
      updateCheckboxInput(session, "inCheckbox", value = x_even)


      # Checkbox group ==================================================
      # Create a list of new options, where the name of the items is something
      # like 'option label x 1', and the values are 'option-x-1'.
      cb_options <- list()
      cb_options[[sprintf("option label %d 1", x)]] <- sprintf("option-%d-1", x)
      cb_options[[sprintf("option label %d 2", x)]] <- sprintf("option-%d-2", x)

      # Change values for input$inCheckboxGroup
      updateCheckboxGroupInput(session, "inCheckboxGroup", choices = cb_options)

      # Can also set the label and select items
      updateCheckboxGroupInput(session, "inCheckboxGroup2",
        label = paste("checkboxgroup label", x),
        choices = cb_options,
        selected = sprintf("option-%d-2", x)
      )

      # Radio group =====================================================
      # Create a list of new options, where the name of the items is something
      # like 'option label x 1', and the values are 'option-x-1'.
      r_options <- list()
      r_options[[sprintf("option label %d 1", x)]] <- sprintf("option-%d-1", x)
      r_options[[sprintf("option label %d 2", x)]] <- sprintf("option-%d-2", x)

      # Change values for input$inRadio
      updateRadioButtons(session, "inRadio", choices = r_options)

      # Can also set the label and select an item
      updateRadioButtons(session, "inRadio2",
        label = paste("Radio label", x),
        choices = r_options,
        selected = sprintf("option-%d-2", x)
      )


      # Select input ====================================================
      # Create a list of new options, where the name of the items is something
      # like 'option label x 1', and the values are 'option-x-1'.
      s_options <- list()
      s_options[[sprintf("option label %d 1", x)]] <- sprintf("option-%d-1", x)
      s_options[[sprintf("option label %d 2", x)]] <- sprintf("option-%d-2", x)

      # Change values for input$inSelect
      updateSelectInput(session, "inSelect", choices = s_options)

      # Can also set the label and select an item (or more than one if it's a
      # multi-select)
      updateSelectInput(session, "inSelect2",
        label = paste("Select label", x),
        choices = s_options,
        selected = sprintf("option-%d-2", x)
      )


      # Tabset input ==========================================================
      # Change the selected tab.
      # The tabsetPanel must have been created with an 'id' argument
      if (x_even) {
        updateTabsetPanel(session, "inTabset", selected = "panel2")
      } else {
        updateTabsetPanel(session, "inTabset", selected = "panel1")
      }
    })
  }
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
