## ----echo=FALSE----------------------------------------------------------
knitr::kable(matrix(c(
  'shiny:connected', 'socket', 'No', 'document',
  'shiny:disconnected', 'socket', 'No', 'document',
  'shiny:busy', '', 'No', 'document',
  'shiny:idle', '', 'No', 'document',
  'shiny:inputchanged', 'name, value, inputType', 'Yes', 'document',
  'shiny:message', 'message', 'Yes', 'document',
  'shiny:conditional', '', 'No', 'document',
  'shiny:bound', 'binding, bindingType', 'No', 'input/output element',
  'shiny:unbound', 'binding, bindingType', 'No', 'input/output element',
  'shiny:value', 'name, value, binding', 'Yes', 'output element',
  'shiny:error', 'name, error, binding', 'Yes', 'output element',
  'shiny:recalculating', '', 'No', 'output element',
  'shiny:recalculated', '', 'No', 'output element',
  'shiny:visualchange', 'visible, binding', 'No', 'output element',
  'shiny:updateinput', 'message, binding', 'Yes', 'input element'
), ncol = 4, byrow = TRUE,
dimnames = list(NULL, c('Name', 'Event Properties', 'Cancelable', 'Target'))))

