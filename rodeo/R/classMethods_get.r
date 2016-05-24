
# Methods to return names
rodeo$methods( namesVars= function() {
  "Returns the names of the model's state variables (vector of mode character)."
  as.character(.self$.vars$name) })
rodeo$methods( namesPars= function() {
  "Returns the names of the model's parameters (vector of mode character)."
  as.character(.self$.pars$name) })
rodeo$methods( namesFuns= function() {
  "Returns the names of functions appearing in the ODE's righthand sides
  (vector of mode character)."
  as.character(.self$.funs$name) })
rodeo$methods( namesPros= function() {
  "Returns the names of simulated processes (vector of mode character)."
  as.character(.self$.pros$name) })

# Methods to return lengths
rodeo$methods( lenVars= function() {
  "Returns the number of state variables appearing in the ODE system (integer)."
  nrow(.self$.vars) })
rodeo$methods( lenPars= function() {
  "Returns the number of parameters appearing in the ODE system (integer)."
  nrow(.self$.pars) })
rodeo$methods( lenFuns= function() {
  "Returns the number of functions appearing in the ODE's righthand sides
  (integer)."
  nrow(.self$.funs) })
rodeo$methods( lenPros= function() {
  "Returns the number of simulated processes (integer)."
  nrow(.self$.pros) })

# Methods to return entire tables
rodeo$methods( getVars= function() {
  "Returns the declaration of the model's state variables (as a data frame)."
  .self$.vars })
rodeo$methods( getPars= function() {
  "Returns the declaration of the model's parameters (as a data frame)."
  .self$.pars })
rodeo$methods( getFuns= function() {
  "Returns the declaration of functions appearing in the ODE's righthand sides
  (as a data frame)."
  .self$.funs })
rodeo$methods( getPros= function() {
  "Returns the declaration of simulated processes (as a data frame)."
  .self$.pros })
rodeo$methods( getStoi= function() {
  "Returns the model's stoichiometric factors (as a data frame)."
  .self$.stoi })

