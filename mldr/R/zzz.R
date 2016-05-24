.onAttach <- function(...) {
  if(interactive())
    packageStartupMessage("Enter mldrGUI() to launch mldr's web-based GUI")
}
