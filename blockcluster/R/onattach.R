.onAttach <-function(lib,pkg){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),
  "Version")
  packageStartupMessage("blockcluster version ", as.character(ver), " loaded\n\n----------------\nCopyright (C)  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>\nPlease post questions and bugs at: <https://gforge.inria.fr/forum/forum.php?forum_id=11190&group_id=3679>\n")
}
