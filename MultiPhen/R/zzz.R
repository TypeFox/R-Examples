.onAttach<- function (lib, pkg){
  if(interactive()){
    pkg.version = packageDescription("MultiPhen", fields = "Version")
    startup.txt = paste("   ==========================\n    MultiPhen", pkg.version, "is loaded\n   ==========================\n\nFor a description of its performance & to cite MultiPhen please refer to:\nO'Reilly et al. 2012. MultiPhen: Joint model of multiple phenotypes can increase discovery in GWAS.\nhttp://dx.plos.org/10.1371/journal.pone.0034861\n\n")
    packageStartupMessage(startup.txt)
  }
}
