# .First.lib <-function (lib, pkg) { 
.onLoad <-function (lib, pkg) { 
library.dynam("timetools", pkg, lib) 
} 
