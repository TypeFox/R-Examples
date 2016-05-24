# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


loadtable = function (...){
  filename = file.choose()
  read.table (filename, ...)
}


