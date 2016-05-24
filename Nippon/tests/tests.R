### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2013-01-23 20:13:15 umusus>
### Tests for kakasi() in Nippon package

suppressPackageStartupMessages(library(Nippon))
data(prefectures)
# fix later
try(kakasi(prefectures$name))




