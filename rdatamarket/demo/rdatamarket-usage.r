oil <- dminfo("17tm")
print(oil)
series <- dmseries(oil, Country="Yemen")
lis <- dmlist(oil, Country=c("Algeria", "Angola"))

