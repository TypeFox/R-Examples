args(write.table)
write.table(ddd,"ddd.txt")
write.csv(ddd,"ddd.csv")
# this system call should work on a Mac or Linux machine
system("head -20 ddd.txt ddd.csv")
