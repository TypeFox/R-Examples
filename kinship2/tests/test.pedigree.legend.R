## Test pedigree.legend function


require(kinship2)


data(minnbreast)
pedMN <- with(minnbreast, pedigree(id, fatherid, motherid, sex,famid=famid,
                         affected=cbind(cancer, bcpc, proband)))

mnf8 <- pedMN['8']
pout <- plot(mnf8)


pedigree.legend(mnf8, radius=.8, location="bottomright")
pedigree.legend(mnf8, radius=.5, location="topright")

pedigree.legend(mnf8, radius=.6, location="topleft")

pedigree.legend(mnf8, new=FALSE, cex=1.5)
