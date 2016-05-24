### R code from vignette source 'poplite.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: poplite.Rnw:33-38
###################################################

library(poplite)
data(clinical)

ls()


###################################################
### code chunk number 2: poplite.Rnw:43-46
###################################################

head(clinical)



###################################################
### code chunk number 3: poplite.Rnw:51-54
###################################################

head(samples)



###################################################
### code chunk number 4: poplite.Rnw:59-62
###################################################

head(dna)



###################################################
### code chunk number 5: poplite.Rnw:69-78
###################################################

sample.tracking <- makeSchemaFromData(clinical, "clinical")

sample.tbsl <- makeSchemaFromData(samples, "samples")

sample.tracking <- append(sample.tracking,  sample.tbsl)

try(dna.tbsl <- makeSchemaFromData(dna, "dna"))



###################################################
### code chunk number 6: poplite.Rnw:84-90
###################################################

new.dna <- correct.df.names(dna)
dna.tbsl <- makeSchemaFromData(new.dna, "dna")

sample.tracking <- append(sample.tracking, dna.tbsl)



###################################################
### code chunk number 7: poplite.Rnw:100-105
###################################################

relationship(sample.tracking, from="clinical", to="samples") <- sample_id~sample_id
relationship(sample.tracking, from="clinical", to="dna") <-sample_id~sample_id

relationship(sample.tracking, from="samples", to="dna") <- sample_id+wave~sample_id+wave


###################################################
### code chunk number 8: poplite.Rnw:111-115
###################################################

sample.tracking.db <- Database(sample.tracking, tempfile())
populate(sample.tracking.db, dna=new.dna, samples=samples, clinical=clinical)



###################################################
### code chunk number 9: poplite.Rnw:127-134
###################################################

select(sample.tracking.db, .tables="dna")

select(sample.tracking.db, sample_id:lab_id, .tables="dna")

select(sample.tracking.db, .tables=c("clinical","dna"))



###################################################
### code chunk number 10: poplite.Rnw:139-146
###################################################

select(sample.tracking.db, sample_id:lab_id)

head(filter(sample.tracking.db, sex == "M" & var_wave_1 > 0))

filter(sample.tracking.db, sample_id == 97 & var_wave_1 > 0)



###################################################
### code chunk number 11: poplite.Rnw:151-154
###################################################

try(filter(sample.tracking.db, sample_id == 97))



###################################################
### code chunk number 12: poplite.Rnw:161-171
###################################################


select(sample.tracking.db, dna.sample_id)

select(sample.tracking.db, dna.sample_id:lab_id)

filter(sample.tracking.db, clinical.sample_id == 97)

filter(sample.tracking.db, clinical.status == 1 & dna.wave==2)



###################################################
### code chunk number 13: poplite.Rnw:177-200
###################################################

#poplite + dplyr
wave.1.samp.pop <- filter(select(sample.tracking.db, .tables=c("samples", "dna")), wave == 1)

#dplyr
src.db <- src_sqlite(dbFile(sample.tracking.db), create = F)
samp.tab <- tbl(src.db, "samples")
dna.tab <- tbl(src.db, "dna")
wave.1.samp.dplyr <- inner_join(filter(samp.tab, wave == 1), dna.tab,
    by=c("sample_id", "wave"))

library(RSQLite)

#RSQLite
samp.db <- dbConnect(SQLite(), dbFile(sample.tracking.db))
wave.1.samp.sql <- dbGetQuery(samp.db, 'SELECT * FROM samples JOIN dna
    USING (sample_id, wave) WHERE wave == 1')
dbDisconnect(samp.db)

all.equal(as.data.frame(wave.1.samp.pop), wave.1.samp.sql)

all.equal(as.data.frame(wave.1.samp.dplyr), wave.1.samp.sql)



###################################################
### code chunk number 14: poplite.Rnw:209-226
###################################################

gender <- data.frame(sex=unique(clinical$sex), stringsAsFactors=F)
gend.tbsl <- makeSchemaFromData(gender, "gender")

sample.tracking <- append(gend.tbsl, sample.tracking)

relationship(sample.tracking, from="gender", to="clinical") <- .~sex

sample.tracking.db <- Database(sample.tracking, tempfile())
populate(sample.tracking.db, dna=new.dna, samples=samples, clinical=clinical, gender=gender)

select(sample.tracking.db, .tables="gender")

head(select(sample.tracking.db, .tables="clinical"))

select(sample.tracking.db, .tables=c("clinical", "gender"))



###################################################
### code chunk number 15: poplite.Rnw:239-243
###################################################

library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")


###################################################
### code chunk number 16: poplite.Rnw:250-264
###################################################

populate.ref.table <- function(vcf.obj)
{
    ref.dta <- cbind(
                    seqnames=as.character(seqnames(vcf.obj)),
                    as.data.frame(ranges(vcf.obj))[,c("start", "end")],
                    ref=as.character(ref(vcf.obj)),
                    stringsAsFactors=FALSE
                    )
    return(ref.dta)
}

vcf.sc <- makeSchemaFromFunction(populate.ref.table, "reference", vcf.obj=vcf[1:5])



###################################################
### code chunk number 17: poplite.Rnw:269-287
###################################################

populate.allele.table <- function(vcf.obj)
{
    exp.obj <- expand(vcf.obj)
    ref.dta <- cbind(
                    seqnames=as.character(seqnames(exp.obj)),
                    as.data.frame(ranges(exp.obj))[,c("start", "end")],
                    ref=as.character(ref(exp.obj)),
                    alt=as.character(alt(exp.obj)),
                    stringsAsFactors=FALSE
                    )
    return(ref.dta)
}

allele.sc <- makeSchemaFromFunction(populate.allele.table, "alleles", vcf.obj=vcf[1:5])

vcf.sc <- poplite::append(vcf.sc, allele.sc)



###################################################
### code chunk number 18: poplite.Rnw:294-317
###################################################

populate.samp.alt.table <- function(vcf.obj)
{
    temp.vrange <- as(vcf.obj, "VRanges")
    
    ret.dta <- cbind(
                        seqnames=as.character(seqnames(temp.vrange)),
                        as.data.frame(ranges(temp.vrange))[,c("start", "end")],
                        ref=ref(temp.vrange),
                        alt=alt(temp.vrange),
                        sample=as.character(sampleNames(temp.vrange)),
                        allele_count=sapply(strsplit(temp.vrange$GT, "\\|"),
                                function(x) sum(as.integer(x), na.rm=T)),
                        stringsAsFactors=F
                        )
    
    return(ret.dta[ret.dta$allele_count > 0,])
}

geno.all.sc <- makeSchemaFromFunction(populate.samp.alt.table, "sample_alleles", vcf.obj=vcf[1:5])

vcf.sc <- poplite::append(vcf.sc, geno.all.sc)



###################################################
### code chunk number 19: poplite.Rnw:325-332
###################################################

relationship(vcf.sc, from="reference", to="alleles") <- .~seqnames+start+end+ref

relationship(vcf.sc, from="reference", to="sample_alleles") <- .~seqnames+start+end+ref

relationship(vcf.sc, from="alleles", to="sample_alleles") <- .~.reference+alt



###################################################
### code chunk number 20: poplite.Rnw:340-371
###################################################

vcf.db <- Database(vcf.sc, tempfile())

populate(vcf.db, vcf.obj=vcf[1:1000])

populate(vcf.db, vcf.obj=vcf[1001:2000])

pop.res <- as.data.frame(poplite::select(vcf.db, .tables=tables(vcf.db)))

vrange.tab <- as(vcf[1:2000], "VRanges")

vrange.dta <- data.frame(seqnames=as.character(seqnames(vrange.tab)),
                         start=start(vrange.tab),
                         end=end(vrange.tab),
                         ref=as.character(ref(vrange.tab)),
                         alt=as.character(alt(vrange.tab)),
                         sample=as.character(sampleNames(vrange.tab)),
                         allele_count=sapply(strsplit(vrange.tab$GT, "\\|"),
                                function(x) sum(as.integer(x), na.rm=T)),
                         stringsAsFactors=F)

vrange.dta <- vrange.dta[vrange.dta$allele_count > 0,]

vrange.dta <- vrange.dta[do.call("order", vrange.dta),]

sub.pop.res <- pop.res[,names(vrange.dta)]

sub.pop.res <- sub.pop.res[do.call("order", sub.pop.res),]

all.equal(sub.pop.res, vrange.dta, check.attributes=F)



###################################################
### code chunk number 21: poplite.Rnw:374-375
###################################################
sessionInfo()


