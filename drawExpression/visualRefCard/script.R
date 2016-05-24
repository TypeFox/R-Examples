#library(visualR);
source("../R/drawExpression.R");

# anatomy of a vector

plotR <- function(expr, file) {
  par(mar=c(0,0,0,0));
  drawExpression(expr, filename=file);
}

# anatomy of a vector
v_char <- c("les", "gli", "los");
v_log <- c(T, F, T);
v_num <- c(10, .32, 2);
plotR("v_char", "v_char.pdf");
plotR("v_log", "v_log.pdf");
plotR("v_num", "v_num.pdf");

plotR("mode(v_char)", "v_char_type.pdf");
plotR("mode(v_log)", "v_log_type.pdf");
plotR("mode(v_num)", "v_num_type.pdf");

plotR("length(v_char)", "v_char_length.pdf");
plotR("length(v_log)", "v_log_length.pdf");
plotR("length(v_num)", "v_num_length.pdf");


v_char_names <- v_char;
names(v_char_names) <- c("fr", "it", "es");
plotR("v_char_names", "v_char_names.pdf");

# pb avec les listes récursives : list.rec, list.complex
plotR("c(1, 7, 3)", "c_1.pdf");
plotR("c(TRUE, F, T)", "c_2.pdf");
plotR("c(\"yes\", \"b\", \"no\")", "c_3.pdf");
plotR("rep(\"yes\", 4)", "rep.pdf");
plotR("seq(1, 4)", "seq.pdf");
plotR("rep(seq(1, 4), 2)", "rep_seq.pdf");

plotR("5:8", "operator_sequence.pdf");

plotR("c(1:3, FALSE, TRUE)", "conversion_1.pdf");
plotR("c(\"oui\", 1:3, FALSE)", "conversion.pdf");

v <- c(3, 102, 14);
plotR("nchar(v)", "conversion_nchar.pdf")
plotR("nchar(as.character(v))", "conversion_as_character.pdf")

# Extraction
v <- c(10, 7, 2);
i <- c(1,3);
j <- c(1,4);
plotR("v[2]", "extract_num.pdf")
plotR("v_char[i]", "extract_nums.pdf")

plotR("v_log[j]", "extract_nums_outsider.pdf")

k <- c(T,F,T);
plotR("v[k]", "extract_logical.pdf")
v_char2 <- c("a", "b", "c", "d");
l <- c(T, F);
plotR("v_char2[l]", "extract_logical_recycling.pdf")

names(v) <- c("a", "b", "c");
plotR("v[2]", "extract_num_names.pdf")
plotR("v[i]", "extract_nums_names.pdf")
plotR("v[\"b\"]", "extract_char_names.pdf")

plotR("v_char2[c(1, 1, 4, 2, 1)]", "extract_repeating.pdf")

# operators
plotR("5 + 6", "operator_plus.pdf");
plotR("5 - 6", "operator_minus.pdf");
plotR("5 * 6", "operator_time.pdf");
plotR("5 / 6", "operator_div.pdf");

plotR("5 < 6", "operator_strict_lt.pdf");
plotR("5 >= 6", "operator_gt_or_equal.pdf");
plotR("5 == 6", "operator_equality_num.pdf");
plotR("5 == 5", "operator_equality_num_TRUE.pdf");
plotR("T == F", "operator_equality_log.pdf");
plotR("\"oui\" == \"non\"", "operator_equality_char.pdf");

x1 <- c(4, 3, 8)
x2 <- c(32, 3, 2)
plotR("x1 + x2", "operator_plus_vectorized.pdf");
plotR("x1 == x2", "operator_equal_vectorized.pdf");

# Recycling

plotR("5:8 > 6", "operator_gt.pdf");
plotR("c(1, 2, 3, 4) + 1", "operator_add.pdf");
v_char_4 <- c("да", "yes", "ja", "si");
v_char_2 <- c("ja", "yes");
plotR("v_char_4 == v_char_2", "operator_recycling.pdf");

plotR("5:8 > 6", "operator_gt.pdf");
plotR("6 < 5:8", "operator_gt_inv.pdf");

# functions
x <- c(2, 1, 3, 4);
plotR("sum(x)", "sum.pdf");
x <- c(2, 1, 3, 4);
y <- c("les", "gli");
plotR("mean(x)", "mean.pdf");
plotR("range(x)", "range.pdf");
plotR("sort(x)", "sort.pdf");
plotR("sort(y)", "sort_char.pdf");
plotR("sort(x, decreasing=TRUE)", "sort2.pdf");
sort.v.names <- c("les", "gli", "los")
names(sort.v.names) <- c("fr", "it", "es");
plotR("sort(sort.v.names)", "sort_char_names.pdf");
plotR("rev(x)", "rev.pdf");
plotR("max(x)", "max.pdf");
plotR("min(x)", "min.pdf");

plotR("crossprod(x)", "crossprod.pdf");
plotR("cumsum(x)", "cumsum.pdf");

plotR("cumprod(x)", "cumprod.pdf");
plotR("var(x)", "var.pdf");
plotR("sd(x)", "sd.pdf");


f1 <- factor(c("det", "adj", "nc", "vb", "adv", "con", "pro", "adj", "vb"))
f2 <- factor(c("s", "s", "", "s", "-", "-", "p", "p", "p"))

plotR("table(f1)", "table_1.pdf");
plotR("table(f1)", "table_1b.pdf");
plotR("table(f1, f2)", "table_2.pdf");

plotR("table(c(1:4, 2:5, 4:7))", "table.pdf");

table.v <- c("a", "b", "c", "a", "c", "d");
#plotR("table(table.v)", "table.pdf");

# String: character vectors

v_char <- c("les", "i", "los");
paste.v.1 <- c("oui", "non");
paste.v.2 <- c("si", "no");
plotR("nchar(v_char)", "nchar.pdf");
plotR("paste(\"oui\", \"non\")", "paste.pdf");
plotR("paste(\"oui\", \"non\", sep=\"\")", "paste2.pdf");
plotR("paste(paste.v.1, paste.v.2)", "paste3.pdf");
plotR("paste(v_char, \"oui\")", "paste4.pdf");

paste.v.names <- c("les", "gli", "los")
names(paste.v.names) <- c("fr", "it", "es");
plotR("paste(\"(\", names(paste.v.names), \") \", paste.v.names, sep=\"\")", "paste5.pdf");



# index

order.v <- c(10, 2, 7, 15);
order.v.2 <- c("les", "gli", "los")
names(order.v.2) <- c("fr", "it", "es");
plotR("order(order.v)", "order.pdf");
plotR("order(order.v.2)", "order2.pdf");
plotR("which.min(x)", "which_min.pdf");
plotR("which.max(x)", "which_max.pdf");
which.v <- c(T, F, F, T, F, T);
plotR("which(which.v)", "which.pdf");

v1 <- c("le", "de", "un", "a");
v2 <- c(60, 100, 40, 30);
plotR("v1", "v_forms.pdf");
plotR("v2", "v_frequencies.pdf");
plotR("v1[which.max(v2)]", "which_max_2.pdf");

# list

plotR("list(1:3, TRUE)[1]", "list_extract_simple.pdf");
plotR("list(1:3, TRUE)[[1]]", "list_extract_double.pdf");
plotR("list(1:3)", "list.pdf");
plotR("list(1:3, TRUE)", "list_twotypes.pdf");
plotR("list(c(1, 2))", "list_onetype.pdf");
plotR("list(1, 2, 3, 4, list(1:4), 1:4)", "list_complex.pdf");
plotR("list(list(1:3))", "list_rec.pdf");
plotR("list(1, 2, 3, 4)", "list_enumerate.pdf");
l <- list(1:3, TRUE, "jour");
plotR("unlist(l)", "unlist.pdf");
plotR("as.list(1:4)", "aslist.pdf");

l.named <- l;
names(l.named) <- c("First", "Second", "Third");
plotR("names(l.named)", "list_names.pdf");

plotR("length(l)", "list_length.pdf");

plotR("mode(l)", "list_mode.pdf");

l.complexe <- list(1, list(1:3), "b")
plotR("l.complexe[[2]][[1]][3]", "list_successive_extraction.pdf");

li <- list(1:3, TRUE, "jour");
plotR("length(li[[1]])", "list_extraction_1.pdf");
plotR("length(li[1])", "list_extraction_2.pdf");

split.v <- 1:8
split.v2 <- c(1, 1, 1, 2, 2, 2, 3, 3);
plotR("split(split.v, split.v2)", "split.pdf");

# matrix

m <- matrix(1:6, 2);
plotR("nrow(m)", "nrow.pdf")
plotR("ncol(m)", "ncol.pdf")
plotR("dim(m)", "dim.pdf")

plotR("length(m)", "length_matrix.pdf")
plotR("mode(m)", "mode_matrix.pdf")
plotR("m[3]", "matrix_extraction_as_vector.pdf")

plotR("matrix(1:6, nrow=3)", "matrix.pdf")
plotR("matrix(1:6, 3)", "matrix_arg.pdf")
plotR("matrix(1:6, 3, byrow=T)", "matrix_byrow.pdf")
plotR("matrix(1:6, ncol=3)", "matrix_nbcol.pdf")

plotR("matrix(T, nrow=2, ncol=3)", "matrix_logical.pdf")

plotR("as.vector(m)", "as_vector.pdf")
plotR("m + 3", "matrix_operator.pdf")

m3 <- matrix(1:6, 2)
rownames(m3) <- LETTERS[1:2]
colnames(m3) <- c("un", "deux", "trois")

m2 <- matrix(11:16, 2)
plotR("cbind(m3, m3)", "cbind.pdf")
plotR("rbind(m3, m3)", "rbind.pdf")

m4 <- matrix(c(5,6,13:16), nrow=2, dimnames=list(c("A", "B"), c("trois", "quatre", "cinq")))
plotR("merge(m3, m4)", "merge.pdf")
plotR("t(m)", "t.pdf")

plotR("m[2,3]", "matrix_extraction.pdf")
plotR("m[2,2:3]", "matrix_extraction2.pdf")
plotR("m[1:2,2:3]", "matrix_extraction3.pdf")
plotR("m3[1,\"deux\"]", "matrix_extraction4.pdf")
b <- c(T,F);
plotR("m3[b,\"deux\"]", "matrix_extraction5.pdf")
plotR("m[2,2, drop=F]", "matrix_extraction_drop.pdf")

plotR("m[2,]", "matrix_extraction_rows.pdf")
plotR("m[,2]", "matrix_extraction_columns.pdf")

a <- 1:2
plotR("m3[1,a]", "matrix_extraction6.pdf")
plotR("m3[a,a]", "matrix_extraction7.pdf")
plotR("m3[2,2]", "matrix_extraction8.pdf")

plotR("rownames(m3)", "rownames.pdf")
plotR("colnames(m3)", "colnames.pdf")

plotR("rowSums(m)", "rowSums.pdf")
plotR("colSums(m)", "colSums.pdf")

m5 <- matrix(1:10, ncol=2)
rowsum.v <- c(2, 1, 2, 1, 3)
plotR("rowsum(m5, rowsum.v)", "rowsum.pdf")

plotR("sum(m)", "matrix_sum.pdf")


plotR("prop.table(m)", "prop_table1.pdf")
plotR("prop.table(m, margin=1)", "prop_table_row.pdf")
plotR("prop.table(m, margin=2)", "prop_table_column.pdf")

plotR("margin.table(m)", "margin_table1.pdf")
plotR("margin.table(m, margin=1)", "margin_table_row.pdf")
plotR("margin.table(m, margin=2)", "margin_table_column.pdf")

m <- matrix(c( 4, 1, 3, 3, 2, 1, 5, 9, 8), 3);
plotR("m[order(m[,1]),]", "matrix_order.pdf")


plotR("list(1:3, matrix(1:4, 2), 2)", "list_matrix.pdf");
plotR("list(1:3, matrix(1:4,2), c(T, F, T))", "ordre.pdf") # illustre les problèmes d'odres dans la résolution des opérateurs
plotR("1:10 + 2", "precedence.pdf")

# regexp

plotR("strsplit(c(\"un\", \"deux\", \"trois\"), \"[aeio]\")", "strsplit.pdf")
plotR("strsplit(c(\"un\", \"deux\", \"trois\"), c(\"u\", \"e\", \"r\"))", "strsplit_2.pdf")

plotR("substr(\"trois\", 2, 3)", "substr.pdf")
plotR("substr(c(\"trois\", \"quatre\"), 1, 3)", "substr_2.pdf")
plotR("substr(c(\"trois\", \"quatre\"), c(2,1), c(3, 4))", "substr_3.pdf")

grep.v <- c("un", "deux", "trois")
plotR("grep(\"[dt]\", grep.v)", "grep.pdf")
plotR("grepl(\"[dt]\", grep.v)", "grepl.pdf")

sub.v <- c("un", "deux", "trois")
plotR("sub(\"[ueaio]\", \"v\", sub.v)", "sub.pdf")
plotR("gsub(\"[ueaio]\", \"v\", sub.v)", "gsub.pdf")

#regexpr.v <- c("un", "deux", "trois")
#plotR("regexpr(\"[drs]\", regexpr.v)", "regexpr.pdf")
#plotR("gregexpr(\"[drs]\", regexpr.v)", "gregexpr.pdf")


#source("/Users/sylvainloiseau/workspace/seeR/R/seeR.R")
#mbx <- matrixBoxGrob(m, draw.index=T, draw.names=T)
#drawDetails(mbx)
#drawDetails(mbx)
#vbx <- vectorBoxGrob(c(oui="oui", non="non"), draw.index=T, draw.names=T)
#drawDetails(vbx)


m2 <- matrix(LETTERS[7:12], 2)
rownames(m2) <- c("row A", "row B")
colnames(m2) <- c("col a", "col b", "col c")
m <- matrix(LETTERS[1:6], 2)
rownames(m) <- c("row 1", "row 2")
colnames(m) <- c("un", "deux", "trois")

plotR("cbind(m, m2)", "cbind_2.pdf")


# factor

chars <- c("bleu", "rouge", "bleu", "vert", "rouge");
plotR("factor(chars)", "factor.pdf")
f <- factor(chars);
plotR("levels(f)", "levels.pdf")

forms1 <- c("pratique", "représentation", "linguistique", "sociale", "Guyane")
pos1 <- c("nc", "nc", "adj", "adj", "npr");
freqs1 <- c(10, 8, 13, 20, 5);
plotR("forms1", "forms1.pdf")
plotR("pos1", "pos1.pdf")
plotR("split(forms1, pos1)", "split2.pdf")
plotR("tapply(freqs1, pos1, mean)", "tapply.pdf")


# Data frame

df1 <- 1:3
df2 <- c("un", "deux", "deux")
df3 <- 1001:1003
df <- data.frame(col1=df1, col2=df2, col3=df3)
plotR("data.frame(col1=df1, col2=df2, col3=df3)", "dataframe.pdf")

plotR("df[1:2, 1:2]", "dataframe_extracting.pdf")
plotR("df[1, ]", "dataframe_extracting_row.pdf")
plotR("df[, 1]", "dataframe_extracting_column.pdf")

plotR("as.matrix(df)", "dataframe_as_matrix.pdf")
m <- matrix(1:6, 2)
plotR("as.data.frame(m)", "dataframe_as_data_frame.pdf")
plotR("as.list(df)", "dataframe_as_list.pdf")
plotR("as.data.frame(list(1:3, 11:13))", "dataframe_as_dataframe_from_list.pdf")


plotR("df$col1", "dataframe_list_like_extraction_1.pdf")
plotR("df[[\"col1\"]]", "dataframe_list_like_extraction_2.pdf")
plotR("df[[1]]", "dataframe_list_like_extraction_3.pdf")
plotR("df[1]", "dataframe_list_like_extraction_4.pdf")







plotR("c(7, 5) + 3", "graphical_conventions.pdf")
