#
#	Rdata.R
#Mon 27 Jun 2005 10:49:06 AM CEST
#system("cd ~/src/Rprivate ; ./exportR.sh");
#system("cd ~/src/Rprivate ; ./exportR.sh"); source("RgenericAll.R"); source("Rgenetics.R"); loadLibraries();

#
#	<§> abstract data functions
#

defined = function(x) exists(as.character(substitute(x)));
defined.by.name = function(name) { class(try(get(name), silent = T)) != 'try-error' }
# equivalent to i %in % v
is.in = function(i, v)(length((1:length(v))[v == i])>0)
rget = function(name, default = NULL, ..., pos = -1, envir = as.environment(pos)) {
	#obj = try(get(name, ...), silent = T);
	#r = if(class(obj) == 'try-error') default else obj;
	r = if (exists(name, where = pos, envir = envir)) get(name, ..., pos = pos, envir = envir) else default;
	r
}
firstDef = function(..., .fdInterpolate = F, .fdIgnoreErrors = F) {
	l = if (.fdInterpolate) c(...) else list(...);
	for (i in l) { if (!is.null(i) && (!.fdIgnoreErrors || class(i) != 'try-error')) return(i)};
	NULL
}
firstDefNA = function(..., .fdInterpolate = F){
	l = if (.fdInterpolate) c(...) else list(...);
	for (i in l) { if (!is.na(i)) return(i)};
	NULL
}
# <N> NULL behaviour
to.list = function(..., .remove.factors = T){
	r = if(is.null(...)) NULL else if (is.list(...)) c(...) else list(...);
	if (.remove.factors) {
		r = sapply(r, function(e)ifelse(is.factor(e), levels(e)[e], e));
	}
	r
}
# pretty much force to vector
#avu = function(v)as.vector(unlist(v))
avu = function(v, recursive = T) {
	r = if (is.list(v)) {
		nls = sapply(v, is.null);	# detects nulls
		# unlist removes NULL values -> NA
		unlist(sapply(1:length(v), function(i)if (nls[[i]]) NA else avu(v[[i]])));
	} else as.vector(v);
	if (!length(r)) return(NULL);
	r
}

assign.list = function(l, pos = -1, envir = as.environment(pos), inherits = FALSE, immediate = TRUE) {
	for (n in names(l)) {
		assign(n, l[[n]], pos, envir, inherits, immediate);
	}
}
eval.text = function(text, envir = parent.frame())eval(parse(text = c[1]), envir= envir);


# replace elements base on list
# l may be a list of lists with elements f (from) and t (to), when f is replaced with t
# if both, f and t arguments are not NULL, l will be ignored and f is replaced with t
vector.replace = function(v, l, regex = F, ..., f = NULL, t = NULL) {
	if (!is.null(f) & !is.null(t)) l = list(list(f = f, t = t));
	# replacments are given in f/t pairs
	if (all(sapply(l, length) == 2)) {
		from = list.key(l, "f");
		to = list.key(l, "t");
	} else {
		from = names(l);
		to = unlist(l);
	}
	for (i in 1:length(from)) {
		if (regex) {
			idcs = which(sapply(v, function(e)(length(fetchRegexpr(from[i], e, ...)) > 0)));
			v[idcs] = sapply(v[idcs], function(e)gsub(from[i], to[i], e));
		} else v[which(v == from[i])] = to[i];
	}
	v
}

vector.with.names = function(v, all_names, default = 0) {
	r = rep(default, length(all_names));
	names(r) = all_names;
	is = which.indeces(names(v), all_names, ret.na = T);
	r[is[!is.na(is)]] = v[!is.na(is)];
	r
}

# dir: direction of selection: 1: select rows, 2: select columns
mat.sel = function(m, v, dir = 1) {
	r = if (dir == 1)
		sapply(1:length(v), function(i)m[v[i], i]) else
		sapply(1:length(v), function(i)m[i, v[i]]);
	r
}

# rbind on list
sapplyId = function(l)sapply(l, function(e)e);

#
#	<§> string manipulation
#

say = function(...)cat(..., "\n");
printf = function(fmt, ...)cat(sprintf(fmt, ...));
join = function(v, sep = " ")paste(v, collapse = sep);
con = function(...)paste(..., sep="");
pastem = function(a, b, ..., revsort = T) {
	if (revsort)
		as.vector(apply(merge(data.frame(a = b), data.frame(b = a), sort = F), 1,
			function(e)paste(e[2], e[1], ...))) else
		as.vector(apply(merge(data.frame(a = a), data.frame(b = b), sort = F), 1,
			function(e)paste(e[1], e[2], ...)))
}
r.output.to.vector.int = function(s) {
	matches = gregexpr("(?<![\\[\\d])\\d+", s, perl=T);
	starts = as.vector(matches[[1]]);
	lengthes = attr(matches[[1]], "match.length");
	v = sapply(1:length(starts), function(i){ substr(s, starts[i], starts[i] + lengthes[i] -1) });
	as.integer(v)
}
r.output.to.vector.numeric = function(s) {
	matches = gregexpr("\\d*\\.\\d+", s, perl=T);
	starts = as.vector(matches[[1]]);
	lengthes = attr(matches[[1]], "match.length");
	v = sapply(1:length(starts), function(i){ substr(s, starts[i], starts[i] + lengthes[i] -1) });
	as.numeric(v)
}
readFile = function(path) { join(scan(path, what = "raw", sep = "\n", quiet = T), sep = "\n") };

Which.max = function(l) {
	if (all(!l)) return(NA);
	r = which.max(l);
	r
}
# capturesN: named captures; for each name in captureN put the captured value assuming names to be ordered
# captures: fetch only first capture per match <!> deprecated
# capturesAll: fetch all caputers for each match
fetchRegexpr = function(re, str, ..., ret.all = F, globally = T, captures = F, captureN = c(),
	capturesAll = F, maxCaptures = 9) {
	if (length(re) == 0) return(c());
	r = if (globally)
		gregexpr(re, str, perl = T, ...)[[1]] else
		regexpr(re, str, perl = T, ...);
	if (all(r < 0)) return(NULL);
	l = sapply(1:length(r), function(i)substr(str, r[i], r[i] + attr(r, "match.length")[i] - 1));
	if (captures) {
		l = sapply(l, function(e)gsub(re, '\\1', e, perl = T, fixed = F));
	} else if (length(captureN) > 0) {
		l = lapply(l, function(e) {
			r = sapply(1:length(captureN), function(i) {
				list(gsub(re, sprintf('\\%d', i), e, perl = T, fixed = F))
			});
			names(r) = captureN;
			r
		});
	} else if (capturesAll) {
		l = lapply(l, function(e) {
			cs = c();	# captures
			# <!> hack to remove zero-width assertions (no nested grouping!)
			#re = gsub('(\\(\\?<=.*?\\))|(\\(\\?=.*?\\))', '', re, perl = T, fixed = F);
			for (i in 1:maxCaptures) {
				n = gsub(re, sprintf('\\%d', i), e, perl = T, fixed = F);
				cs = c(cs, n);
			}
			cs
		});
		# trim list
		maxEls = maxCaptures - min(
			c(maxCaptures + 1, sapply(l, function(e)Which.max(rev(e != ''))))
		, na.rm = T) + 1;
		l = lapply(l, function(e)(if (maxEls > 0) e[1:maxEls] else NULL));
	}
	if (!ret.all) l = l[l != ""];
	l
}

splitString = function(re, str, ...) {
	r = gregexpr(re, str, perl = T, ...)[[1]];
	if (r[1] < 0) return(str);
	l = sapply(1:(length(r) + 1), function(i) {
		substr(str, ifelse(i == 1, 1, r[i - 1] + attr(r, "match.length")[i - 1]),
			ifelse(i > length(r), nchar(str), r[i] - 1))
	});
	l
}
quoteString = function(s)sprintf('"%s"', s)

mergeDictToString = function(d, s, valueMapper = function(s)
	ifelse(is.na(d[[n]]), '{\\bf Value missing}', d[[n]]),
	iterative = F, re = F, maxIterations = 100, doApplyValueMap = T, doOrderKeys = T, maxLength = 1e7) {
	ns = names(d);
	# proceed in order of decreasing key lengthes
	if (doOrderKeys) ns = ns[rev(order(sapply(ns, nchar)))];
	for (i in 1:maxIterations) {
		s0 = s;
		for (n in ns) {
			# counteract undocumented string interpolation
			subst = if (doApplyValueMap)
				gsub("[\\\\]", "\\\\\\\\", valueMapper(d[[n]]), perl = T)
				else d[[n]];
			# <!> quoting
			if (!re) n = sprintf("\\Q%s\\E", n);
			s = gsub(n, firstDef(subst, ""), s, perl = T, fixed = F);
			# <A> if any substitution was made, it is nescessary to reiterate ns to preserver order
			#	of substitutions
			if (iterative && s != s0) break;
		}
		if (!iterative || s == s0 || nchar(s) > maxLength) break;
	}
	s
}

mergeDictToVector = function(d, v) { unlist(ifelse(is.na(names(d[v])), v, d[v])) }

mergeDictToDict = function(dMap, dValues, ...) {
	r = lapply(dValues, function(v)mergeDictToString(dMap, v, ...));
	r
}

#r = getPatternFromStrings(DOC, '(?:\\nDOCUMENTATION_BEGIN:)([^\\n]+)\\n(.*?)(?:\\nDOCUMENTATION_END\\n)');
getPatternFromStrings = function(strings, pattern, keyIndex = 1) {
	r = lapply(strings, function(s) {
		ps = fetchRegexpr(pattern, s, capturesAll = T);
		listKeyValue(sapply(ps, function(e)e[[keyIndex]]), sapply(ps, function(e)e[-keyIndex]));
	});
	r
}

getPatternFromFiles = function(files, locations = NULL, ...) {
	strings = sapply(files, function(f)readFile(f, prefixes = locations));
	getPatternFromStrings(strings, ...);
}

#
#	hex strings
#

asc = function(x)strtoi(charToRaw(x), 16L);
character.as.characters = function(str) {
	sapply(str, function(s) sapply(1:nchar(s), function(i)substr(str, i, i)));
}

# bit_most_sig in bits
hex2int = function(str, bit_most_sig = 32) {
	cs = rev(sapply(character.as.characters(tolower(str)), asc));
	cms = bit_most_sig / 4;	# character containing most significant bit
	is = ifelse(cs >= asc('a'), cs - asc('a') + 10, cs - asc('0'));
	flipSign = (length(is) >= cms && is[cms] >= 8);
	if (flipSign) is[cms] = is[cms] - 8;
	r = sum(sapply(1:length(is), function(i)(is[i] * 16^(i-1))));
	if (flipSign) r = r - 2^(bit_most_sig - 1);
	r = if (r == - 2^(bit_most_sig - 1)) NA else as.integer(r);
	r
}

# chunk_size in bits
hex2ints = function(str, chunk_size = 32) {
	l = nchar(str);
	csc = chunk_size / 4;	# chunk_size in characters
	chunks = (l + csc - 1) %/% csc;
	r = sapply(1:chunks, function(i)hex2int(substr(str, (i - 1)*csc + 1, min(l, i*csc))));
	r
}

#
#	<§> binary numbers/n-adic numbers
#

ord2base = dec2base = function(o, digits = 5, base = 2) {
	sapply(1:digits, function(i){(o %/% base^(i-1)) %% base})
}
base2ord = base2dec = function(v, base = 2) {
	sum(sapply(1:length(v), function(i)v[i] * base^(i-1)))
}

ord2bin = dec.to.bin = function(number, digits = 5) ord2base(number, digits, base = 2);
bin2ord = bin.to.dec = function(bin) base2ord(bin, base = 2);

#
#	<Par> sequences
#

#
#	counts is a vector of lengthes of blocks and converted to pairs of indeces indicating the
#	first and last element in an appropriate vector
#
count2blocks = function(counts) {
	ccts = cumsum(counts);
	fidcs = c(1, ccts[-length(ccts)] + 1);
	blks = as.vector(rbind(fidcs, fidcs + counts - 1));
	blks
}

#
#	expand a block list - for example as from count2blocks - to a list of integers
#
expandBlocks = function(blks) {
	apply(matrix(blks, ncol = 2, byrow = T), 1, function(r) { r[1]:r[2] } )
}


splitListIndcs = function(M, N = 1, .compact = F, .truncate = T) {
	if (.truncate & M < N) N = M;
	if (.compact) {
		n = rep(ceiling(M / N), N);	# size of parts
		idcs = c(0, cumsum(n));
		idcs = idcs[idcs < M];
		idcs = c(idcs, M);
	} else {
		n = rep(floor(M / N), N);		# size of parts
		R = M - n[1] * N;
		n = n + c(rep(1, R), rep(0, N - R));
		idcs = c(0, cumsum(n));
	}
	idcs = cbind(idcs + 1, c(idcs[-1], 0))[-length(idcs), ];	# from:to in a row
	# <!> usual R degeneracy
	if (!is.matrix(idcs)) idcs = matrix(idcs, nrow = 1);
	idcs
}
splitListEls = function(l, N, returnElements = F) {
	idcs = splitListIndcs(length(l), N);
	li = apply(idcs, 1, function(r)(if (returnElements) l[r[1]:r[2]] else r[1]:r[2]));
	# <!> R ambiguity of apply return type
	if (is.matrix(li)) li = lapply(1:(dim(li)[2]), function(i)li[, i]);
	if (is.vector(li)) li = as.list(li);;
	li
}

# splitting based on fractions
# voting percentages to seats
#	simple algorithm based on size of residuals
splitSeatsForFractions = function(Nseats, fractions) {
	# number of parties
	Nparties = length(fractions);
	# fractional seats
	Nseats0 = fractions * Nseats;
	# garuantee one seat, otherwise round to nearest
	Nseats1 = ifelse (Nseats0 < 1, 1, round(Nseats0));
	# mismatch
	diff = sum(Nseats1) - Nseats;
	# redistribute deficit/overshoot
	if (diff != 0) {
		Nresid = sapply(Nseats0 - Nseats1, function(i)ifelse(i < 0, 1, i));
		subtr = order(Nresid, decreasing = diff < 0)[1:abs(diff)];
		# assume one round of correction is always sufficient <!>
		Nseats1[subtr] = Nseats1[subtr] - sign(diff);
	}
	Nseats1
}

# tranform number of elements (as from splitSeatsForFractions) into from:to per row in a matrix
counts2idcs = function(counts) {
	idcs = c(0, cumsum(counts));
	idcs = cbind(idcs + 1, c(idcs[-1], 0))[-length(idcs), ];
	idcs
}

# N is partitioned into fractions from p, where each element of p partitions the remaining part of N
# procedure makes sure to leave space for length(p) elements
cumpartition = function(N, p) {
	I = c();	# indeces within 1:N
	for (i in 1:length(p)) {
		# partition remaining space (ifelse), leave room for subsequent indeces
		Ii = floor(p[i] * (ifelse(i == 1, N, N - I[i - 1]) - (length(p) - i))) + 1;
		I = c(I, ifelse(i == 1, Ii, I[i - 1] + Ii));
	}
	as.integer(I)
}

#
#	<§> vector functions
#

# does the position exists in vector v
exists.pos = function(v, i)(is.vector(v) && !is.na(v[i]))

#
#	<par> lists
#

merge.lists = function(..., ignore.nulls = TRUE, listOfLists = F) {
	lists = if (listOfLists) c(...) else list(...);
	l1 = lists[[1]];
	if (length(lists) > 1) for (i in 2:length(lists)) {
		l2 = lists[[i]];
		for(n in names(l2)) {
			if (is.null(n)) print("Warning: tried to merge NULL key");
			if (!is.null(n) & (!ignore.nulls | !is.null(l2[[n]]))) l1[[n]] = l2[[n]];
		}
	}
	l1
}

# use.names preserves names and concatenates with lower level names
# reset sets names to top level names
unlist.n = function(l, n = 1, use.names = T, reset = F) {
	if (n > 0) for (i in 1:n) {
		ns = names(l);
		#names(l) = rep(NULL, length(l));	# <!> untested removal Tue Oct 19 17:11:53 2010
		l = unlist(l, recursive = F, use.names = use.names);
		if (reset) names(l) = ns;
	}
	l
}

# <N> obsolete, better: with(l, { ...})
instantiate.list = function(l, n = 1) {
	for (nm in names(l)) {
 		eval.parent(parse(file = "", text = sprintf("%s = %s", nm, deparse(l[[nm]]))), n = n);
# 		if (is.integer(l[[nm]])) {
# 			eval.parent(parse(file = "", text = sprintf("%s = %d", nm, l[[nm]])), n = n);
# 		} else if (is.numeric(l[[nm]])) {
# 			eval.parent(parse(file = "", text = sprintf("%s = %f", nm, l[[nm]])), n = n);
# 		} else {
# 			eval.parent(parse(file = "", text = sprintf("%s = \"%s\"", nm, l[[nm]])), n = n);
# 		};
	}
}
# assume a list of lists (aka vector of dicts) and extract a certain key from each of the lists
list.key = function(v, key, unlist = T, template = NULL, null2na = F) {
	l = lapply(v, function(i){
		if (is.list(i)) {
			if (is.null(i[[key]])) { if (null2na) NA else NULL } else i[[key]]
		} else template});
	if (unlist) l = unlist(l);
	l
}
# extract key path from list, general, recursive version
#	key path recursive worker
list.kprw = function(l, keys, unlist.pats, template, null2na, carryNames) {
	key = keys[1];
	# <p> extract key
	r = if (key != "*") {
		if (is.list(l)) {
			index = fetchRegexpr("\\A\\[\\[(\\d+)\\]\\]\\Z", key, captures = T);
			if (length(index) > 0) key = as.integer(index[[1]]);
			r = if (is.null(l[[key]])) { if (null2na) NA else NULL } else l[[key]];
			if (length(keys) > 1)
				list.kprw(r, keys[-1], unlist.pats[-1], template, null2na, carryNames) else r;
		} else return(template)
	} else {
		if (length(keys) > 1)
			lapply(l, function(sl)
				list.kprw(sl, keys[-1], unlist.pats[-1], template, null2na, carryNames)
			) else l;
	}
	# <p> unlisting
	if (!is.null(unlist.pats)) if (unlist.pats[1]) r = unlist.n(r, 1, reset = carryNames);
	r
}
# wrapper for list.kprw
# keyPath obeys EL1 $ EL2 $ ..., where ELn is '*' or a literal
# unlist.pat is pattern of truth values TR1 $ TR2 $..., where TRn is in 'T|F' and specifies unlist actions
# carryNames determines names to be carried over from the top level in case of unlist
list.kpr = function(l, keyPath, do.unlist = F, template = NULL,
	null2na = F, unlist.pat = NULL, carryNames = T, as.matrix = F) {
	keys = fetchRegexpr("[^$]+", keyPath);
	unlist.pats = if (!is.null(unlist.pat)) as.logical(fetchRegexpr("[^$]+", unlist.pat)) else NULL;

	r = list.kprw(l, keys, unlist.pats, template, null2na, carryNames);
	if (do.unlist) { r = unlist(r); }
	if (as.matrix) r = t(sapply(r, function(e)e));
	r
}
# extract key path from list
# <!> interface change: unlist -> do.unlist (Wed Sep 29 18:16:05 2010)
list.kp = function(l, keyPath, do.unlist = F, template = NULL, null2na = F) {
	r = list.kpr(l, sprintf("*$%s", keyPath), do.unlist = do.unlist, template, null2na);
	r
}

list.keys = function(l, keys, default = NA) {
	l = as.list(l);
	r = lapply(unlist(keys), function(key) if (is.null(l[[key]])) default else l[[key]]);
	r
}


# return list without listed keys
list.min  = function(l, keys) {
	l[-which.indeces(keys, names(l))]
}
# list generation on steroids (wraps other functions)
.list = function(l, .min = NULL) {
	if (!is.null(.min)) l = list.min(l, .min);
	l
}
# get apply
gapply = function(l, key, unlist = F)list.key(l, key, unlist)
# construct list as a dictionary for given keys and values
listKeyValue = function(keys, values) {
	if (length(keys) != length(values))
		stop("listKeyValues: number of provided keys does not match that of values");

	l = as.list(values);
	names(l) = keys;
	l
}
listInverse = function(l)listKeyValue(avu(l), names(l));

# name the list elements by the iterated vector elements ns (names)
nlapply = function(ns, f, ...) {
	if (is.list(ns)) ns = names(ns);
	r = lapply(ns, f, ...);
	names(r) = ns;
	r
}
# USE.NAMES logic reversed for sapply
sapplyn = function(l, f, ...)sapply(l, f, ..., USE.NAMES = F);
list.with.names = function(..., .key = 'name') {
	l = list(...);
	ns = names(l);
	r = nlapply(l, function(n) c(l[[n]], listKeyValue(.key, n)));
	r
}

#
#	<par> data type conversions
#

# assure m has at least 1 column
to.col = function(m) { if (is.null(dim(m))) t(t(m)) else m }
col.frame = function(l, col.name = 'value', minus = NULL, ignore.null = TRUE,
	do.paste = NULL, do.format = T, digits = 3, plus = NULL) {
	if (ignore.null) { for (n in names(l)) { if (is.null(l[[n]])) l[[n]] = NULL; } }
	if (!is.null(minus)) { for (n in minus) { l[[n]] = NULL; } }
	my.names = if (!is.null(plus)) plus else names(l);
	digits = if (length(digits) > 1) digits else rep(digits, length(l));
	if (!is.null(do.paste)) {
		if (do.format) {
			i = 1;
			for (n in my.names) { if (is.vector(l[[n]])) {
				l[[n]] = paste(sapply(l[[n]],
						function(e){if (is.numeric(e)) sprintf("%.*f", digits[i], e) else e}
					), collapse = do.paste)
				i = i + 1;
			}}
		} else {
			for (n in my.names) { if (is.vector(l[[n]])) l[[n]] = paste(l[[n]], collapse = do.paste) }
		}
	}
	f = as.data.frame(l);
	if (dim(f)[2] > length(col.name) && length(col.name) == 1)
		row.names(f) = paste(col.name, 1:dim(f)[1], sep = "")
	else row.names(f) = c(col.name);
	t(f)
}

# <i> collect recursively until list or data.frame
# convert list of lists to data frame (assuming identical keys for each sub list)
#	also works on list of vectors
listOfLists2data.frame = function(l, idColumn = "id", .names = NULL) {
	# collect keys
	keys = if (is.list(l[[1]]))
		sort(unique(as.vector(unlist(sapply(l, function(e)names(e)))))) else 1:length(l[[1]]);
	if (is.null(.names)) .names = keys;
	# row names
	rows = names(l);
	if (is.null(rows)) rows = 1:length(l);
	# build df

	#df = t(sapply(rows, function(r) { unlist(l[[r]][keys]) }));
	df = t(sapply(rows, function(r)list2df(l[[r]], keys)));
	df = if (!is.null(idColumn)) {
		data.frame.types(data.frame(..idColumn.. = rows, df),
			row.names = 1:length(rows), names = c(idColumn, .names));
	} else {
		data.frame.types(df, row.names = rows, names = .names);
	}
	df
}

# resetColNames: reset column names to names of first data frame
# colsFromFirstDf: take columns from the first data frame
# <i> improved algorithm: unlist everything, bind together: cave: data types,
#	strictly valid only for matrices
# Use cases:
#	list with named vectors: get data frame that contains all vectors with all possible names represented
#		listOfDataFrames2data.frame(cfs, colsFromUnion = T, do.transpose = T, idColumn = NULL);
listOfDataFrames2data.frame = function(l, idColumn = "id", do.unlist = T, direction = rbind,
	resetColNames = T, colsFromFirstDf = F, colsFromUnion = F, do.transpose = F) {
	# row names
	# <!> 2009-11-20 changed from: rows = firstDef(names(l), list(1:length(l)));
	rows = firstDef(names(l), 1:length(l));
	# columns
	ns = NULL;
	if (colsFromUnion) {
		ns = unique(unlist(sapply(l, names)));
		# get data.frame names
		ns = names(do.call(data.frame, listKeyValue(ns, rep(NA, length(ns)))));
		resetColNames = F;	# <!> mutually exclusive
	}
	# build df
	df = NULL;
	for (i in 1:length(rows)) {
		if (is.null(l[[i]])) next;	# ignore empty entries
		# <p> force to data frame
		df0 = if (do.transpose) as.data.frame(t(l[[i]])) else as.data.frame(l[[i]]);
		# <p> homogenize columns
		if (colsFromUnion) {
			# add missing columns
			ns0 = setdiff(ns, names(df0));
			df0 = do.call(data.frame, c(list(df0), listKeyValue(ns0, rep(NA, length(ns0)))));
			# correct order of columns
			df0 = df0[, ns];
		}
		if (!is.null(df)) {
			if (colsFromFirstDf) df0 = df0[, names(df)] else
			if (resetColNames) {
				names(df0) = if (is.null(idColumn)) names(df) else names(df)[-1];
			}
		}
		# <p> add id column
		df0 = if (is.null(idColumn)) df0 else cbind(rep(rows[i], dim(df0)[1]), df0);
		# <A> case differentiation should not me necessary
		df = if (i == 1) df0 else direction(df, df0);
	}
	if (!is.null(idColumn)) names(df)[1] = idColumn;
	if (do.unlist) for (n in names(df)) { df[[n]] = unlist(df[[n]]); }
	row.names(df) = NULL;
	df
}
cbindDataFrames = function(l, do.unlist = F) {
	listOfDataFrames2data.frame(l, idColumn = NULL, do.unlist, direction = cbind)
}
rbindDataFrames = function(l, do.unlist = F, useDisk = F, idColumn = NULL, transpose = F,
	resetColNames = F, colsFromFirstDf = F) {
	r = if (useDisk) {
		tempTable = tempfile();
		for (i in 1:length(l)) {
			d0 = l[[i]];
			if (class(d0) != 'data.frame') d0 = as.data.frame(d0);
			if (transpose) d0 = t(d0);
			if (!is.null(idColumn)) {
				d0 = data.frame(idColumn = names(l)[i], d0);
				names(d0)[1] = idColumn;
			}
			write.table(d0, file = tempTable, col.names = i == 1, append = i != 1, row.names = F);
		}
		read.table(tempTable, header = T, as.is = T);
	} else {
		listOfDataFrames2data.frame(l, idColumn = idColumn, do.unlist = do.unlist,
			direction = rbind, resetColNames = resetColNames, colsFromFirstDf = colsFromFirstDf)
	}
	r
}

# names2col assigns names of the list to a column of the data frame and values to the valueCol
list2df = function(l, cols = names(l), row.name = NULL, names2col = NULL, valueCol = 'value') {
	idcs = if (is.null(cols)) 1:length(l) else
		if (all(is.integer(cols))) cols else which.indeces(names(l), cols);
	if (is.null(cols) || all(is.integer(cols))) cols = paste('C', 1:length(l), sep = '');
	r = as.list(rep(NA, length(cols)));
	names(r) = cols;
	r[idcs] = l;
	r = as.data.frame(r, stringsAsFactors = F);
	if (!is.null(row.name)) row.names(r)[1] = row.name;
	if (!is.null(names2col)) {
		r = data.frame(name = names(r), value = unlist(r[1, ]), row.names = NULL, stringsAsFactors = F);
		names(r) = c(names2col, valueCol);
	}
	r
}

be.numeric = function(v)
	sapply(v, function(e)grepl('^-?\\d*(\\.\\d+)?(e-?\\d+)?$', e, ignore.case = T, perl = T));

list2df.print = function(l, valueCol = 'value', names2col = NULL, ..., digits = 3, scientific = 3) {
	l1 = list2df(l, valueCol = valueCol, names2col = names2col, ...);
	numericRows = be.numeric(l1[[valueCol]]);
	numbers = as.numeric(l1[[valueCol]][numericRows]);
	log10range = max(floor(log10(numbers))) - min(floor(log10(numbers)));
	#fmt = if (log10range > digits + 1) '%.*e' else '%.*f';
	numbers = sprintf(ifelse(abs(floor(log10(numbers))) > scientific, '%.*e', '%.*f'), digits, numbers);
	#numbers = sapply(numbers, function(n)sprintf(fmt, digits, n));
	separators = as.vector(names(l) == '' & is.na(l));
	l1[separators, names2col] = '-';
	l1[separators, valueCol] = '';
	l1[numericRows, valueCol] = numbers;
	print(l1);
}


rbind.list2df = function(d, l, row.name = NULL) {
	d = as.data.frame(d);
	r = list2df(l, names(d), row.name);
	r0 = rbind(d, r);
	r0
}

# d: data frame, l: list with names corresponding to cols, values to be searched for in columns
searchDataFrame = function(d, l, .remove.factors = T) {
	ns = names(l);
	d = d[, ns, drop = F];
	if (.remove.factors) {
		l = sapply(l, function(e)ifelse(is.factor(e), levels(e)[e], e));
		#d = apply(d, 2, function(col)(if (is.factor(col)) levels(col)[col] else col));
	}
	rs = which(as.vector(apply(apply(d, 1, function(r)(r == l)), 2, all)));
	rs
}

.df.cols = which.cols = function(d, cols, regex = F) {
	cols[is.numeric(cols)] = as.integer(cols[is.numeric(cols)]);
	cols[is.character(cols)] = which.indeces(cols[is.character(cols)], names(d), regex = regex);
	as.integer(cols)
}
# select columns by name
.df = function(d, names, regex = T, as.matrix = F) {
	cols = which.indeces(names, names(d), regex = regex);
	d0 = d[, cols, drop = F];
	# <t> simpler version:
	# d0 = d[, .df.cols(d, names, regex)];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}
.df.reorder = function(d, names, regex = T) {
	cols = .df.cols(d, names, regex);
	d0 = d[, c(cols, setdiff(1:dim(d)[2], cols))];
	d0
}
# remove columns by name
.dfm = function(d, names, regex = F, as.matrix = F) {
	cols = if (all(is.numeric(names))) as.integer(names) else which.indeces(names, names(d), regex = regex);
	d0 = d[, -cols, drop = F];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}
# remove rows by name
.dfrmr = function(d, names, regex = F, as.matrix = F) {
	rows = if (all(is.numeric(names)))
		as.integer(names) else
		which.indeces(names, row.names(d), regex = regex);
	d0 = d[-rows, , drop = F];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}

# remove rows/columns by name
.dfrm = function(d, rows = NULL, cols = NULL, regex = F, as.matrix = F) {
	d = as.data.frame(d);	# enforce data frame
	rows = if (is.null(rows)) 1:dim(d)[1] else
		-(if (all(is.numeric(rows))) as.integer(rows) else which.indeces(rows, row.names(d), regex = regex));
	cols = if (is.null(cols)) 1:dim(d)[2] else 
		-(if (all(is.numeric(cols))) as.integer(cols) else which.indeces(cols, names(d), regex = regex));
	d0 = d[rows, cols, drop = F];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}

# convert strings to data frame names
#	<i> create a data frame and extract names
.dfns = function(ns)gsub(':', '.', ns);

# manipulate list of vectors
# vectors i = 1,.., n with entries v_ij are represented as vector v_11, ..., v_n1, v_21, ...
meshVectors = function(...) {
	l = list(...);
	if (length(l) == 1) l = l[[1]];
	v = as.vector(t(sapply(l, function(v)unlist(v))));
	v
}

is.sorted = function(...)(!is.unsorted(...))
is.ascending = function(v) {
	if (length(v) < 2) return(T);
	for (i in 2:length(v)) if (v[i] <= v[i - 1]) return(F);
	return(T);
}

# pad a vector to length N
pad = function(v, N, value = NA)c(v, rep(value, N - length(v)));

#
#	<par> number sequences
#

rep.each = function(l, n) as.vector(sapply(l, function(e)rep(e, n)));
rep.each.row = function(m, n) matrix(rep.each(m, n), ncol = dim(m)[2])
rep.list = function(l, n) lapply(1:length(l), function(e)l);

# produce indeces for indeces positioned into blocks of blocksize of which count units exists
# example: expand.block(2, 10, 1:2) == c(1, 2, 11, 12)
expand.block = function(count, blocksize, indeces) {
	as.vector(apply(to.col(1:count), 1,
		function(i){ (i - 1) * blocksize + t(to.col(indeces)) }
	));
}

search.block = function(l, s) {
	b.sz = length(s);
	which(sapply(
		1:(length(l)/b.sz), function(i){all(l[((i - 1) * b.sz + 1):(i * b.sz)] == s)}
	));
}

#
#	<par> matrix functions
#

which.row = function(m, row) {
	cols = names(as.list(row));
	if (is.null(cols)) cols = 1:length(row);
	rows = 1:(dim(m)[1]);
	rows.found = rows[sapply(rows, function(i){ all(m[i, cols] == row) })];
	rows.found
}

# lsee:	list with searchees
# lsed:	list with searched objects
# inverse: lsed are regexes matched against lsee; pre-condition: length(lsee) == 1
# <!><t> cave: semantics changed as of 17.8.2009: return NA entries for unfound lsee-entries
# <!> match multi only implemented for merge = T
which.indeces = function(lsee, lsed, regex = F, ret.na = F, merge = T, match.multi = F, ...,
	inverse = F) {
	if (!length(lsed) || !length(lsee)) return(c());
	v = if (is.list(lsed)) names(lsed) else lsed;
	idcs = if (regex) {
		which(sapply(lsed, function(e)(
			if (inverse) length(fetchRegexpr(e, lsee, ...)) > 0 else
				any(sapply(lsee, function(see)(length(fetchRegexpr(see, e, ...)) > 0)))
		)))
	} else if (merge) {
		d0 = merge(data.frame(d = lsed, ix = 1:length(lsed)),
			data.frame(d = lsee, iy = 1:length(lsee)), all.y = T);
		idcs = if (match.multi) { d0$ix[unlist(sapply(lsee, function(e)which(d0$d == e)))]
		} else d0$ix[unlist(sapply(lsee, function(e)which(d0$d == e)[1]))];
#		} else d0$ix[order(d0$iy)]
		if (!ret.na) idcs = idcs[!is.na(idcs)];
		idcs
	} else {
		unlist(as.vector(sapply(lsee, function(e){
			w = which(e == v);
			if (!ret.na) return(w);
			ifelse(length(w), w, NA)
		})))
	};
	as.integer(idcs)
}

grep.vector = function(lsee, lsed, regex = F, ret.na = F, merge = T, match.multi = F, ..., inverse = F) {
	lsed[which.indeces(lsee, lsed, regex, ret.na, merge, match.multi, ..., inverse = inverse)]
}
grep.infixes = function(lsee, lsed, ...) {
	r = grep.vector(sapply(lsee, function(v)sprintf('^%s.*', v)), lsed, regex = T, inverse = F, ... );
	r
}

# force structure to be matrix (arrange vector into a row)
MR = function(m) {
	if (!is.matrix(m)) m = matrix(m, byrow = T, ncol = length(m));
	m
}
# force structure to be matrix (arrange vector into a columns)
MC = function(m) {
	if (!is.matrix(m)) m = matrix(m, byrow = F, nrow = length(m));
	m
}

#
#	<par> data processing
#

# like table but produce columns for all numbers 1..n (not only for counts > 0)
# cats are the expected categories
table.n = function(v, n, min = 1, categories = NULL) {
	if (is.null(categories)) categories = min:n;
	t = as.vector(table(c(categories, v)) - rep(1, length(categories)));
	t
}

#
#	<par> data types
#

to.numeric = function(x) { suppressWarnings(as.numeric(x)) }

# set types for columns: numeric: as.numeric
data.frame.types = function(df, numeric = c(), character = c(), factor = c(), integer = c(),
	do.unlist = T, names = NULL, row.names = NULL, reset.row.names = F, do.rbind = F, do.transpose = F,
	stringsAsFactors = F) {
	if (do.rbind) {
		#old code: df = t(sapply(df, function(e)e));
		lengthes = sapply(df, length);
		maxL = max(lengthes);
		df = t(sapply(1:length(df), function(i)c(df[[i]], rep(NA, maxL - lengthes[i]))));
	}
	if (do.transpose) df = t(df);
	df = as.data.frame(df, stringsAsFactors = stringsAsFactors);
	# set or replace column names
	if (!is.null(names)) {
		if (class(names) == "character") names(df)[1:length(names)] = names;
		if (class(names) == "list") names(df) = vector.replace(names(df), names);
	}
	if (do.unlist) for (n in names(df)) { df[[n]] = unlist(df[[n]]); }
	for (n in numeric) { df[[n]] = as.numeric(df[[n]]); }
	for (n in integer) { df[[n]] = as.integer(df[[n]]); }
	for (n in character) { df[[n]] = as.character(df[[n]]); }
	for (n in factor) { df[[n]] = as.factor(df[[n]]); }
	if (reset.row.names) row.names(df) = NULL;
	if (length(row.names) > 0) row.names(df) = row.names;
	df
}

Df_ = function(df0, headerMap = NULL, names = NULL, min_ = NULL) {
	r = df0;
	if (!is.null(names)) {
		if (class(names) == 'character') names(r)[1:length(names)] = names;
		if (class(names) == 'list') names(r) = vector.replace(names(r), names);
	}
	if (!is.null(headerMap)) names(r) = vector.replace(names(r), headerMap);
	if (!is.null(min_)) r = r[, -which.indeces(min_, names(r))];
	r
}

Df = function(..., headerMap = NULL, names = NULL, min_ = NULL) {
	r = data.frame(...);
	Df_(r, headerMap = headerMap, names = names, min_ = min_);
}

List_ = .List = function(l, min_ = NULL, rm.null = F) {
	if (!is.null(min_)) {
		i = which.indeces(min_, names(l));
		if (length(i) > 0) l = l[-i];
	}
	if (rm.null) {
		l = l[-which(sapply(l, is.null))];
	}
	l
}
List = function(..., min_ = NULL) {
	l = eval(list(...), envir = parent.frame(n = 1));
	.List(l, min_ = min_);
}

#
#	<par> sets and permutations
#

# this is the identity
inverseOrder = inversePermutation = function(p) {
# 	o = order(p);
# 	i = rep(NA, length(o));
# 	for (j in 1:length(o)) { i[o[j]] = j};
# 	i
	which.indeces(1:length(p), order(p))
}

# permutation is in terms of elements of l (not indeces)

applyPermutation = function(l, perm, from = 'from', to = 'to', returnIndeces = T) {
	# 1. bring perm[[from]] in the same order as l
	# 2. apply this order to perm[[to]]
	r0 = perm[[to]][order(perm[[from]])[inverseOrder(l)]];
	# 3. determine permutation going from l to r0
	r = order(l)[inverseOrder(r0)]
	if (!returnIndeces) r = l[r];
	r
}

order.df = function(df, cols = NULL, decreasing = F, na.last = F) {
	if (is.null(cols)) cols = 1:ncol(df);
	if (!is.numeric(cols)) cols = which.indeces(cols, names(df));
	orderText = sprintf("order(%s, decreasing = %s, na.last = %s)",
		paste(sapply(cols, function(i) { sprintf("df[, %d]", i) }), collapse = ", "
		), as.character(decreasing), as.character(na.last)
#		paste(sapply(cols, function(i) {
#			if (is.numeric(i)) sprintf("df[, %d]", i) else sprintf("df$%s", i) }), collapse = ", "
#		), as.character(decreasing), as.character(na.last)
	);
	o = eval(parse(text = orderText));
	#print(list(text = orderText, order = o, df=df));
	o
}

order.df.maps = function(d, maps, ..., regex = F) {
	cols = NULL;
	for (i in 1:length(maps)) {
		m = names(maps)[i];
		map = maps[[i]];
		keys = names(map);
		cols = c(cols, if (is.list(map)) {
			tempColName = sprintf("..order.df.maps.%04d", i);
			col = if (regex)
				sapply(d[[m]], function(e){ j = which.indeces(e, keys, regex = T, inverse = T)
					if (length(j) == 0) NA else map[[j]]
				}) else	as.character(map[d[[m]]]);
			col[col == "NULL"] = NA;
			d = data.frame(col, d, stringsAsFactors = F);
			names(d)[1] = tempColName;
		} else { m });
	}
	o = order.df(d, cols, ...);
	o
}

data.frame.union = function(l) {
	dfu = NULL;
	for (n in names(l)) {
		df = l[[n]];
		factor = rep(n, dim(df)[1]);
		dfu = rbind(dfu, cbind(df, factor));
	}
	dfu
}

Union = function(...) {
	l = list(...);
	r = NULL;
	for (e in l) { r = union(r, e); }
	r
}

# row bind of data.frames/matrices with equal number of cols
lrbind = function(l, as.data.frame = F, names = NULL) {
	d = dim(l[[1]])[2];
	v = unlist(sapply(l, function(m) unlist(t(m))));
	m = matrix(v, byrow = T, ncol = d);
	dimnames(m) = list(NULL, names(l[[1]]));
	if (as.data.frame) {
		m = data.frame(m);
		if (!is.null(names)) names(m) = names;
	}
	m
}

#
#	logic arrays
#

# same as in Rlab
count = function(v, na.rm = T) {
	if (na.rm) v = v[!is.na(v)];
	sum(v)	# old version: length((1:length(v))[v])
}
# v assumed to be logical
fraction = function(v, na.rm = T){
	if (na.rm) v = v[!is.na(v)];
	(sum(v)/length(v)) 	# old version: { length(v[v]) / length(v) }
}
# treat v as set
set.card = function(v)count(unique(v))

# null is false
#nif = function(b)(!(is.null(b) | is.na(b) | !b))
#nif = function(b)sapply(b, function(b)(!(is.null(b) || is.na(b) || !b)))
nif = function(b) {
	if (length(b) == 0) return(F);
	!(is.null(b) | is.na(b) | !b)
}
# null is true
#nit = function(b)(is.null(b) | is.na (b) | b)
#nit = function(b)sapply(b, function(b)(is.null(b) || is.na (b) || b))
nit = function(b) {
	if (length(b) == 0) return(T);
	is.null(b) | is.na (b) | b
}
# null is zero
#niz = function(e)ifelse(is.null(e) | is.na(e), 0, e)
niz = function(e)ifelse(is.null(e) | is.na(e), 0, e)

#
#	<p> complex structures
#

#
# Averaging a list of data frames per entry over list elements
#

# meanMatrices = function(d) {
# 	df = as.data.frame(d[[1]]);
# 	ns = names(df);
# 	# iterate columns
# 	dfMean = sapply(ns, function(n) {
# 		m = sapply(d, function(e)as.numeric(as.data.frame(e)[[n]]));
# 		mn = apply(as.matrix(m), 1, mean, na.rm = T);
# 		mn
# 	});
# 	dfMean
# }
meanMatrices = function(d) {
	dm = dim(d[[1]]);
	m0 = sapply(d, function(e)as.vector(e));
	m1 = apply(m0, 1, mean, na.rm = T);
	r = matrix(m1, ncol = dm[2]);
	r
}
meanVectors = function(d) {
	ns = names(d[[1]]);
	mn = apply(as.matrix(sapply(d, function(e)e)), 1, mean, na.rm = T);
	mn
}
meanList = function(l)mean(as.numeric(l));

meanStructure = function(l) {
	r = list();
	ns = names(l[[1]]);
	for (n in ns) {
		meanFct = if (is.matrix(l[[1]][[n]])) meanMatrices else
			if (length(l[[1]][[n]]) > 1) meanVectors else meanList;
		r[[n]] = meanFct(lapply(l, function(e)(e[[n]])));
	}
	r
}

#
#	<p> combinatorial functions
#

# form all combinations of input arguments as after being constraint to lists
# .first.constant designates whether the first list changes slowest (T) or fastest (F)
#	in the resulting data frame,
#	i.e. all other factors are iterated for a fixed value of l[[1]] (T) or not
# .constraint provides a function to filter the resulting data frame
merge.multi.list = function(l, .col.names = NULL, .col.names.prefix = "X",
	.return.lists = F, .first.constant = T, stringsAsFactors = F, .cols.asAre = F, .constraint = NULL) {
	# <p> determine column names of final data frame
	.col.names.generic = paste(.col.names.prefix, 1:length(l), sep = "");
	if (is.null(.col.names)) .col.names = names(l);
	if (is.null(.col.names)) .col.names = .col.names.generic;
	.col.names[.col.names == ""] = .col.names.generic[.col.names == ""];
	names(l) = .col.names;		# overwrite names
	# <p> construct combinations
	if (.first.constant) l = rev(l);
	df0 = data.frame();
	if (length(l) >= 1) for (i in 1:length(l)) {
		newNames = if (.cols.asAre) names(l[[i]]) else names(l)[i];
		# <p> prepare data.frame: handle lists as well as data.frames
		dfi = if (is.list(l[[i]])) unlist(l[[i]]) else l[[i]];
		df1 = data.frame.types(dfi, names = newNames, stringsAsFactors = stringsAsFactors);
		# <p> perform merge
		df0 = if (i > 1) merge(df0, df1) else df1;
	}
	if (.first.constant) df0 = df0[, rev(names(df0)), drop = F];
	if (.return.lists) df0 = apply(df0, 1, as.list);
	if (!is.null(.constraint)) {
		df0 = df0[apply(df0, 1, function(r).do.call(.constraint, as.list(r))), ];
	}
	df0
}

# analysis pattern using merge.multi.list
# i needs not to be an argument to f as .do.call strips excess arguments
iterateModels_old = function(modelList, f, ...,
	.constraint = NULL, .clRunLocal = T, .resultsOnly = F, .unlist = 0, lapply__ = clapply) {
	models = merge.multi.list(modelList, .constraint = .constraint);

	r = lapply__(1:dim(models)[1], function(i, ..., f__, models__) {
		args = c(list(i = i), as.list(models__[i, , drop = F]), list(...));
		.do.call(f__, args)
	}, ..., f__ = f, models__ = models);
	r = if (.resultsOnly) r else list(models = models, results = r);
	r = unlist.n(r, .unlist);
	r
}

# list of list, vector contains index for each of these lists to select elements from
#	these elements are merged and return
#	if sub-element is not a list, take name of sub-element and contruct list therefrom
#	namesOfLists controls whether, if a selected element is a list, its name is used instead
#		can be used to produce printable summaries
merge.lists.takenFrom = function(listOfLists, v) {
	l = list();
	ns = names(listOfLists);
	for (i in 1:length(v)) {
		new = if (!is.list(listOfLists[[i]]))
			listKeyValue(ns[i], listOfLists[[i]][v[i]]) else {
				t = listOfLists[[i]][[v[i]]];
				# list of vectors
				t = (if (!is.list(t)) listKeyValue(names(listOfLists[[i]])[v[i]], list(t)) else t);
				t
			}
		l = merge.lists(l, new);
	}
	l
}

# take indeces given by v from a nested list
# namesOfLists: take the name of the list at the position in v
#	if null, take first element or leave aggregation to the function aggregator
# aggregator: called with the final result, should flatten existing lists into characters
lists.splice = function(listOfLists, v, namesOfLists = F, aggregator = NULL) {
	ns = names(listOfLists);
	l = lapply(1:length(ns), function(i) {
		name = ns[i];
		e = listOfLists[[i]][v[i]];
		r = if (!is.list(e)) e else {
			f = if (namesOfLists) {
				g = names(e)[1];
				# handle name == NULL
				if (is.null(g)) {
					# make an attempt later to print element
					if (!is.null(aggregator)) e[[1]] else e[[1]][[1]]
				} else g
			} else e[[1]];
		}
		r
	});
	if (!is.null(aggregator)) l = aggregator(listKeyValue(ns, l), v, l);
	l
}

merge.multi.list.symbolic = function(modelList, ..., symbolizer = NULL) {
	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize, ...);
	r = data.frame.types(sapply(1:dim(models)[1], function(i, ...) {
		r = lists.splice(modelList, unlist(models[i, ]), namesOfLists = T, aggregator = symbolizer);
		r
	}), do.transpose = T, names = names(modelList));
	r
}

# <!> should be backwards compatible with iterateModels_old, not tested
iterateModels = function(modelList, f, ...,
	.constraint = NULL, .clRunLocal = T, .resultsOnly = F, .unlist = 0,
	lapply__ = Lapply, callWithList = F, symbolizer = NULL) {
	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize);
	models_symbolic = merge.multi.list.symbolic(modelList, symbolizer = symbolizer);
	if (!is.null(.constraint)) {
		sel = apply(models_symbolic, 1, function(r).do.call(.constraint, as.list(r)));
		models = models[sel, ];
		models_symbolic = models_symbolic[sel, ];
	}

	r = lapply__(1:dim(models)[1], function(i, ...) {
		modelPars = merge.lists.takenFrom(modelList, unlist(models[i, ]));
		if (callWithList) f(i, modelPars, ...) else {
			args = c(list(i = i), modelPars, list(...));
			.do.call(f, args)
		}
	}, ...);
	r = if (.resultsOnly) r else list(
		models = models,
		results = r,
		models_symbolic = models_symbolic
	);
	r = unlist.n(r, .unlist);
	r
}

iterateModelsExpand = function(modelList, .constraint = NULL) {
	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize, .constraint = .constraint);
	r = list(
		models = models,
		models_symbolic = merge.multi.list.symbolic(modelList, .constraint = .constraint)
	);
	r
}

# reverse effect of .retern.lists = T
#	list.to.df(merge.multi.list(..., .return.lists = T)) === merge.multi.list(..., .return.lists = F)
list.to.df = function(l)t(sapply(l, function(e)e))

merge.multi = function(..., .col.names = NULL, .col.names.prefix = "X",
	.return.lists = F, stringsAsFactors = F, .constraint = NULL) {
	merge.multi.list(list(...), .col.names = .col.names, .return.lists = .return.lists,
		stringsAsFactors = stringsAsFactors, .constraint = .constraint)
}

merge.multi.dfs = function(l, .first.constant = T, all = T, stringsAsFactors = F) {
	if (.first.constant) l = rev(l);
	if (length(l) >= 1) for (i in 1:length(l)) {
		df1 = data.frame.types(l[[i]], stringsAsFactors = stringsAsFactors);
		df0 = if (i > 1) merge(df0, df1, all = all) else df1;
	}
	if (.first.constant) df0 = df0[, rev(names(df0)), drop = F];
	df0
}

Merge = function(x, y, by = intersect(names(x), names(y)), ..., safemerge = T) {
	if (safemerge && length(by) == 0) {
		stop(sprintf('Merge: safemerge triggered. No common columns between "%s" and "%s"',
			join(names(x), sep = ','), join(names(y), sep = ',')))
	}
	r = merge(x = x, y = y, by = by, ...);
	r
}

# ids: variables identifying rows in final table
# vars: each combination of vars gets transformed to an own column
# <!> not tested for length(ids) > 1 || ength(rvars) > 1
# blockVars: should the repeated vars go in blocks or be meshed for vars
reshape.wide = function(d, ids, vars, blockVars = F, reverseNames = F, sort.by.ids = T) {
	# remaining vars
	rvars = setdiff(names(d), union(ids, vars));
	# levels of variables used in the long expansion
	levls = lapply(vars, function(v)unique(as.character(d[[v]])));
	# combinations at the varying vars as passed to vars
	cbs = merge.multi.list(levls, .col.names = vars, .first.constant = !blockVars);
	# repvars: repeated variables
	repvars = merge.multi.list(c(list(rvars), levls),
		.first.constant = !blockVars, .col.names = c("..var", vars));
	varnames = apply(repvars, 1, function(r)join(if (reverseNames) rev(r) else r, "."));

	r0 = data.frame.types(unique(d[, ids], drop = F), names = ids);
	r1 = data.frame.types(apply(r0, 1, function(r) {
		# <p> isolate rows which match to current id columns
		ids = which(apply(d[, ids, drop = F], 1, function(id)all(id == r)));
		d1 = d[ids, ];
		# <p> construct vector of repeated values
		vs = sapply(1:dim(cbs)[1], function(i) {
			# <A> should be equal to one
			row = which(apply(d1[, vars, drop = F], 1, function(r)all(r == cbs[i, ])));
			v = if (length(row) != 1) rep(NA, length(rvars)) else d1[row, rvars];
			v
		});
		# heed blockVars
		vs = as.vector(unlist(if (!blockVars) t(vs) else vs));
		vs
	}), do.transpose = T, names = varnames);
	r = data.frame(r0, r1);
	if (sort.by.ids) r = r[order.df(r, ids), ];
	row.names(r) = NULL;
	r
}

# factors: provide factor combinations explicitly for vars (otherwise split by '.', <i>)
#	valueColumn: name of new column
#	vars: vars to be transferred to the new valueColumn
#	factors: values in the new factorColumn column; by default the names of the variables transformed
# Example:
#	d0 = reshape.long(d, vars = 2:9, factors = c('case', 'ctr'), factorColumn = 'group',
#		valueColumn = c('AA', 'AG', 'GG', 'tot'));
#	reshape variables 2:9 (forming two groups: case/ctr), value of which is named 'group'
#		the shortened columns will get names valueColumn
reshape.long = function(d, vars = NULL, factorColumn = 'factor', valueColumn = 'value',
	factors = as.factor(vars), useDisk = F) {
	if (is.null(vars)) vars = names(d);
	# indeces of columns vars
	Ivars = .df.cols(d, vars);
	# remaining vars
	rvars = setdiff(1:length(names(d)), Ivars);
	# names thereof
	Nrvars = names(d)[rvars];

	# how wide are the blocks?
	S = length(vars) / length(factors);
	# columns of intermediate data.frame
	N = length(rvars);
	# create list of data frames
	dfs = lapply(1:nrow(d), function(i) {
		st = d[i, rvars];	# start of the new row
		df0 = data.frame(factors, value =  matrix(d[i, vars], nrow = length(factors), byrow = T));
		df1 = data.frame(st, df0, row.names = NULL);
		names(df1) = c(Nrvars, factorColumn, valueColumn);
		df1
	});
	r = rbindDataFrames(dfs, do.unlist = T, useDisk = useDisk);
	r
}

#
# <p> string functions
#

uc.first = firstUpper = function(s) {
	paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = "");
}

#
#	<p> factor transformations for data frames
#

dataExpandedNames = function(data) {
	dnames = unlist(lapply(names(data), function(v){
		if (is.factor(data[[v]])) paste(v, 1:(length(levels(data[[v]])) - 1), sep = "") else v;
	}));
	dnames
}
# model.matrix removes missing columns and could not be tweaked into working
dataExpandFactors = function(data, vars =  NULL) {
	if (is.null(vars)) vars = names(data);
	d0 = lapply(vars, function(v) {
		if (is.factor(data[[v]])) {
			ls = levels(data[[v]]);
			dcNa = rep(NA, length(ls) - 1);	# missing data coding
			dc = rep(0, length(ls) - 1);	# dummy coding
			sapply(data[[v]], function(e) {
				if (is.na(e)) return(dcNa);
				i = which(e == ls);
				if (i == 1) return(dc);
				dc[i - 1] = 1;
				return(dc);
			});
		} else data[[v]];
	});
	d0names = dataExpandedNames(data[, vars]);
	# re-transform data
	d1 = data.frame(matrix(unlist(lapply(d0, function(e)t(e))), ncol = length(d0names), byrow = F));
	names(d1) = d0names;
	d1
}
coefficientNamesForData = function(vars, data) {
	lnames = dataExpandedNames(data);	# names of levels of factors
	cnames = lnames[unlist(sapply(vars, function(v)which.indeces(v, lnames, regex = T)))];
	cnames
}

#
# <p> statistic oriented data frame manipulation
#

variableIndecesForData = function(d, vars, varsArePrefixes = T) {
	if (varsArePrefixes) vars = sapply(vars, function(e)sprintf('%s.*', e));
	which.indeces(vars, names(d), regex = T, match.multi = T)
}
variablesForData = function(d, vars, varsArePrefixes = T) {
	names(d)[variableIndecesForData(d, vars, varsArePrefixes)]
}

subData = function(d, vars, varsArePrefixes = T) {
	dfr = d[, variableIndecesForData(d, vars, varsArePrefixes), drop = F];
	dfr
}

subDataFromFormula = function(d, formula, responseIsPrefix = T, covariateIsPrefix = T) {
	resp = formula.response(formula);
	cov = formula.covariates(formula);
	ns = names(d);
	r = list(
		response = subData(d, resp, responseIsPrefix),
		covariate = subData(d, cov, covariateIsPrefix)
	);
	r
}

#
#	<p> graph functions
#

sub.graph.merge = function(df, leader, follower) {
	# next transitive step
	r0 = merge(df, data.frame(leader = leader, follower = follower), by = 'follower');
	# add new connections
	r1 = rbind(df, data.frame(follower = r0$leader.y, leader = r0$leader.x, cluster = r0$cluster));
	# symmetric closure
	r1 = rbind(r1, data.frame(follower = r1$leader, leader = r1$follower, cluster = r1$cluster))
	# form clusters by selecting min cluster number per connection
	r1 = r1[order(r1$cluster), ];
	row.names(r1) = 1:dim(r1)[1];
	r2 = unique(r1[, c('leader', 'follower')]);
	# select unique rows (first occurunce selects cluster)
	r = r1[as.integer(row.names(r2)), ];
	# pretty sort data frame
	r = r[order(r$cluster), ];
	r
}
# form clusters from a relationally defined hierarchy
sub.graph = function(df) {
	df = as.data.frame(df);
	names(df)[1:2] = c('follower', 'leader');
	df = df[order(df$follower), ];
	# seed clusters
	ids = sort(unique(df$follower));
	idsC = as.character(ids);
	counts = lapply(ids, function(id)sum(df$follower == id));
	names(counts) = idsC;
	clusters = unlist(sapply(idsC, function(id){ rep(as.integer(id), counts[[id]]) }));

	df = cbind(df, data.frame(cluster = rep(clusters, 2)));
	df = unique(rbind(df, data.frame(follower = df$leader, leader = df$follower, cluster = df$cluster)));
	# receiving frame
	df0 = df;
	# results with clusters
	i = 1;
	repeat {
		Nrows = dim(df0)[1];
		cls = df0$clusters;
		# add transitive connections
		df0 = sub.graph.merge(df0, follower = df0$leader, leader = df0$follower);
		if (dim(df0)[1] == Nrows && all(cls == df0$clusters)) break();
	}
	df0 = df0[order(df0$cluster), ];
	cIds = unique(df0$cluster);
	cls = lapply(cIds, function(id)unique(avu(df0[df0$cluster == id, c('follower', 'leader')])));
	cls
}

#
#	<p> formulas
#

# formula: formula as a character string with wildcard character '%'
# 	<!>: assume whitespace separation in formula between terms
#	<!>: write interaction with spaces <!> such as in:
#		f = 'MTOTLOS_binair ~ ZRES% + sq(ZRes%) + ( ZRES% )^2';
formula.re = function(formula, data, ignore.case = F) {
	vars = names(data);
	#regex = '(?:([A-Za-z_.]+[A-Za-z0-9_.]*)[(])?([A-Za-z.]+[%][A-Za-z0-9.%_]*)(?:[)])?';
	#			function names				(    regex						   )
	regex = '(?:([A-Za-z_.]+[A-Za-z0-9_.]*)[(])?([A-Za-z%.]+[A-Za-z0-9.%_]*)(?:[)])?';
	patterns = unique(fetchRegexpr(regex, formula, ignore.case = ignore.case));
	subst = nlapply(patterns, function(p) {
		comps = fetchRegexpr(regex, p, captureN = c('fct', 'var'), ignore.case = ignore.case)[[1]];
		p = sprintf("^%s$", gsub('%', '.*', comps$var));
		mvars = vars[sapply(vars, function(v)regexpr(p, v, perl = T, ignore.case = ignore.case)>=0)];
		if (comps$fct != '') {
			varf = sprintf('%s', paste(sapply(mvars, function(v)sprintf('%s(%s)', comps$fct, v)),
				collapse = " + "));
		} else {
			varf = sprintf('%s', paste(mvars, collapse = " + "));
		}
	});
	formulaExp = as.formula(mergeDictToString(subst, formula));
	formulaExp
}

formula.response = function(f) {
	#r = fetchRegexpr('[^\\s~][^~]*?(?=\\s*~)', if (is.formula(f)) deparse(f) else f);
	f = if (class(f) == 'formula') join(deparse(f), '') else f;
	r = as.character(fetchRegexpr('^\\s*([^~]*?)(?:\\s*~)', f, captures = T));
	# <p> version 2
	#fs = as.character(as.formula(as.character(f)));	# "~" "response" "covs"
	#r = fs[2];
	# <p> version 1
	#f = as.formula(f);
	#r = all.vars(f)[attr(terms(f), "response")];	# fails to work on 'response ~ .'
	r
}
formula.rhs = function(f)as.formula(
	fetchRegexpr('([~].*)', if (!is.character(f)) formula.to.character(f) else f, captures = T)
);
formula.covariates = function(f) {
	covs = all.vars(formula.rhs(f));
	#covs = setdiff(all.vars(as.formula(f)), formula.response(f));
	covs
}
formula.vars = function(f)union(formula.response(f), formula.covariates(f));

formula.nullModel = function(f) {
	r = formula.response(f);
	fn = as.formula(sprintf("%s ~ 1", r));
	fn
}
formula.to.character = function(f)join(deparse(f), '');

# <i> use terms.formula from a (a + ... + z)^2 formula
# <i> merge.multi.list(rep.list(covs, 2), .constraint = is.ascending)
covariatePairs = function(covs) {
	pairs = merge(data.frame(c1 = 1:length(covs)), data.frame(c2 = 1:length(covs)));
	pairs = pairs[pairs[, 1] > pairs[ ,2], ];
	df = data.frame(c1 = covs[pairs[, 1]], c2 = covs[pairs[, 2]]);
	df
}

formulaWith = function(repsonse = "y", covariates = "x")
	as.formula(sprintf("%s ~ %s", repsonse,  paste(covariates, collapse = "+")))
