# Do a read-write-read cycle tests on all test files

files <- list.files("data", pattern="test", full.names=TRUE);

for (f in files) {
	fn <- as.filename(f);
	context(toupper(fn$ext));
	test_read_write_read(f);
}
