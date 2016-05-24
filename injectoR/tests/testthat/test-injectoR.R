# Changes of tests assertive of API integrity must necessitate bump in major version

describe ("Binder factory", {
  it ("Should be a function accepting a parent a callback function", {
    expect_true (is.function (binder));
    expect_equal (names (formals (binder)), c ('parent', 'callback'));
  });

  it ("Should create new binder", expect_true (is.environment (binder ())));

  it ("Should create new child binder", {
    b <- binder ();
    b2 <- binder (parent = b);
    expect_equal (parent.env (b2), b);
  });

  it ("Should inject callback", expect_true (binder (callback = function (b) {
    expect_true (is.environment (b));
    TRUE;
  })));
});

describe ("Singleton scope", {
  it ("Should be a function accepting provider", {
    expect_true (is.function (singleton));
    expect_equal (names (formals (singleton)), 'provider');
  });

  it ("Should provision", expect_equal (singleton (function () function () 'foo') () (), 'foo'));

  it ("Should cache provision", {
    called <- FALSE;
    scoped <- singleton (function (c = 0) {
                           expect_false (called); called <<- TRUE;
                           function () c <<- c + 1;
                         });
    injected <- scoped ();
    expect_equal (injected (), 1);
    expect_equal (injected (), 2);
    expect_equal (scoped () (), 3);
  });
});

describe ("Default scope", {
  it ("Should be a function accepting provider", {
    expect_true (is.function (default));
    expect_equal (names (formals (default)), 'provider');
  });

  it ("Should provide on each injection", {
    called <- 0;
    scoped <- default (function (c = 0) {
      called <<- called + 1;
      function () c <<- c + 1;
    });
    injected <- scoped ();
    expect_equal (injected (), 1);
    expect_equal (injected (), 2);
    expect_equal (scoped () (), 1);
    expect_equal (called, 2);
  });
});

describe ("Binding definition", {
  it ("Should be a function accepting ..., scope, and binder", {
    expect_true (is.function (define));
    expect_equal (names (formals (define)), c ('...', 'scope', 'binder'));
  });

  it ("Should define bindings", {
    b <- binder ();
    define (foo = function () 'foo', binder = b);
    expect_equal (ls (b), 'foo');
    expect_equal (b[[ 'foo' ]] (), 'foo');
  });

  it ("Should define multiple bindings", {
    b <- binder ();
    define (foo = function () 'foo', bar = function () 'bar', binder = b);
    expect_equal (c ('foo', 'bar') %in% ls (b), c (TRUE, TRUE))
    expect_equal (c (b$foo (), b$bar ()), c ('foo', 'bar'));
  });
});

describe ("Multibinder definition", {
  it ("Should be a function accepting key, scope, combine function, and binder", {
    expect_true (is.function (multibind));
    expect_equal (names (formals (multibind)), c ('key', 'scope', 'combine', 'binder'));
  });

  it ("Multibinder aggregator should be a function accepting ... and scope", {
    b <- binder ();
    a <- multibind ('foo', binder = binder ());
    expect_true (is.function (a));
    expect_equal (names (formals (a)), c ('...', 'scope'));
  });

  it ("Should define an unnamed multibinding", {
    b <- binder ();
    multibind ('foo', binder = b) (function () 'foo');
    expect_equal (ls (b), 'foo');
    expect_equal (b[[ 'foo' ]] (), list ('foo'));
  });

  it ("Should define a named multibinding", {
    b <- binder ();
    multibind ('bar', binder = b) (bar = function () 'bar');
    expect_equal (ls (b), 'bar');
    expect_equal (b[[ 'bar' ]] (), list (bar = 'bar'));
  });
});

describe ("Shim binding", {
  it ("Should be a function acepting ..., library.paths, callback", {
    expect_true (is.function (shim));
    expect_equal (names (formals (shim)), c ('...', 'library.paths', 'callback', 'binder'));
  });

  it ("Should shim package", {
    b <- binder ();
    shim ('injectoR', binder = b);
    expect_true ('define' %in% ls (b));
  });

  it ("Should shim named package", {
    b <- binder ();
    shim (i = 'injectoR', binder = b);
    expect_true ('i.define' %in% ls (b));
  });

  it ("Should shim and inject callback", expect_true (shim (i = 'injectoR', callback = function (i.inject) {
      expect_equal (names (formals (i.inject)), c ('callback', 'binder'));
      TRUE;
    }, binder = binder ())));

  it ("Shims should function",
      expect_equal (shim ('injectoR',
                          callback = function (define, inject, binder)
                            inject (function (f) f (5),
                              define (f = function (f) function (x) if (x < 1) 1 else x * f (x - 1),
                                      binder = binder ())),
                          binder = binder ()), 120));

  it ("Should throw on non string argument in library name list", tryCatch ({
    shim ('injectoR', function (inject) inject, binder = binder ());
    fail ('Did not throw error on non string argument to library list');
  }, error = function (e) expect_equal (e$message, 'Library name list must consist of strings only')));
});

describe ("Injection", {
  it ("Should be a function accepting callback and binder", {
    expect_true (is.function (inject));
    expect_equal (names (formals (inject)), c ('callback', 'binder'));
  });

  it ("Should inject indepentent callback", expect_equal (inject (function () 1), 1));

  it ("Should inject defined indepentent beans into callback", {
    b <- binder ();
    define (foo = function () 'foo', binder = b);
    define (bar = function () 'bar', binder = b);
    expect_equal (inject (function (foo, bar) list (foo, bar), b), list ('foo', 'bar'));
  });

  it ("Should inject defined beans with transitive dependencies into callback", {
    b <- binder ();
    define (three = function () 3, binder = b);
    define (power = function () p <- function (x, n) if (n < 1) 1 else x * p (x, n - 1), binder = b);
    define (cube = function (three, power) function (x) power (x, three), binder = b);
    expect_equal (inject (function (cube) cube (2), b), 8);
  });

  it ("Should inject multibound beans", {
    b <- binder ();

    multibind ('foo', binder = b) (one = function () 1)
    multibind ('foo', binder = b) (two = function () 2);
    multibind ('foo', binder = b) (three = function () 3);
    multibind ('foo', binder = b) (four = function () 4);
    expect_equal (inject (function (foo) foo, b), list (one = 1, two = 2, three = 3, four = 4));
    
    multibind ('bar', binder = b) (one = function () 1, two = function () 2);
    multibind ('bar', binder = b) (three = function () 3, four = function () 4);
    expect_equal (inject (function (bar) bar, b), list (one = 1, two = 2, three = 3, four = 4));
  });

  it ("Should injected multibound beans with respected scopes", {
    b <- binder ();
    c <- function () { c <- 0; function () c <<- c + 1; };
    multibind ('foo', binder = b) (p = c);
    multibind ('foo', binder = b) (s = c, scope = singleton);
    expect_equal (inject (function (foo) list (p = foo$p (), s = foo$s ()), b), list (p = 1, s = 1));
    expect_equal (inject (function (foo) list (p = foo$p (), s = foo$s ()), b), list (p = 1, s = 2));
  });

  it ("Should allow circular dependencies", {
    b <- binder ();
    define (f = function (f) function (x) if (x < 1) 1 else x * f (x - 1), binder = b);
    expect_equal (inject (function (f) f (6), b), 720);
  });

  it ("Should allow optional dependencies", {
    b <- binder ();
    define (g = function (n = 'stranger', t) paste ('Hello', n), binder = b);
    expect_equal (inject (function (g) g, b), 'Hello stranger');
    define (n = function () 'Bob', binder = b);
    expect_equal (inject (function (g) g, b), 'Hello Bob');
  });

  it ("Should throw error on access unbound dependency", {
    b <- binder ();
    define (t = function (h) h, binder = b);
    tryCatch ({
      inject (function (t) t, b);
      fail ('Did not throw on access to missing variable');
    }, error = function (e) {
      expect_equal (e$message, "Unbound dependency on h");
    });
  });

  it ("Should not throw error on no access to unbound dependency", {
    b <- binder ();
    define (t = function (h) 3, binder = b);
    expect_equal (inject (function (t) t, b), 3);
  });

  it ("missing() should return true for unbound dependencies", {
    expect_true (inject (function (q) missing (q)));
    expect_true (inject (function (w = 2) missing (w)));
  });

  it ("missing() should return false for bound dependencies", {
    b <- binder ();
    define (u = function () 4, binder = b);
    expect_false (inject (function (u) missing (u), b));
  });

  it ("missing() should revert back to original functionality in new environment", {
    b <- binder ();
    define (g = function () {
      t <- function (r) missing (r);
      expect_true (t ());
      expect_false (t (4));
      TRUE;
    }, binder = b);
    expect_true (inject (function (g) g, b));
  });

  it ("Should preserve injected function environment", {
    g <- 4;
    expect_equal (inject (function (d) d + g, define (d = function () 5, binder = binder ())), 9);
  });
});