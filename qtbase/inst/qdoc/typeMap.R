## low-level type conversion
typeMapBase <- c("QByteArray" = "raw", "const unsigned char \\*" = "raw",
                 "unsigned char" = "raw[1]", "signed char" = "character[1][1]",
                 "const char \\* const \\*" = "character",
                 "const char \\*" = "character[1]", char = "character[1][1]",
                 bool = "logical[1]", "unsigned int" = "numeric[1]",
                 "QList<int>" = "integer", "QList<unsigned int>" = "numeric",
                 "QVector<int>" = "integer",
                 "QVector<unsigned int>" = "numeric",
                 "QVector<double>" = "numeric", "QList<double>" = "numeric",
                 int = "integer[1]", "unsigned short" = "integer[1]",
                 short = "integer[1]", double = "numeric[1]",
                 "unsigned long long" = "numeric[1]",
                 "long long" = "numeric[1]",
                 "unsigned long" = "integer[1]", long = "integer[1]",
                 float = "numeric[1]", "void \\*" = "externalptr",
                 "QList<QString>" = "character",
                 "QMap<QString, " = "named list<",
                 "QHash<QString, " = "named list<",
                 "QPair" = "list[2]",
                 QString = "character[1]", QVariant = "ANY",
                 "QStringList" = "character",
                 QList = "list", "QGenericMatrix" = "matrix",
                 size_t = "numeric[1]", GLenum = "integer[1]",
                 GLfloat = "numeric[1]", GLint = "integer[1]",
                 GLuint = "numeric[1]", GLbitfield = "numeric[1]",
                 qreal = "numeric[1]")

## Goal: only replace when not part of a symbol
typeMap <- structure(paste("\\1", typeMapBase, "\\2", sep = ""),
                     names = names(typeMapBase))
## do these separate to get replacement numbers right
typeMap <- c(typeMap, "QMap<(.*?)>" = "\\1list<list[2]<\\2>>\\3",
             "QHash<(.*?)>" = "\\1list<list[2]<\\2>>\\3")
names(typeMap) <- paste("([^A-Za-z_])", names(typeMap), "([^A-Za-z_])",
                        sep = "")


