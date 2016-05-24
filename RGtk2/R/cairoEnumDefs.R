CairoStatus<-c("success" = 0,
	"no-memory" = 1,
	"invalid-restore" = 2,
	"invalid-pop-group" = 3,
	"no-current-point" = 4,
	"invalid-matrix" = 5,
	"invalid-status" = 6,
	"null-pointer" = 7,
	"invalid-string" = 8,
	"invalid-path-data" = 9,
	"read-error" = 10,
	"write-error" = 11,
	"surface-finished" = 12,
	"surface-type-mismatch" = 13,
	"pattern-type-mismatch" = 14,
	"invalid-content" = 15,
	"invalid-format" = 16,
	"invalid-visual" = 17,
	"file-not-found" = 18,
	"invalid-dash" = 19,
	"invalid-dsc-comment" = 20,
	"invalid-index" = 21,
	"clip-not-representable" = 22)
storage.mode(CairoStatus) <- 'integer'
class(CairoStatus) <- 'enums' 

CairoOperator<-c("clear" = 0,
	"source" = 1,
	"over" = 2,
	"in" = 3,
	"out" = 4,
	"atop" = 5,
	"dest" = 6,
	"dest-over" = 7,
	"dest-in" = 8,
	"dest-out" = 9,
	"dest-atop" = 10,
	"xor" = 11,
	"add" = 12,
	"saturate" = 13)
storage.mode(CairoOperator) <- 'integer'
class(CairoOperator) <- 'enums' 

CairoFillRule<-c("winding" = 0,
	"even-odd" = 1)
storage.mode(CairoFillRule) <- 'integer'
class(CairoFillRule) <- 'enums' 

CairoLineCap<-c("butt" = 0,
	"round" = 1,
	"square" = 2)
storage.mode(CairoLineCap) <- 'integer'
class(CairoLineCap) <- 'enums' 

CairoLineJoin<-c("miter" = 0,
	"round" = 1,
	"bevel" = 2)
storage.mode(CairoLineJoin) <- 'integer'
class(CairoLineJoin) <- 'enums' 

CairoFontSlant<-c("normal" = 0,
	"italic" = 1,
	"oblique" = 2)
storage.mode(CairoFontSlant) <- 'integer'
class(CairoFontSlant) <- 'enums' 

CairoFontWeight<-c("normal" = 0,
	"bold" = 1)
storage.mode(CairoFontWeight) <- 'integer'
class(CairoFontWeight) <- 'enums' 

CairoAntialias<-c("default" = 0,
	"none" = 1,
	"gray" = 2,
	"subpixel" = 3)
storage.mode(CairoAntialias) <- 'integer'
class(CairoAntialias) <- 'enums' 

CairoSubpixelOrder<-c("default" = 0,
	"rgb" = 1,
	"bgr" = 2,
	"vrgb" = 3,
	"vbgr" = 4)
storage.mode(CairoSubpixelOrder) <- 'integer'
class(CairoSubpixelOrder) <- 'enums' 

CairoHintStyle<-c("default" = 0,
	"none" = 1,
	"slight" = 2,
	"medium" = 3,
	"full" = 4)
storage.mode(CairoHintStyle) <- 'integer'
class(CairoHintStyle) <- 'enums' 

CairoHintMetrics<-c("default" = 0,
	"off" = 1,
	"on" = 2)
storage.mode(CairoHintMetrics) <- 'integer'
class(CairoHintMetrics) <- 'enums' 

CairoPathDataType<-c("move-to" = 0,
	"line-to" = 1,
	"curve-to" = 2,
	"close-path" = 3)
storage.mode(CairoPathDataType) <- 'integer'
class(CairoPathDataType) <- 'enums' 

CairoFormat<-c("argb32" = 0,
	"rgb24" = 1,
	"a8" = 2,
	"a1" = 3,
	"rgb16-565" = 4)
storage.mode(CairoFormat) <- 'integer'
class(CairoFormat) <- 'enums' 

CairoExtend<-c("none" = 0,
	"repeat" = 1,
	"reflect" = 2)
storage.mode(CairoExtend) <- 'integer'
class(CairoExtend) <- 'enums' 

CairoContent<-c("color" = 4096,
	"alpha" = 8192,
	"color-alpha" = 12288)
storage.mode(CairoContent) <- 'integer'
class(CairoContent) <- 'enums' 

CairoFilter<-c("fast" = 0,
	"good" = 1,
	"best" = 2,
	"nearest" = 3,
	"bilinear" = 4,
	"gaussian" = 5)
storage.mode(CairoFilter) <- 'integer'
class(CairoFilter) <- 'enums' 

CairoSurfaceType<-c("image" = 0,
	"pdf" = 1,
	"ps" = 2,
	"xlib" = 3,
	"xcb" = 4,
	"glitz" = 5,
	"quartz" = 6,
	"win32" = 7,
	"beos" = 8,
	"directfb" = 9,
	"svg" = 10)
storage.mode(CairoSurfaceType) <- 'integer'
class(CairoSurfaceType) <- 'enums' 

CairoFontType<-c("toy" = 0,
	"ft" = 1,
	"win32" = 2,
	"atsui" = 3)
storage.mode(CairoFontType) <- 'integer'
class(CairoFontType) <- 'enums' 

CairoPatternType<-c("solid" = 0,
	"surface" = 1,
	"linear" = 2,
	"radial" = 3)
storage.mode(CairoPatternType) <- 'integer'
class(CairoPatternType) <- 'enums' 

CairoSvgVersion<-c("1-1" = 0,
	"1-2" = 1)
storage.mode(CairoSvgVersion) <- 'integer'
class(CairoSvgVersion) <- 'enums' 

CairoPsLevel<-c("2" = 0,
	"3" = 1)
storage.mode(CairoPsLevel) <- 'integer'
class(CairoPsLevel) <- 'enums' 

CairoTextClusterFlags<-c("backward" = 0)
storage.mode(CairoTextClusterFlags) <- 'integer'
class(CairoTextClusterFlags) <- 'enums' 

