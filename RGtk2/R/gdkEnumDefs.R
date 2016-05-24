GdkCursorType<-c("x_cursor" = 0,
	"arrow" = 2,
	"based_arrow_down" = 4,
	"based_arrow_up" = 6,
	"boat" = 8,
	"bogosity" = 10,
	"bottom_left_corner" = 12,
	"bottom_right_corner" = 14,
	"bottom_side" = 16,
	"bottom_tee" = 18,
	"box_spiral" = 20,
	"center_ptr" = 22,
	"circle" = 24,
	"clock	" = 26,
	"coffee_mug" = 28,
	"cross" = 30,
	"cross_reverse" = 32,
	"crosshair" = 34,
	"diamond_cross" = 36,
	"dot" = 38,
	"dotbox" = 40,
	"double_arrow" = 42,
	"draft_large" = 44,
	"draft_small" = 46,
	"draped_box" = 48,
	"exchange" = 50,
	"fleur" = 52,
	"gobbler" = 54,
	"gumby" = 56,
	"hand1" = 58,
	"hand2" = 60,
	"heart" = 62,
	"icon" = 64,
	"iron_cross" = 66,
	"left_ptr" = 68,
	"left_side" = 70,
	"left_tee" = 72,
	"leftbutton" = 74,
	"ll_angle" = 76,
	"lr_angle" = 78,
	"man" = 80,
	"middlebutton" = 82,
	"mouse" = 84,
	"pencil" = 86,
	"pirate" = 88,
	"plus" = 90,
	"question_arrow" = 92,
	"right_ptr" = 94,
	"right_side" = 96,
	"right_tee" = 98,
	"rightbutton" = 100,
	"rtl_logo" = 102,
	"sailboat" = 104,
	"sb_down_arrow" = 106,
	"sb_h_double_arrow" = 108,
	"sb_left_arrow" = 110,
	"sb_right_arrow" = 112,
	"sb_up_arrow" = 114,
	"sb_v_double_arrow" = 116,
	"shuttle" = 118,
	"sizing" = 120,
	"spider		" = 122,
	"spraycan" = 124,
	"star" = 126,
	"target" = 128,
	"tcross" = 130,
	"top_left_arrow" = 132,
	"top_left_corner" = 134,
	"top_right_corner" = 136,
	"top_side" = 138,
	"top_tee" = 140,
	"trek" = 142,
	"ul_angle" = 144,
	"umbrella" = 146,
	"ur_angle" = 148,
	"watch" = 150,
	"xterm" = 152,
	"last-cursor" = 153,
	"gdk-cursor-is-pixmap" = -1)
storage.mode(GdkCursorType) <- 'integer'
class(GdkCursorType) <- 'enums' 

GdkDragProtocol<-c("motif" = 0,
	"xdnd" = 1,
	"rootwin" = 2,
	"none" = 3,
	"win32-dropfiles" = 4,
	"ole2" = 5,
	"local" = 6)
storage.mode(GdkDragProtocol) <- 'integer'
class(GdkDragProtocol) <- 'enums' 

GdkFilterReturn<-c("continue" = 0,
	"translate" = 1,
	"remove" = 2)
storage.mode(GdkFilterReturn) <- 'integer'
class(GdkFilterReturn) <- 'enums' 

GdkEventType<-c("nothing" = -1,
	"delete" = 0,
	"destroy" = 1,
	"expose" = 2,
	"motion-notify" = 3,
	"button-press" = 4,
	"2button-press" = 5,
	"3button-press" = 6,
	"button-release" = 7,
	"key-press" = 8,
	"key-release" = 9,
	"enter-notify" = 10,
	"leave-notify" = 11,
	"focus-change" = 12,
	"configure" = 13,
	"map" = 14,
	"unmap" = 15,
	"property-notify" = 16,
	"selection-clear" = 17,
	"selection-request" = 18,
	"selection-notify" = 19,
	"proximity-in" = 20,
	"proximity-out" = 21,
	"drag-enter" = 22,
	"drag-leave" = 23,
	"drag-motion" = 24,
	"drag-status" = 25,
	"drop-start" = 26,
	"drop-finished" = 27,
	"client-event" = 28,
	"visibility-notify" = 29,
	"no-expose" = 30,
	"scroll" = 31,
	"window-state" = 32,
	"setting" = 33,
	"owner-change" = 34,
	"grab-broken" = 35,
	"gdk-damage" = 36)
storage.mode(GdkEventType) <- 'integer'
class(GdkEventType) <- 'enums' 

GdkVisibilityState<-c("unobscured" = 0,
	"partial" = 1,
	"fully-obscured" = 2)
storage.mode(GdkVisibilityState) <- 'integer'
class(GdkVisibilityState) <- 'enums' 

GdkScrollDirection<-c("up" = 0,
	"down" = 1,
	"left" = 2,
	"right" = 3)
storage.mode(GdkScrollDirection) <- 'integer'
class(GdkScrollDirection) <- 'enums' 

GdkNotifyType<-c("ancestor" = 0,
	"virtual" = 1,
	"inferior" = 2,
	"nonlinear" = 3,
	"nonlinear-virtual" = 4,
	"unknown" = 5)
storage.mode(GdkNotifyType) <- 'integer'
class(GdkNotifyType) <- 'enums' 

GdkCrossingMode<-c("normal" = 0,
	"grab" = 1,
	"ungrab" = 2)
storage.mode(GdkCrossingMode) <- 'integer'
class(GdkCrossingMode) <- 'enums' 

GdkPropertyState<-c("new-value" = 0,
	"delete" = 1)
storage.mode(GdkPropertyState) <- 'integer'
class(GdkPropertyState) <- 'enums' 

GdkSettingAction<-c("new" = 0,
	"changed" = 1,
	"deleted" = 2)
storage.mode(GdkSettingAction) <- 'integer'
class(GdkSettingAction) <- 'enums' 

GdkFontType<-c("font" = 0,
	"fontset" = 1)
storage.mode(GdkFontType) <- 'integer'
class(GdkFontType) <- 'enums' 

GdkCapStyle<-c("not-last" = 0,
	"butt" = 1,
	"round" = 2,
	"projecting" = 3)
storage.mode(GdkCapStyle) <- 'integer'
class(GdkCapStyle) <- 'enums' 

GdkFill<-c("solid" = 0,
	"tiled" = 1,
	"stippled" = 2,
	"opaque-stippled" = 3)
storage.mode(GdkFill) <- 'integer'
class(GdkFill) <- 'enums' 

GdkFunction<-c("copy" = 0,
	"invert" = 1,
	"xor" = 2,
	"clear" = 3,
	"and" = 4,
	"and-reverse" = 5,
	"and-invert" = 6,
	"noop" = 7,
	"or" = 8,
	"equiv" = 9,
	"or-reverse" = 10,
	"copy-invert" = 11,
	"or-invert" = 12,
	"nand" = 13,
	"nor" = 14,
	"set" = 15)
storage.mode(GdkFunction) <- 'integer'
class(GdkFunction) <- 'enums' 

GdkJoinStyle<-c("miter" = 0,
	"round" = 1,
	"bevel" = 2)
storage.mode(GdkJoinStyle) <- 'integer'
class(GdkJoinStyle) <- 'enums' 

GdkLineStyle<-c("solid" = 0,
	"on-off-dash" = 1,
	"double-dash" = 2)
storage.mode(GdkLineStyle) <- 'integer'
class(GdkLineStyle) <- 'enums' 

GdkSubwindowMode<-c("clip-by-children" = 0,
	"include-inferiors" = 1)
storage.mode(GdkSubwindowMode) <- 'integer'
class(GdkSubwindowMode) <- 'enums' 

GdkImageType<-c("normal" = 0,
	"shared" = 1,
	"fastest" = 2)
storage.mode(GdkImageType) <- 'integer'
class(GdkImageType) <- 'enums' 

GdkExtensionMode<-c("none" = 0,
	"all" = 1,
	"cursor" = 2)
storage.mode(GdkExtensionMode) <- 'integer'
class(GdkExtensionMode) <- 'enums' 

GdkInputSource<-c("mouse" = 0,
	"pen" = 1,
	"eraser" = 2,
	"cursor" = 3)
storage.mode(GdkInputSource) <- 'integer'
class(GdkInputSource) <- 'enums' 

GdkInputMode<-c("disabled" = 0,
	"screen" = 1,
	"window" = 2)
storage.mode(GdkInputMode) <- 'integer'
class(GdkInputMode) <- 'enums' 

GdkAxisUse<-c("ignore" = 0,
	"x" = 1,
	"y" = 2,
	"pressure" = 3,
	"xtilt" = 4,
	"ytilt" = 5,
	"wheel" = 6,
	"last" = 7)
storage.mode(GdkAxisUse) <- 'integer'
class(GdkAxisUse) <- 'enums' 

GdkPropMode<-c("replace" = 0,
	"prepend" = 1,
	"append" = 2)
storage.mode(GdkPropMode) <- 'integer'
class(GdkPropMode) <- 'enums' 

GdkFillRule<-c("even-odd-rule" = 0,
	"winding-rule" = 1)
storage.mode(GdkFillRule) <- 'integer'
class(GdkFillRule) <- 'enums' 

GdkOverlapType<-c("in" = 0,
	"out" = 1,
	"part" = 2)
storage.mode(GdkOverlapType) <- 'integer'
class(GdkOverlapType) <- 'enums' 

GdkRgbDither<-c("none" = 0,
	"normal" = 1,
	"max" = 2)
storage.mode(GdkRgbDither) <- 'integer'
class(GdkRgbDither) <- 'enums' 

GdkByteOrder<-c("lsb-first" = 0,
	"msb-first" = 1)
storage.mode(GdkByteOrder) <- 'integer'
class(GdkByteOrder) <- 'enums' 

GdkGrabStatus<-c("success" = 0,
	"already-grabbed" = 1,
	"invalid-time" = 2,
	"not-viewable" = 3,
	"frozen" = 4)
storage.mode(GdkGrabStatus) <- 'integer'
class(GdkGrabStatus) <- 'enums' 

GdkVisualType<-c("static-gray" = 0,
	"grayscale" = 1,
	"static-color" = 2,
	"pseudo-color" = 3,
	"true-color" = 4,
	"direct-color" = 5)
storage.mode(GdkVisualType) <- 'integer'
class(GdkVisualType) <- 'enums' 

GdkWindowClass<-c("output" = 0,
	"only" = 1)
storage.mode(GdkWindowClass) <- 'integer'
class(GdkWindowClass) <- 'enums' 

GdkWindowType<-c("root" = 0,
	"toplevel" = 1,
	"child" = 2,
	"dialog" = 3,
	"temp" = 4,
	"foreign" = 5)
storage.mode(GdkWindowType) <- 'integer'
class(GdkWindowType) <- 'enums' 

GdkWindowTypeHint<-c("normal" = 0,
	"dialog" = 1,
	"menu" = 2,
	"toolbar" = 3,
	"splashscreen" = 4,
	"utility" = 5,
	"dock" = 6,
	"desktop" = 7)
storage.mode(GdkWindowTypeHint) <- 'integer'
class(GdkWindowTypeHint) <- 'enums' 

GdkGravity<-c("north-west" = 0,
	"north" = 1,
	"north-east" = 2,
	"west" = 3,
	"center" = 4,
	"east" = 5,
	"south-west" = 6,
	"south" = 7,
	"south-east" = 8,
	"static" = 9)
storage.mode(GdkGravity) <- 'integer'
class(GdkGravity) <- 'enums' 

GdkWindowEdge<-c("north-west" = 0,
	"north" = 1,
	"north-east" = 2,
	"west" = 3,
	"east" = 4,
	"south-west" = 5,
	"south" = 6,
	"south-east" = 7)
storage.mode(GdkWindowEdge) <- 'integer'
class(GdkWindowEdge) <- 'enums' 

GdkPixbufAlphaMode<-c("bilevel" = 0,
	"full" = 1)
storage.mode(GdkPixbufAlphaMode) <- 'integer'
class(GdkPixbufAlphaMode) <- 'enums' 

GdkColorspace<-c("b" = 0)
storage.mode(GdkColorspace) <- 'integer'
class(GdkColorspace) <- 'enums' 

GdkPixbufError<-c("corrupt-image" = 0,
	"insufficient-memory" = 1,
	"bad-option-value" = 2,
	"unknown-type" = 3,
	"unsupported-operation" = 4,
	"failed" = 5)
storage.mode(GdkPixbufError) <- 'integer'
class(GdkPixbufError) <- 'enums' 

GdkPixbufRotation<-c("none" = 0,
	"counterclockwise" = 90,
	"upsidedown" = 180,
	"clockwise" = 270)
storage.mode(GdkPixbufRotation) <- 'integer'
class(GdkPixbufRotation) <- 'enums' 

GdkInterpType<-c("nearest" = 0,
	"tiles" = 1,
	"bilinear" = 2,
	"hyper" = 3)
storage.mode(GdkInterpType) <- 'integer'
class(GdkInterpType) <- 'enums' 

GdkOwnerChange<-c("new-owner" = 0,
	"destroy" = 1,
	"close" = 2)
storage.mode(GdkOwnerChange) <- 'integer'
class(GdkOwnerChange) <- 'enums' 

GdkDragAction<-c("default" = 1,
	"copy" = 2,
	"move" = 4,
	"link" = 8,
	"private" = 16,
	"ask" = 32)
storage.mode(GdkDragAction) <- 'numeric'
class(GdkDragAction) <- 'flags' 

GdkEventMask<-c("exposure-mask" = 2,
	"pointer-motion-mask" = 4,
	"pointer-motion-hint-mask" = 8,
	"button-motion-mask" = 16,
	"button1-motion-mask" = 32,
	"button2-motion-mask" = 64,
	"button3-motion-mask" = 128,
	"button-press-mask" = 256,
	"button-release-mask" = 512,
	"key-press-mask" = 1024,
	"key-release-mask" = 2048,
	"enter-notify-mask" = 4086,
	"leave-notify-mask" = 8192,
	"focus-change-mask" = 16384,
	"structure-mask" = 32768,
	"property-change-mask" = 65536,
	"visibility-notify-mask" = 131072,
	"proximity-in-mask" = 262144,
	"proximity-out-mask" = 524288,
	"substructure-mask" = 1048576,
	"scroll-mask" = 2097152,
	"all-events-mask" = 4194302)
storage.mode(GdkEventMask) <- 'numeric'
class(GdkEventMask) <- 'flags' 

GdkWindowState<-c("withdrawn" = 1,
	"iconified" = 2,
	"maximized" = 4,
	"sticky" = 8,
	"fullscreen" = 16,
	"above" = 32,
	"below" = 64)
storage.mode(GdkWindowState) <- 'numeric'
class(GdkWindowState) <- 'flags' 

GdkModifierType<-c("shift-mask" = 1,
	"lock-mask" = 2,
	"control-mask" = 4,
	"mod1-mask" = 8,
	"mod2-mask" = 16,
	"mod3-mask" = 32,
	"mod4-mask" = 64,
	"mod5-mask" = 128,
	"button1-mask" = 256,
	"button2-mask" = 512,
	"button3-mask" = 1024,
	"button4-mask" = 2048,
	"button5-mask" = 4096,
	"release-mask" = 8192,
	"modifier-mask" = 16384)
storage.mode(GdkModifierType) <- 'numeric'
class(GdkModifierType) <- 'flags' 

GdkWindowAttributesType<-c("title" = 1,
	"x" = 2,
	"y" = 4,
	"cursor" = 8,
	"colormap" = 16,
	"visual" = 32,
	"wmclass" = 64,
	"noredir" = 128)
storage.mode(GdkWindowAttributesType) <- 'numeric'
class(GdkWindowAttributesType) <- 'flags' 

GdkWindowHints<-c("pos" = 1,
	"min-size" = 2,
	"max-size" = 4,
	"base-size" = 8,
	"aspect" = 16,
	"resize-inc" = 32,
	"win-gravity" = 64,
	"user-pos" = 128,
	"user-size" = 256)
storage.mode(GdkWindowHints) <- 'numeric'
class(GdkWindowHints) <- 'flags' 

GdkWMDecoration<-c("all" = 1,
	"border" = 2,
	"resizeh" = 4,
	"title" = 8,
	"menu" = 16,
	"minimize" = 32,
	"maximize" = 64)
storage.mode(GdkWMDecoration) <- 'numeric'
class(GdkWMDecoration) <- 'flags' 

GdkWMFunction<-c("all" = 1,
	"resize" = 2,
	"move" = 4,
	"minimize" = 8,
	"maximize" = 16,
	"close" = 32)
storage.mode(GdkWMFunction) <- 'numeric'
class(GdkWMFunction) <- 'flags' 

