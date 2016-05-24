gdkColormap <- gdkColormapNew

gdkDisplay <- gdkDisplayOpen

gdkDragContext <- gdkDragContextNew

gdkWindow <- gdkWindowNew

gdkPixmap <-
function(drawable, data, transparent.color, height, width, depth, filename, colormap, fg, bg)
{
  if (!missing(filename)) {
    if (!missing(colormap)) {
      gdkPixmapColormapCreateFromXpm(drawable, colormap, transparent.color, filename)
    }
    else {
      gdkPixmapCreateFromXpm(drawable, transparent.color, filename)
    }
  }
  else {
    if (!missing(width)) {
      if (!missing(depth)) {
        gdkPixmapNew(drawable, width, height, depth)
      }
      else {
        gdkBitmapCreateFromData(drawable, data, width, height)
      }
    }
    else {
      if (!missing(height)) {
        gdkPixmapCreateFromData(drawable, data, height, depth, fg, bg)
      }
      else {
        if (!missing(colormap)) {
          gdkPixmapColormapCreateFromXpmD(drawable, colormap, transparent.color, data)
        }
        else {
          gdkPixmapCreateFromXpmD(drawable, transparent.color, data)
        }
      }
    }
  }
}

gdkGC <- gdkGCNew

gdkImage <- gdkImageNew

gdkPixbuf <-
function(width = -1, height = -1, filename, colorspace, has.alpha, bits.per.sample, preserve.aspect.ratio = 1, data, stream, cancellable = NULL, rowstride, .errwarn = TRUE)
{
  if (!missing(width)) {
    if (!missing(colorspace)) {
      if (!missing(data)) {
        gdkPixbufNewFromData(data, colorspace, has.alpha, bits.per.sample, width, height, rowstride)
      }
      else {
        gdkPixbufNew(colorspace, has.alpha, bits.per.sample, width, height)
      }
    }
    else {
      if (!missing(filename)) {
        if (!missing(preserve.aspect.ratio)) {
          gdkPixbufNewFromFileAtScale(filename, width, height, preserve.aspect.ratio, .errwarn)
        }
        else {
          gdkPixbufNewFromFileAtSize(filename, width, height, .errwarn)
        }
      }
      else {
        gdkPixbufNewFromStreamAtScale(stream, width, height, preserve.aspect.ratio, cancellable, .errwarn)
      }
    }
  }
  else {
    if (!missing(filename)) {
      gdkPixbufNewFromFile(filename, .errwarn)
    }
    else {
      if (!missing(.errwarn)) {
        gdkPixbufNewFromStream(stream, cancellable, .errwarn)
      }
      else {
        gdkPixbufNewFromXpmData(data)
      }
    }
  }
}

gdkPixbufAnimation <- gdkPixbufAnimationNewFromFile

gdkPixbufSimpleAnim <- gdkPixbufSimpleAnimNew

gdkPixbufLoader <-
function(image.type, mime.type, .errwarn = TRUE)
{
  if (!missing(image.type)) {
    gdkPixbufLoaderNewWithType(image.type, .errwarn)
  }
  else {
    if (!missing(.errwarn)) {
      gdkPixbufLoaderNewWithMimeType(mime.type, .errwarn)
    }
    else {
      gdkPixbufLoaderNew()
    }
  }
}

gdkPangoRenderer <- gdkPangoRendererNew

gdkAppLaunchContext <- gdkAppLaunchContextNew

