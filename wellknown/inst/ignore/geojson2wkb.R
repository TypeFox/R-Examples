# '\x00': The first byte of any WKB string. Indicates big endian byte
# ordering for the data.
BIG_ENDIAN = '\x00'
# '\x01': The first byte of any WKB string. Indicates little endian byte
# ordering for the data.
LITTLE_ENDIAN = '\x01'

# WKB_2D = {
#   'Point': b'\x00\x00\x00\x01',
#   'LineString': b'\x00\x00\x00\x02',
#   'Polygon': b'\x00\x00\x00\x03',
#   'MultiPoint': b'\x00\x00\x00\x04',
#   'MultiLineString': b'\x00\x00\x00\x05',
#   'MultiPolygon': b'\x00\x00\x00\x06',
#   'GeometryCollection': b'\x00\x00\x00\x07',
# }
#
# #: Mapping of GeoJSON geometry types to the "Z" 4-byte binary string
# #: representation for WKB. "Z" indicates that the geometry is 3-dimensional,
# #: with X, Y, and Z components.
# #: NOTE: Byte ordering is big endian.
# WKB_Z = {
#   'Point': b'\x00\x00\x03\xe9',
#   'LineString': b'\x00\x00\x03\xea',
#   'Polygon': b'\x00\x00\x03\xeb',
#   'MultiPoint': b'\x00\x00\x03\xec',
#   'MultiLineString': b'\x00\x00\x03\xed',
#   'MultiPolygon': b'\x00\x00\x03\xee',
#   'GeometryCollection': b'\x00\x00\x03\xef',
# }
#
# #: Mapping of GeoJSON geometry types to the "M" 4-byte binary string
# #: representation for WKB. "M" indicates that the geometry is 2-dimensional,
# #: with X, Y, and M ("Measure") components.
# #: NOTE: Byte ordering is big endian.
# WKB_M = {
#   'Point': b'\x00\x00\x07\xd1',
#   'LineString': b'\x00\x00\x07\xd2',
#   'Polygon': b'\x00\x00\x07\xd3',
#   'MultiPoint': b'\x00\x00\x07\xd4',
#   'MultiLineString': b'\x00\x00\x07\xd5',
#   'MultiPolygon': b'\x00\x00\x07\xd6',
#   'GeometryCollection': b'\x00\x00\x07\xd7',
# }
#
# # Mapping of GeoJSON geometry types to the "ZM" 4-byte binary string
# # representation for WKB. "ZM" indicates that the geometry is 4-dimensional,
# # with X, Y, Z, and M ("Measure") components.
# # NOTE: Byte ordering is big endian.
# WKB_ZM <- list(
#   'Point' = '\x00\x00\x0b\xb9',
#   'LineString' = '\x00\x00\x0b\xba',
#   'Polygon' = '\x00\x00\x0b\xbb',
#   'MultiPoint' = '\x00\x00\x0b\xbc',
#   'MultiLineString': = '\x00\x00\x0b\xbd',
#   'MultiPolygon' = '\x00\x00\x0b\xbe',
#   'GeometryCollection' = '\x00\x00\x0b\xbf',
# )
#
# #: Mapping of dimension types to maps of GeoJSON geometry type -> 4-byte binary
# #: string representation for WKB.
# wkb <- function(x) {
#   switch(x,
#          '2D' = WKB_2D,
#          'Z' = WKB_Z,
#          'M' = WKB_M,
#          'ZM' = WKB_ZM
#   )
# }
#
# int2dimlabel <- function(x) {
#   switch(x,
#          `2`='2D',
#          `3`='Z',
#          `4`='ZM')
# }
#
# #' Convert GeoJSON-like POINT object to WKT.
# #'
# #' @inheritParams geojson2wkt
# #' @keywords internal
# dump_point <- function(obj, fmt = 16){
#   coords <- obj$coordinates
#   sprintf('POINT (%s)', paste0(format(coords, nsmall = fmt), collapse = ""))
# }
#
#
# def _dump_point(obj, big_endian):
#   """
#     Dump a GeoJSON-like `dict` to a point WKB string.
#     :param dict obj:
#         GeoJson-like `dict` object.
#     :param bool big_endian:
#         If `True`, data values in the generated WKB will be represented using
#         big endian byte order. Else, little endian.
#     :returns:
#         A WKB binary string representing of the Point ``obj``.
#     """
#     coords = obj['coordinates']
#     num_dims = len(coords)
#
#     wkb_string, byte_fmt, _ = _header_bytefmt_byteorder(
#       'Point', num_dims, big_endian
#     )
#
#     wkb_string += struct.pack(byte_fmt, *coords)
#     return wkb_string
#
# def _header_bytefmt_byteorder(geom_type, num_dims, big_endian):
#     """
#     Utility function to get the WKB header (endian byte + type header), byte
#     format string, and byte order string.
#     """
#     dim = _INT_TO_DIM_LABEL.get(num_dims)
#     if dim is None:
#         pass  # TODO: raise
#
#     type_byte_str = _WKB[dim][geom_type]
#
#     if big_endian:
#         header = BIG_ENDIAN
#         byte_fmt = b'>'
#         byte_order = '>'
#     else:
#         header = LITTLE_ENDIAN
#         byte_fmt = b'<'
#         byte_order = '<'
#         # reverse the byte ordering for little endian
#         type_byte_str = type_byte_str[::-1]
#
#     header += type_byte_str
#     byte_fmt += b'd' * num_dims
#
#     return header, byte_fmt, byte_order
