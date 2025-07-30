#' Global point density
#'
#' Calculates the 'global' density of a vector of point coordinates in relation 
#' to the total size of the window or study region (e.g. the site or trench).
#'
#' @inheritParams coord_matrix
#' @param area Total area of the window. If `NULL`, inferred as the area of the
#'   bounding box of the coordinates specified by `x` and `y`.
#'
#' @return Number of points divided by area
#' @export
#'
#' @examples
#' global_density(AZ_A1020_BLM, area = 2409)
#' global_density(AZ_A1020_BLM) # area estimated from bounding box
global_density <- function(x, area = NULL) {
  coords <- coord_matrix(x)
  if (is.null(area)) area <- area_bbox(coords) # TODO: message with area?
  nrow(x) / area
}

#' Calculate area of the bounding box of a set of coordinates
#'
#' @noRd
#' @keywords internal
area_bbox <- function(x) {
  (max(x[,1]) - min(x[,1])) * (max(x[,2]) - min(x[,2]))
}
