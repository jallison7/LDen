#' Local Density Analysis
#'
#' `lda()` calculates the local density coefficient within and between
#' all the types in a dataset. The local density coefficient from
#' Type A to Type B is the mean density of points of Type B within
#' the specified radius of points of Type A, divided by the global
#' density of points of Type B (i.e., the number of points of
#' Type B divided by the total area)
#'
#' @inheritParams local_count
#' @inheritParams global_density
#' @param type Vector the same length as `x` giving the types or other 
#'   categorical variable to construct sets for pairwise comparison.
#' 
#' @return A data.frame containing the local density coefficients.
#'
#' @export
#'
#' @examples
#' lda(AZ_A1020_BLM, type = AZ_A1020_BLM$type, radius = 2, area = 2049)
lda <- function(x, type, radius, area = NULL) {
  coords <- coord_matrix(x)
  dat <- data.frame(x = coords[,1], y = coords[,2], type = type)
  dat <- split(dat, ~type)

  # For each combination of types...
  res <- lapply(dat, function(subdat, dat, r, a) {
    lapply(dat, function(daty, datx, r, a) {
      # Calculate two-sample local densities
      ldens <- local_density(datx, daty, radius = r)
      # Return mean density divided by global density of the second type
      mean(ldens) / global_density(daty, area = a)
    }, datx = subdat, r = r, a = a)
  }, dat = dat, r = radius, a = area)

  # Flatten to matrix
  res <- do.call(cbind, lapply(res, unlist))

  # Recalculate diagonal using local_count - 1 (removing the origin point from
  # its own neighbourhood)
  # TODO: questioning, see https://github.com/jallison7/LDen/issues/10
  diag(res) <- vapply(dat, function(subdat, r, a) {
    ldens <- local_density(subdat, radius = r)
    mean(ldens) / global_density(subdat, area = a)
  }, numeric(1), r = radius, a = area)

  # TODO: return as long data frame?
  res
}
