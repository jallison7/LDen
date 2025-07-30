#' Local point counts and density
#'
#' For a vector of points specified by x and y coordinates, `local_count()`
#' calculates the number of points within the given radius. `local_density()`
#' returns the density of points, i.e. the count divided by the area of a circle
#' with the given radius.
#'
#' `local_count2()` and `local_density2()` are two-sample versions of the above,
#' calculating the count of density of points in y within the specified radius
#' of points in x.
#'
#' @inheritParams coord_matrix
#' @param y For the two-sample case, coordinates of points to count, in the same
#'   format as `x`. Set `NULL` (the default) for the one-sample case.
#' @param radius Size of the neighbourhood
#'
#' @return Numeric vector the same length as `x` with the number or density of 
#' points within the specified radius of each point. In the one sample case,
#' the count includes the origin point (and thus is always >= 1) but the 
#' density does not.
#'
#' @export
#'
#' @examples
#' local_count(AZ_A1020_BLM, radius = 2)
#' local_density(AZ_A1020_BLM, radius = 2)
local_count <- function(x, y = NULL, radius) {
  stopifnot(is.numeric(radius))

  # Distance matrix for the one-sample case, cross-distance matrix for the
  # two-sample case
  if (is.null(y)) {
    dist_matrix <- proxy::dist(coord_matrix(x), method = "euclidean")
  }
  else {
    dist_matrix <- proxy::dist(coord_matrix(x), coord_matrix(y), 
                               method = "euclidean")
  }

  rowSums(as.matrix(dist_matrix) <= radius)
}

#' @rdname local_count
#' @export
local_density <- function(x, y = NULL, radius) {
  count <- local_count(x, y, radius)

  # Don't count origin point in density for the one-sample case
  # See https://github.com/jallison7/LDen/issues/10
  if (is.null(y)) count <- count - 1

  count / (pi * radius^2)
}

