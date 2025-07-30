#' Coordinates as matrix
#'
#' @param x Coordinates specified as a two-column data frame, matrix, or list 
#'   of equal-length vectors
#'
#' @keywords internal
coord_matrix <- function(x) {
  if (is.list(x)) x <- data.frame(x) # also ensures equal length
  stopifnot(is.numeric(x[,1]) & is.numeric(x[,2]))
  stopifnot(all(!is.na(c(x[,1], x[,2]))))
  matrix(c(x[,1], x[,2]), ncol = 2)
}
