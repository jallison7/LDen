#' Hodder and Okell's dispersion ratio and index of association
#' 
#' `association()` calculates Hodder and Okell's (1978) 'A' index of association 
#' between two sets of points, which is the product of the ratios of the average 
#' distance between points within each set to the average distance between the 
#' sets.
#'
#' `dispersion()` calculates the ratio of the average distance between points in
#' two sets (i.e. the upper term of A).
#' 
#' @inheritParams coord_matrix
#' @param y Same as `x`.
#' 
#' @return 
#' `association()` returns Hodder and Okell's A index of association between 
#' `x` and `y`. For the interpretation of this value see Orton (1980, p. 154):
#' "When the distribution of types A and B are randomly mingled, with no
#' association or dissociation, A will have a value of about 1: closely
#' packed but separate distributions have low values of A, and if As and Bs
#' tend to occur together, then A will be greater than 1."
#'
#' `dispersion()` returns the dispersion ratio, which is 1 for identical
#' distributions. Note that this is not symmetrical, i.e. `dispersion(x, y)`
#' does not necessarily equal `dispersion(y, x)`.
#'
#' @references
#' * Hodder, I. and E. Okell. 1978. A new method for assessing the association 
#'   between dismbutions of points in archaeology. In I. Hodder (ed), 
#'   *Simulation Studies in Archaeology*, p. 97-107. Cambridge: Cambridge 
#'   University Press.
#' * Orton, C. 1980. *Mathematics in archaeology*. London: Collins.
#' 
#' @export
#' @examples
#' ceramics <- AZ_A1020_BLM[AZ_A1020_BLM$type == "ceramics",]
#' stone_tools <- AZ_A1020_BLM[AZ_A1020_BLM$type == "stone tool",]
#'
#' association(ceramics, stone_tools)
#' dispersion(ceramics, stone_tools)
association <- function(x, y) {
  coordx <- coord_matrix(x)
  coordy <- coord_matrix(y)

  dxx <- mean(proxy::dist(coordx, method = "euclidean"))
  dyy <- mean(proxy::dist(coordy, method = "euclidean"))
  dxy <- mean(proxy::dist(coordx, coordy, method = "euclidean"))

  (dxx / dxy) * (dyy / dxy)
}

#' @export
#' @rdname association
dispersion <- function(x, y) {
  coordx <- coord_matrix(x)
  coordy <- coord_matrix(y)

  dxx <- mean(proxy::dist(coordx, method = "euclidean"))
  dyy <- mean(proxy::dist(coordy, method = "euclidean"))

  dxx / dyy
}
