#' AZ_A1020_BLM Test Data
#'
#' This data set ontains point location and type data for 358 "diagnostic"
#' surface artifacts from an Ancestral Pueblo site on the
#' Shivwits Plateau in northwestern Arizona. In the Arizona site
#' numbering system, the site number is formally "AZ A:10:20 (BLM)".
#'
#' This is real data, although somewhat simplified here.
#' There are five 'types' (three of them rare), and the two most common
#' 'types' group artifacts that have more detailed identifications
#' ("ceramics" = either painted potsherds or sherds from vessel rims;
#' "stone tool" = any flaked stone tool except for projectile points).
#'
#' @format ## `AZ_A1020_BLM`
#' A data frame with 358 rows and 3 columns:
#' \describe{
#'   \item{x}{Easting of location in EPSG ...}
#'   \item{y}{Northing of location in EPSG ...}
#'   \item{type}{Type of Evidence at location, see description}
#' }
"AZ_A1020_BLM"
