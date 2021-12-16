#' Convert an (i, j) index to a linear index.
#'
#' Converts an (i, j) index to a linear index. Converts the double index of
#' a square matrix to the corresponding linear one. This is column-major
#' as it is default in R.
#'
#' @param i i index, i.e. row position; indexing starts at 1.
#' @param j j index, i.e. column position; indexing starts at 1.
#' @param n size of the square matrix.
#' @return Linear position.
#' @keywords internal
ij2k <- function(i, j, n) (j - 1) * n + i
