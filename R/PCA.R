#' Title
#'
#' @param X a matrix
#' @param standardize a boolean declaring whether to standardize
#'
#' @return X a centred matrix
#' @export
#'
#' @examples X <- matrix(seq(1:9))
#' centre(X)
centre <- function(X, standardize = F) {
  if (standardize) {
    X <- apply(X, MARGIN  = 2, FUN = function(x) {x <- x - mean(x)
                                          return(x/sd(x)) })
  }

  else {
    X <- apply(X, MARGIN  = 2, FUN = function(x) x - mean(x))
  }
  return(X)
}
