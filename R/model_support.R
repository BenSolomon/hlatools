#' contr.dummy from kknn
#'
#' @param n test
#' @param contrasts test
#'
#' @return
#'
# contr.dummy <- function (n, contrasts = TRUE)
# {
#   if (length(n) <= 1) {
#     if (is.numeric(n) && length(n) == 1 && n > 1)
#       levels <- 1:n
#     else stop("contrasts are not defined for 0 degrees of freedom")
#   }
#   else levels <- n
#   lenglev <- length(levels)
#   cont <- array(0, c(lenglev, lenglev), list(levels, levels))
#   cont[col(cont) == row(cont)] <- 1
#   cont
# }
