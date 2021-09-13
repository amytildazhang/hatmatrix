
#' Calculate hat matrix from model matrix and variance parameters
#' @description
#' \loadmathjax
#'
#'
#' @template xf-arg
#' @template phi-arg
#' @template xr-arg
#' @template sigma-arg
#' @template ijs-args
#'
#' @return Numeric matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' library(hatmatrix)
#' Xf <- model.matrix(mpg ~ disp, data = mtcars)
#' Xr <- model.matrix(mpg ~ -1 + factor(cyl), data = mtcars)
#' Sigma_inv <- diag(rep(3, ncol(Xr)))
#' Phi_inv <- rep(2, nrow(Xf))
#' hatmatrix(Xf, Phi_inv, Xr, Sigma_inv)
#' hatmatrix(Xf, Phi_inv, Xr, Sigma_inv, i = 1:4) # returns 4 x N matrix
#' hatmatrix(Xf, Phi_inv, Xr, Sigma_inv, i = 1:4, subset = 1:10) # returns 4 x 10 matrix
#' hatmatrix(Xf, Phi_inv, Xr, Sigma_inv, i = 1:4, subset = 1:10, j = 11) # throws an error because
#' hatmatrix(Xf, Phi_inv, Xr, Sigma_inv, i = 1:4, subset = 1:10, j = 5:8) # returns 4 x 4 matrix
#' }
hatmatrix <- function(Xf, Phi_inv, Xr = NULL, Sigma_inv = NULL,
                      i = 1:nrow(Xf), j = subset, subset = 1:nrow(Xf)) {
  hm <- hatmatrix_calculator(Xf, Phi_inv, Xr, Sigma_inv)
  hm$calc(i = i, j = j, subset = subset)
}







HatMatrixExtractor <- R6::R6Class(
    "HatMatrixExtractor",
    private = list(),
    public = list(
        initialize = function(obj) {
            #
        },
        extract_xfr = function(obj) {

        },
        extract_phi = function(obj) {

        },
        extract_sigma = function(obj) {

        },
        calc_hatmatrix = function(obj) {
            X <- self$extract_xfr(obj)
            Xf <- X$Xf; Xr <- X$Xr;
            Phi_inv <- self$extract_phi(obj)
            Sigma_inv <- self$extract_sigma(obj)

            hatmatrix(Xf, Phi_inv, Xr, Sigma_inv)

        }
    )

)

#' Generic for getting the hat matrix from a fitted model object
#'
#' @param x Fitted model object
#'
#' @return The hat matrix
#' @export
#' @family {hat matrix helpers}
#' @examples
hatmatrix_from <- function(x) {
  UseMethod("hatmatrix_from")
}




hatmatrix_from.formula <- function(x, ftype, Sigma_inv, Phi_inv, i, j, ...) {
  # call build_x
}

#
# hatmatrix_from.hm_model_matrix <- function(x, Sigma_inv, Phi_inv, i, j) {
#   validate_hm_model_matrix(x)
#   Xf <- hmm_to_xf(x)
#   Xr <- hmm_to_xr()
#   hatmatrix(Xf, Xr, Sigma_inv, Phi_inv, i, j)
# }


LMHatMatrixExtractor <- R6::R6Class(
    "LMHatMatrixExtractor",
    inherit = HatMatrixExtractor,
    private = list(),
    public = list(
        initialize = function(obj) {
            #
        },
        extract_xfr = function(obj) {

            # list(Xr = NULL, Xf = )
        },
        extract_phi = function(obj) {

        },
        extract_sigma = function(obj) {
            return(NULL)
        }
    )
)


#' Title
#'
#' @param obj
#'
#' @return
#'
#' @examples
#' fit <- lm(mpg ~ wt, data = mtcars)
hatmatrix_from.lm <- function(obj) {
    ext <- LMHatMatrixExttractor$new(obj)
    ext$calc_hatmatrix()
}


#' Title
#'
#' @param mod
#'
#' @return
#'
#' @examples
hatmatrix_from.lmer <- function(mod) {
  # TODO
}
