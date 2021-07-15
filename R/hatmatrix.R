
#' Calculate hat matrix from model matrix and variance parameters
#' @description
#' \loadmathjax
#'
#' @details
#' For the regression model with response vector \mjeqn{Y \in R^{N}}{Y in R^N} and
#' vector of coefficients \mjeqn{\beta \in R^P}{} such that
#' \mjdeqn{Y| \beta, \Phi \sim N(X_f\beta_f + X_r\beta_r, \Phi),}{}
#' \mjdeqn{\beta_f \sim N(0, C), \quad \beta_r \sim N(0, \Sigma),}{}
#' where \mjeqn{\beta_f \in R^{P_f}}{} is the vector of fixed effect coefficients,
#' \mjeqn{\beta_r \in R^{P_r}}{} the vector of random effect coefficients,
#' \mjeqn{X_f}{} and \mjeqn{X_r}{} are likewise the design matrix corresponding to the
#' fixed and random effects, respectively, and \mjeqn{\Sigma}{} is the random effect
#' covariance matrix. The hat matrix, \mjeqn{W}{}, is then
#'
#'
#'\mjdeqn{W = XVX'\Phi^{-1}}{},
#'\mjdeqn{V = \left(X'\Phi^{-1}X + \begin{bmatrix} C^{-1} & 0 \\
#'0 & \Sigma^{-1}\end{bmatrix}\right)^{-1}}{}.
#'\mjeqn{C^{-1}}{C^-1} is taken to be the matrix of 0s.
#'
#' @param Xf  \mjeqn{X_f \text{ in `Details'.}}{}
#'            The \mjeqn{{N \times P_1}}{} model matrix for the fixed effects,
#'            where \mjeqn{N}{} is the amount of data and \mjeqn{P_1}{} the dimension
#'            of fixed effect coefficients. Must be full-rank.
#' @param Phi_inv \mjeqn{\Phi^{-1} \text{ in `Details'.}}{} The \mjeqn{N}{}-length numeric
#'            vector with positive values
#'            corresponding to the precision of \mjeqn{Y | X\beta}{}.
#' @param Xr \mjeqn{X_r \text{ in `Details'.}}{} \mjeqn{N \times P_2}{} model matrix for the random effects,
#'            where \mjeqn{P_2}{} is the dimension of random effect coefficients.
#'            Defaults to `NULL`.
#'
#' @param Sigma_inv \mjeqn{\Sigma^{-1} \text{ in `Details'.}}{} Sigma_inv \mjeqn{P_2 \times P_2}{} precision matrix for the random
#'            effects (\mjeqn{\Sigma^{-1}}{}). Must be positive-definite.
#'            Defaults to `NULL`.
#'
#' @param i Integer vector. Selects rows of \mjeqn{W}{} to output.
#'          Entries must be between \mjeqn{1}{} and \mjeqn{N}{}.
#'          `thm$calc(i = 1)` is equivalent to `thm$calc()[1, , drop = F]`.
#'          Defaults to all rows.
#'
#' @param j Integer vector. Selects columns of \mjeqn{W}{} to output.
#'          Entries must be between \mjeqn{1}{} and \mjeqn{N}{},
#'          and are restricted to values present in `subset`. Defaults to `subset`.
#'          `thm$calc(j = 1)` is equivalent to `thm$calc()[, 1, drop = F]`.
#' @param subset Integer vector or logical vector with length \mjeqn{N}{}.
#'          Selects rows of \mjeqn{X}{} to use for calculation of \mjeqn{V}{} (see `Details'), e.g.
#'          if in the case of training data, `subset` specifies the indices of the
#'          training data. Defaults to all rows.
#'
#'
#' @return Numeric matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' library(hatmatrix)
#' Xf <- model.matrix(mpg ~ disp, data = mtcars)
#' Xr <- model.matrix(mpg~-1 + factor(cyl), data = mtcars)
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


hatmatrix_from.hm_model_matrix <- function(x, Sigma_inv, Phi_inv, i, j) {
    validate_hm_model_matrix(x)
    Xf <- hmm_to_xf(x)
    Xr <- hmm_to_xr()
    hatmatrix(Xf, Xr, Sigma_inv, Phi_inv, i, j)
}


#' Title
#'
#' @param x
#'
#' @return
#'
#' @examples
#' fit <- lm(mpg ~ wt, data = mtcars)
hatmatrix_from.lm <- function(mod) {
    # validate lm object
    #TODO
}


#' Title
#'
#' @param mod
#'
#' @return
#'
#' @examples
hatmatrix_from.lmer <- function(mod) {
    #TODO

}



