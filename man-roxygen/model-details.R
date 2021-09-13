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
#' \mjdeqn{W = XVX'\Phi^{-1}}{},
#' \mjdeqn{V = \left(X'\Phi^{-1}X + \begin{bmatrix} C^{-1} & 0 \\
#' 0 & \Sigma^{-1}\end{bmatrix}\right)^{-1}}{}.
#' \mjeqn{C^{-1}}{C^-1} is taken to be the matrix of 0s.
