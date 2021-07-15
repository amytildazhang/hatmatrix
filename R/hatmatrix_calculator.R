

# HatMatrix -----------------------------------------------------------------

#' Create a new [`HatMatrixCalculator`] object.
#'
#' @inheritParams hatmatrix
#' @return A [`HatMatrixCalculator`] object.
#' @export
#'
#' @examples
hatmatrix_calculator <- function(Xf, Phi_inv, Xr = NULL, Sigma_inv = NULL) {
    mc <- match.call()
    m <- match(c("Xf", "Xr", "Phi_inv", "Sigma_inv"), table = names(mc), nomatch = 0L)
    mc <- mc[c(1L, m)]
    mc[[1L]] <- quote(HatMatrix$new)
    eval(mc, parent.frame())
}




#' R6 Class representing a HatMatrixCalculator
#'
#' @description A `HatMatrixCalculator` object is an [R6][R6::R6Class] object created
#'   by the [`hatmatrix_calculator()`] function. Internally, it stores and validates inputs to
#' [`hatmatrix()`] or [`hatmatrix_from()`] and contains the methods to calculate the hat matrix,
#' sub-matrix of the hat matrix,
#' or the hat matrix based on a subset of the data.
#'
#'
#'
#'
#' @export
HatMatrixCalculator <- R6::R6Class(
    "HatMatrixCalculator",
    private = list(
        .n = 1,
        .i = NULL,
        .j = NULL,
        .subset = NULL,
        .Xf = matrix(),
        .Xr = matrix(),
        .Sigma_inv = matrix(),
        .Phi_inv = c(),
        .W = matrix(),
        .V = matrix(),

        .save_values = function(Xf, Xr, Sigma_inv, Phi_inv) {
            private$.Xf <- Xf
            private$.Xr <- Xr
            private$.n <- nrow(Xf)
            private$.Sigma_inv <- Sigma_inv
            private$.Phi_inv <- Phi_inv
        },
        .subset_idx = function(subset) {
            if (is.logical(subset)) {
                return(which(subset))
            } else {
                return(as.integer(subset))
            }
        }
    ),
    public = list(
        initialize = function(Xf, Phi_inv, Xr = NULL, Sigma_inv = NULL) {
            HatMatrixValidator$new(Xf, Phi_inv, Xr, Sigma_inv, n = nrow(Xf))
            private$.save_values(Xf, Xr, Sigma_inv, Phi_inv)
        }
    ),
    active = list(
        #' @field number of data points (rows in \mjeqn{X_f}{X_f})
        n = pryr::unenclose(read_only(.n)),
        #' @field rows selected from \mjeqn{W}{W} in the most recent call to [$calc()]
        i =  pryr::unenclose(read_only(.i)),
        #' @field columns selected from \mjeqn{W}{W} in the most recent call to [$calc()]
        j =  pryr::unenclose(read_only(.j)),
        #' @field supplied design matrix for fixed effects
        Xf =  pryr::unenclose(read_only(.Xf)),
        #' @field supplied design matrix for random effects
        Xr =  pryr::unenclose(read_only(.Xr)),
        #' @field supplied \mjeqn{\Sigma^{-1}}{Sigma^-1}
        Sigma_inv =  pryr::unenclose(read_only(.Sigma_inv)),
        #' @field supplied \mjeqn{\Phi^{-1}}{Phi^-1}
        Phi_inv =  pryr::unenclose(read_only(.Phi_inv)),
        #' @field hat matrix from most recent call to [$calc()]
        W =  pryr::unenclose(read_only(.W)),
        #' @field V matrix from most recent call to [$calc()]
        V = pryr::unenclose(read_only(.V))
    )
)


#' Title
#'

#' @export
#'
#' @examples
expanded_sigma <- function() {
    fixef_idx <- 1:ncol(self$Xf)
    Sigma_inv <- self$Sigma_inv

    P <- ncol(Sigma_inv) + ncol(self$Xf)
    sigi <- matrix(0, ncol = P, nrow = P)
    sigi[-fixef_idx, -fixef_idx] <- Sigma_inv
    sigi
}
HatMatrixCalculator$set("public", name = "expanded_sigma", value = expanded_sigma)

#' Title
#'
#' @description
#' Convenience function to create a diagonal matrix of noise precision.
#' Created for the edge case where only one data point is used.
#' @param idx The indices of precision values to select from the stored Phi^-1.
#'
#' @return
#' @export
#'
#' @examples
diag_phi <- function(idx) {
    diag(self$Phi_inv[idx], nrow = length(idx), ncol = length(idx))
}
HatMatrixCalculator$set("public", name = "diag_phi", value = diag_phi)



#' Title
#'
#' @description \loadmathjax
#' @inheritParams hatmatrix
#' @return
#' @export
#'
#' @examples
calc <- function(i = 1:self$n, j = 1:self$n, return_V = FALSE, subset = NULL) {
    i <- as.integer(i); j <- as.integer(j)

    HatMatrixValidator$new(i = i,j = j, n = self$n)
    if (!is.null(subset)) {
        HatMatrixValidator$new(subset = subset, j = j, n = self$n)
        subset <- private$.subset_idx(subset)
    } else {
        subset <- 1:self$n
    }
    private$.i <- i
    private$.j <- j
    private$.subset <- subset

    X <- cbind(self$Xf, self$Xr)
    X_i <- X[i, , drop = F]
    X_j <- X[j, , drop = F]
    X_tr <- X[subset, , drop = F]

    Sigma_i <- self$expanded_sigma()

    V <- solve(t(X_tr) %*% self$diag_phi(subset) %*% X_tr + Sigma_i)
    private$.V <- V
    if (return_V) {
        return(V)
    }
    W <-  X_i %*% V %*% t(X_j) %*% self$diag_phi(j)
    private$.W <- W
    HatMatrixValidator$new(W = W)
    W
}
HatMatrixCalculator$set("public", name = "calc", value = calc)


