Validator <- R6::R6Class(
  "Validator",
  private = list(
    .print_error = function(message, e, ...) {
      if (is.character(e)) {
        stop(sprintf(message, e, ...))
      }
    }
  )
)


HatMatrixValidator <- R6::R6Class(
  "HatMatrixValidator",
  inherit = Validator,
  public = list(
    initialize = function(Xf = NULL, Phi_inv = NULL, Xr = NULL, Sigma_inv = NULL,
                          P_2 = NULL, n = NULL, i = NULL, j = NULL, subset = NULL,
                          W = NULL) {
      if (!is.null(Xf)) self$validate_x(Xf, Xr)
      if (!is.null(Sigma_inv)) self$validate_sigma(Sigma_inv, P_2)
      if (!is.null(Phi_inv)) self$validate_phi(Phi_inv, n)
      if (!is.null(i)) self$validate_ij(i, n, "i")
      if (!is.null(j)) self$validate_ij(j, n, "j")
      if (!is.null(subset)) self$validate_subset(subset, j, n)
      if (!is.null(W)) self$validate_w(W)

      invisible(TRUE)
    },
    validate_subset = function(subset, j, n) {
      checkmate::assert(
        checkmate::check_logical(subset, len = n, any.missing = FALSE),
        checkmate::check_integer(subset, min.len = 1, any.missing = F, lower = 1, upper = n)
      )

      if (is.logical(subset)) {
        subset <- which(subset)
      }
      if (!all(j %in% subset)) {
        dif <- setdiff(j, subset)
        stop(sprintf(
          "Input error in subset or j vector: Column indices (j) are not present in training data (subset). Extra elements in j are %s",
          dif
        ))
      }
    },
    validate_ij = function(vec, n, name) {
      msg <- checkmate::check_integer(vec, lower = 1, upper = n)
      private$.print_error(
        "Input error for %s. Based on supplied Xf with %s rows, entries must be between 1 and %s.",
        ifelse(is.character(msg), sprintf("vector %s: %s", name, msg), TRUE),
        n, n
      )
    },
    validate_phi = function(phi, n) {
      msg <- checkmate::check_numeric(phi,
        lower = 0,
        any.missing = F, len = n,
        finite = TRUE
      )
      private$.print_error(
        "Input error for Phi^-1: %s. Based on supplied Xf with %s rows, Phi^-1 must be a %s-length vector with finite, positive values.",
        msg, n, n
      )
    },
    validate_sigma = function(sigma, P_2) {
      msg <- checkmate::check_matrix(sigma,
        mode = "numeric",
        nrows = P_2, ncols = P_2
      )
      private$.print_error(
        "Input erorr for Sigma^-1: %s. Based on supplied Xr with %s columns, numeric matrix Sigma^-1 must have %s rows and columns.",
        msg, P_2, P_2
      )
      msg <- ifelse(Matrix::rankMatrix(sigma) == P_2, TRUE, "Input error for Sigma^-1: Sigma^-1 is not positive-definite")
      if (is.character(msg)) stop(msg)
    },
    validate_x = function(Xf, Xr) {
      msg_xf <- checkmate::check_matrix(Xf,
        mode = "numeric",
        any.missing = F,
        min.rows = 2, min.cols = 1
      )
      if (is.logical(msg_xf)) {
        P_1 <- ncol(Xf)
        msg_xf <- ifelse(
          Matrix::rankMatrix(Xf) == P_1, TRUE, "Xf is not full-rank"
        )
      }
      private$.print_error("Input error for Xf: %s.", msg_xf)
      if (!is.null(Xr)) {
        msg_xr <- checkmate::check_matrix(Xr,
          mode = "numeric",
          any.missing = F,
          nrows = nrow(Xf), min.cols = 1
        )
        private$.print_error("Input error for Xr: %s.", msg_xr)
      }
    },
    validate_inputs = function(Xf, Xr, Sigma_inv, Phi_inv) {
      self$validate_x(Xf, Xr)
      if (!is.null(Xr)) {
        P_2 <- ncol(Xr)
        self$validate_sigma(Sigma_inv, P_2)
      }
      self$validate_phi(Phi_inv, nrow(Xf))
    },
    validate_w = function(w) {
      # TODO
    }
  )
)
