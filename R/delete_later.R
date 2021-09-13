#'
#' # UNUSED
#'
#' #' Conditional to throw an error or return FALSE
#' #'
#' #' Convenience function to avoid repeating the same structure.
#' #'
#' #' @param errormsg If throwing an error, the message to display.
#' #' @param throw_error Whether to throw an erorr or display the message.
#' #'
#' #' @return If `throw_error = FALSE`, returns `FALSE`.
#' #' @export
#' #'
#' #' @examples
#' .conditional_error <- function(errormsg, throw_error = TRUE) {
#'     if (throw_error) {
#'         stop(errormsg)
#'     } else {
#'         return(FALSE)
#'     }
#' }
#'
#'
#' #' Check if column exists in dataframe
#' #'
#' #' @param df Dataframe
#' #' @param col Column index (as a number) or column name.
#' #' @param throw_error Whether to throw an error or return FALSE
#' #'
#' #' @return Boolean value based on check.
#' #' @export
#' #'
#' #' @examples
#' #' @seealso .conditional_error
#' .check_col_in_df <- function(df, col, throw_error = TRUE) {
#'     if (is.numeric(col) & ncol(df) < col) {
#'         .conditional_error(
#'             sprintf("Dataframe has only %s columns and
#'                                    the %sth column was requested",
#'                     ncol(df), col)
#'         )
#'     } else if (is.character(col) & !(col %in% colnames(df))) {
#'         .conditional_error(
#'             sprintf("There is no column named %s in dataframe", col)
#'         )
#'     }
#'     invisible(TRUE)
#' }
#'
#'
#' HMData <- R6::R6Class(
#'     "HMData",
#'     private = list(
#'         .data = NULL,
#'         # metadata, derived in private$update
#'         # @field .d Stores calculated number of parameters
#'         .d = 0,
#'
#'         # @field .y Stores number of rows in `data`
#'         .n = 0,
#'         # @description Function that includes all checks to validate
#'         #    the `data`
#'         .check_data = function(data) {
#'             checkmate::assert_data_frame(data, any.missing = F)
#'         },
#'
#'         # @description Convenience function that calls .check_data,
#'         #    adds the id__ column to the data if not specified, and then
#'         #    updates private$.data
#'         .update_data = function(data, id_col = NULL) {
#'             private$.check_data(data)
#'
#'             if (is.null(id_col)) {
#'                 data <- dplyr::mutate(data, id__ = 1:nrow(data))
#'             } else {
#'                 data$id__ <- data[[id_col]]
#'                 data[[id_col]] <- NULL
#'             }
#'
#'             private$.data <- data
#'             private$.n <- nrow(data)
#'
#'             invisible(self)
#'         },
#'
#'         .x_function = NULL,
#'         .mask = NULL,
#'
#'         # calculated via private$update(x_function, mask)
#'         .X = NULL,
#'         ### objects which should only be changed internally
#'         # @field .y Stores supplied vector of responses for each row in
#'         #   `data`
#'         .y = NULL,
#'
#'         # @description Function that includes all checks to validate
#'         #    the numeric vector y (or yhat)
#'         .check_y = function(y, data) {
#'             if (length(y) == 1) {
#'                 .check_col_in_df(df, y, throw_error = TRUE)
#'                 y <- df[, y]
#'             }
#'             checkmate::assert_atomic_vector(y, len = nrow(data))
#'
#'         },
#'
#'         # @description Convenience function that calls .check_y and then
#'         #    updates private$.y
#'         .update_y = function(y, data) {
#'             if (length(y) == 1) {
#'
#'             }
#'             private$.check_y(y, yhat)
#'             if (yhat) {
#'                 private$.yhat <- y
#'             } else {
#'                 private$.y <- y
#'             }
#'         }
#'
#'     ),
#'     public = list(
#'         print = function() {
#'             #TODO
#'         },
#'         initialize = function(data, y = NULL, id_col = NULL) {
#'
#'             ## store initialization values
#'             # store immutable objects
#'             if (!is.null(y)) {
#'                 private$.update_y(ycol, data)
#'             }
#'             private$.update_data(data, id_col = id_col)
#'
#'
#'         }
#'     )
#' )
#'
#'
#'
#' #' Conditional to throw an error or return FALSE
#' #'
#' #' Convenience function to avoid repeating the same structure.
#' #'
#' #' @param errormsg If throwing an error, the message to display.
#' #' @param throw_error Whether to throw an erorr or display the message.
#' #'
#' #' @return If `throw_error = FALSE`, returns `FALSE`.
#' #' @export
#' #'
#' #' @examples
#' .conditional_error <- function(errormsg, throw_error = TRUE) {
#'     if (throw_error) {
#'         stop(errormsg)
#'     } else {
#'         return(FALSE)
#'     }
#' }
#'
#'
#' #' Check if column exists in dataframe
#' #'
#' #' @param df Dataframe
#' #' @param col Column index (as a number) or column name.
#' #' @param throw_error Whether to throw an error or return FALSE
#' #'
#' #' @return Boolean value based on check.
#' #' @export
#' #'
#' #' @examples
#' #' @seealso .conditional_error
#' .check_col_in_df <- function(df, col, throw_error = TRUE) {
#'     if (is.numeric(col) & ncol(df) < col) {
#'         .conditional_error(
#'             sprintf("Dataframe has only %s columns and
#'                                    the %sth column was requested",
#'                     ncol(df), col)
#'         )
#'     } else if (is.character(col) & !(col %in% colnames(df))) {
#'         .conditional_error(
#'             sprintf("There is no column named %s in dataframe", col)
#'         )
#'     }
#'     invisible(TRUE)
#' }
#'
#'
#' validate_mask <- function(mask) {
#'     checkmate::assert_logical(mask, any.missing = FALSE,
#'                                    min.len = 2)
#' }
#'
#'
#' HMModelMatrix <- R6::R6Class(
#'     "HMModelMatrix",
#'     private = list(
#'         .fixef_idx = c(),
#'         .X = matrix(),
#'         .x_function = function(mask) {
#'             .masked_X(private$.X, mask)
#'         },
#'         .masked_X = function(X, mask) {
#'             masked_X <- matrix(0, nrow = nrow(X) + sum(!mask), ncol = ncol(X))
#'             masked_X[mask, ] <- X
#'             masked_X
#'         },
#'         .check_x_function = function(x_function) {
#'             checkmate::assert_function(x_function,
#'                                       args = c("mask"),
#'                                       ordered = T)
#'         },
#'         .update_x_function = function(x_function) {
#'             private$.check_xfunction(x_function)
#'
#'             private$.x_function <- x_function
#'         },
#'
#'         .check_X = function(X, n) {
#'             check_out <- checkmate::check_matrix(
#'                 X,
#'                 mode = "numeric",
#'                 any.missing = FALSE,
#'                 min.cols = 1, nrows = n)
#'             # TODO: check that all columns have variable names
#'             if (is.character(check_out)) {
#'                 stop(sprintf("Supplied x_function %s",
#'                              stringr::str_to_lower(check_out)))
#'             }
#'
#'         }
#'     ),
#'     public = list(
#'         print = function(){
#'             #TODO
#'         },
#'         xfun_args = list(),
#'         initialize = function(Xf = NULL, Xr = NULL,
#'                               x_function = NULL,  ...) {
#'             #TODO: check both XF and XR supplied OR x_function supplied
#'             if (is.null(x_function)) {
#'                 X <- cbind(Xf, Xr)
#'                 private$.fixef_idx <- 1:ncol(Xf)
#'                 private$.check_X(X)
#'             } else {
#'                 private$.check_xfunction(x_function)
#'
#'             }
#'             # can i save ... as arguments to pass to x_function?
#'             self$xfun_args <- as.list(...) #TODO
#'         },
#'         X = function(mask = rep(T, private$.n)) {
#'             validate_mask(mask)
#'             Xmats <- do.call(private$.x_function,
#'                          c(list(mask = mask), xfun_args))
#'             X <- private$.masked_X(X, mask)
#'             private$.check_X(X, length(mask))
#'             X
#'         }
#'     ),
#'     active = list(
#'         #' @field
#'         #' Return the fixed effect column indices in X
#'         fixef_idx =  pryr::unenclose(read_only(.fixef_idx))
#'     )
#' )
#'
#' HMRprec <- R6::R6Class(
#'     "HMRprec",
#'     private = list(
#'         .precision = matrix(),
#'         .check_prec = function(X, n) {
#'             # TODO
#'         }
#'     ),
#'     public = list(
#'         print = function(){
#'             #TODO
#'         },
#'         initialize = function(prec = NULL, sigmalist = NULL) {
#'         },
#'         Precision = function(fixef_idx = NULL,  ) {
#'             if (checkmate::check_matrix(rprec, "numeric", any.missing = F)) {
#'                 iSig <- self$rprec
#'             } else {
#'                 sigi_diag <- rep(0, ncol(X))
#'                 for (partition in names(rprec)) {
#'                     sigi_diag[stringr::str_detect(colnames(X), partition)] <- 1/rprec[[partition]]^2
#'
#'                 }
#'
#'                 iSig <- diag(sigi_diag)
#'             }
#'         }
#'
#'     ),
#'     active = list(
#'         #' @field
#'         #' Return the fixed effect column indices in X
#'         fixef_idx =  pryr::unenclose(read_only(.fixef_idx))
#'     )
#' )
#'
#' #' HatMatrix R6 object
#' #'
#' #' @description
#' #'     Along with the fields and methods listed here, `HatMatrix` inherits from BasicHatMatrix
#' #'     and thus includes (or overwrites, if the names are the same) all methods
#' #'     and fields in BasicHatMatrix.
#' #'
#' #'
#' #' @details
#' #' here
#'
#' HatMatrix <- R6::R6Class(
#'     "HatMatrix",
#'     # TODO: make cloneable -- to use same data but different models/masks
#'     cloneable = TRUE,
#'
#'     private = list(
#'            # @field .data Stores supplied dataframe
#'         .data = NULL,
#'
#'         # @field .W Stores supplied matrix of weight W for rows in `data`
#'         .W = NULL,
#'
#'         #@ field .W_rowid For each row in W, stores its corresponding index
#'         #    in `data` such that `data[.W_rowid, ]` produces rows of `data`
#'         #    that match W
#'         .W_rowid = NULL,
#'
#'         #@ field .W_rowid For each column in W, stores its corresponding index
#'         #    in `data` such that `data[.W_colid, ]` produces rows of `data`
#'         #    that match W
#'         .W_colid = NULL,
#'
#'
#'         # @field .yhat Stores supplied vector of data-level estimates for
#'         #    each row in `data`
#'         .yhat = NULL,
#'
#'         # @field .clusters Stores supplied character vector of cluster labels
#'         #    for each row in `data`
#'         .clusters = NULL,
#'
#'         # @field .contains_new Whether the weight matrix W contains weights for new data; default is F
#'         .contains_new = F,
#'
#'         # @description Centralized method to check and produce error if supplied object
#'         #     does not have length or rows matching the number of datapoints.
#'         # @param obj The object to check
#'         # @param name The name of the object, as a string
#'         # @param type The type of object, as a string. Available options are
#'         #     "vector" and "matrix".
#'         .assert_length_matches_rows = function(obj, name, type = c("vector", "matrix")) {
#'             if (type == "vector") {
#'                 condition <- length(obj) != self$n
#'                 description <- "Length of vector"
#'             } else {
#'                 condition <- nrow(obj) != self$n
#'                 description <- "Number of rows in matrix"
#'             }
#'
#'             if (condition) {
#'                 stop(sprintf("%s `%s` must match number of rows in `data`",
#'                              description, name))
#'             }
#'         },
#'
#'         # functions checking and related to W
#'
#'         # @description Interior method allowing for higher-level functions
#'         #     to remain the same across all Slice objects.
#'         # @return The \eqn{N \times N} matrix of weights `W` relevant
#'         #     to the Slice object (e.g. for BasicSlice, the supplied W; for
#'         #     Slice, the calculated W; for CVSlice, the cv_W calculated over
#'         #     cross-validation folds.)
#'         .get_W = function() {
#'             private$.W
#'         },
#'
#'
#'         # @description Convenience function that checks whether the rows
#'         #    of a matrix W sum to 1 and throws an error if they do not.
#'         # @param W A numerc matrix
#'         .check_rowsums_W = function(W) {
#'             sums <- rowSums(W)
#'             dif <- round(sums - 1, digits = 6)
#'             if (any(dif != 0)) {
#'                 stop("Rows of W do not sum to 1. This is possible if W contains weights for new data; if so, call Slice$new again and set contains_new to TRUE.")
#'             }
#'         },
#'
#'
#'         # @description Function that includes all checks to validate
#'         #    the weight matrix W
#'         .check_W = function(W) {
#'
#'             if (!private$.contains_new) {
#'                 private$.assert_length_matches_rows(W, "W", "matrix")
#'                 private$.check_rowsums_W(W)
#'             }
#'
#'             checkmate::assert_matrix(W, mode = "numeric",
#'                                      any.missing = FALSE,
#'                                      ncols = length(private$.W_colid),
#'                                      nrows = length(private$.W_rowid))
#'
#'             # check symmetry
#'             # checkmate::assert_true(W == t(W))
#'
#'             # check the rows (and thus columns) sum to 1
#'
#'             # check positive along diagonal
#'             # checkmate::assert_true(all(diag(W) > 0))
#'
#'
#'
#'         },
#'
#'         # @description Convenience function that calls .check_W and then
#'         #    updates private$.W
#'         .update_W = function(W) {
#'
#'             private$.check_W(W)
#'             private$.W <- W
#'
#'         },
#'
#'
#'         # @description Function that includes all checks to validate
#'         #    the numeric vector y (or yhat)
#'         .check_y = function(y, yhat = F) {
#'             if (private$.contains_new) {
#'                 checkmate::assert_atomic_vector(y)
#'             } else {
#'                 checkmate::assert_atomic_vector(y, any.missing = FALSE)
#'             }
#'
#'             # for more readable output
#'             private$.assert_length_matches_rows(y, ifelse(yhat, "yhat", "y"), "vector")
#'
#'         },
#'
#'         # @description Convenience function that calls .check_y and then
#'         #    updates private$.y
#'         .update_y = function(y, yhat = F) {
#'             private$.check_y(y, yhat)
#'             if (yhat) {
#'                 private$.yhat <- y
#'             } else {
#'                 private$.y <- y
#'             }
#'         },
#'
#'
#'
#'
#'         # @description Function that includes all checks to validate
#'         #    supplied logical vector of training mask
#'         .check_mask = function(mask) {
#'             checkmate::assert_logical(mask,
#'                                       any.missing = FALSE,
#'                                       all.missing = FALSE)
#'             private$.assert_length_matches_rows(mask, "mask", "vector")
#'
#'             if (sum(mask) == 0) {
#'                 stop("Supplied vector for `mask` must include at least `TRUE`` value.")
#'             }
#'
#'         },
#'
#'
#'         .rprec = NULL,
#'         .Phi_inv = NULL,
#'
#'
#'         .V = NULL,
#'
#'
#'         .vars = list(),
#'         #
#'
#'         .check_x_function = function(x_function) {
#'             checkmate::assert(
#'                 checkmate::check_matrix(x_function, "numeric", any.missing = F),
#'                 checkmate::check_function(x_function,
#'                                           args = c("data", "mask"),
#'                                           ordered = T)
#'             )
#'
#'         },
#'         .update_x_function = function(x_function) {
#'             private$.check_xfunction(x_function)
#'
#'             if (checkmate::test_matrix(x_function)) {
#'                 private$.X <- x_function
#'                 x_function <- function(data, train_mask) {
#'                     private$.X
#'                 }
#'             }
#'             private$.x_function <- x_function
#'         },
#'
#'         .update_mask = function(mask) {
#'             private$.check_mask(mask)
#'             private$.mask <- mask
#'         },
#'
#'
#'         .check_X = function(X) {
#'             check_out <- checkmate::check_matrix(
#'                 X,
#'                 mode = "numeric",
#'                 any.missing = FALSE,
#'                 min.cols = 1, nrows = self$n)
#'             # TODO: check that all columns have variable names
#'             if (is.character(check_out)) {
#'                 stop(sprintf("Supplied x_function %s",
#'                              stringr::str_to_lower(check_out)))
#'             }
#'
#'         },
#'
#'         # individual update calculations and checks for X, V, W
#'         .update_X = function() {
#'             X <- self$xfunction(self$data, self$mask)
#'
#'             # check that x has correct number of rows and is named
#'
#'             private$.X <- X
#'         },
#'
#'         .train_X = function(train_mask = self$mask) { # zero-out datapoints not in the mask
#'             X <- matrix(0, nrow = nrow(private$.X), ncol = ncol(private$.X))
#'             X[train_mask,] <- private$.X[train_mask, ]
#'             colnames(X) <- colnames(private$.X)
#'             X
#'         },
#'
#'         .check_Phiinv = function() {
#'
#'         },
#'
#'         .train_Phiinv = function(train_mask = self$mask) {
#'             private$.check_Phiinv()
#'             private$.Phi_inv[train_mask]
#'         }
#'
#'         .update_V = function() {
#'             X <- private$.train_X()
#'             sig_i <- self$Sigma_inv(X, private$.rprec)
#'
#'             stopifnot(all(sig_i >= 0), any(sig_i > 0))
#'             private$.V <-  solve(t(X) %*% diag(private$.train_Phiinv()) %*% X + sig_i)
#'
#'         },
#'
#'         .update_W = function() {
#'             masked_X <- private$.train_X()
#'
#'             W <- pc_W(masked_X, Phi_inv = private$.train_Phiinv(), V = self$V)
#'             if (sum(!self$mask) > 0) {
#'                 # no weight on test data points
#'                 stopifnot(all(W[, !self$mask, drop = FALSE] == 0))
#'             }
#'
#'             # weights sum to 1
#'             stopifnot(all(round(rowSums(private$.X %*% W) - 1, digits = 6) == 0))
#'
#'             private$.W <- private$.X %*% W
#'
#'
#'
#'         },
#'
#'
#'         .update_vars = function() {
#'             # TODO
#'
#'             return(NULL)
#'
#'         },
#'
#'         # @description
#'         #
#'         # @details Intentionally over-built function to allow for changing
#'         #    the mask so it can be easily extended to CVHatMatrix.
#'         .update = function(x_function = NULL, mask = NULL) {
#'             # TODO: is there a common way to check supplied variables?
#'             stopifnot(!is.null(x_function) | !is.null(mask))
#'
#'             if (!is.null(x_function)) {
#'                 private$.update_xfunction(x_function)
#'             }
#'
#'             if (!is.null(mask)) {
#'                 private$.update_mask(mask)
#'             }
#'
#'             private$.update_X()
#'             private$.update_Phi()
#'             private$.update_V()
#'             private$.update_W()
#'
#'             private$.d <- ncol(private$.X)
#'
#'             invisible(self)
#'         }
#'     ),
#'     public = list(
#'
#'         clusters = NULL,
#'
#'         print = function() {
#'
#'             # TODO: brief summary of object and what can be done/how to access fields
#'         },
#'         #' @description Create a new HatMatrix object by supplying
#'         #'
#'         #' @details
#'         #' @param x_function Either the model matrix `X` of the form used to
#'         #'     calculate the conditional covariance (as described in [[CITE]]) or
#'         #'     a function which takes as parameters `data` and `mask` and produces
#'         #'     the model matrix `X`.
#'         #' @param data A dataframe with \code{N} rows containing all covariates
#'         #'     used in the model. It may include additional columns if [[CONTINUE]]
#'         #' @param rprec Either the random effect precision matrix or a list of
#'         #'     the precision values (if random effects are i.i.d). The name of
#'         #'     each list entry must match the name(s) of the
#'         #'     column(s) in `X` for which the precision applies.
#'         #' @param id_col [[CP from HatMatrix_basic]]
#'         #' @param Phi_inv
#'         #' @param glm
#'         #' @param y
#'         #' @param mask
#'         #'
#'         #' @return A HatMatrix object.
#'         initialize =  function(x_function,  rprec,  Phi_inv = NULL,
#'                                data = NULL, mask = NULL, id_col = NULL) {
#'
#'             # TODO: allow for family name and link function
#'             # i.e. supply glm = list(family = ###, link = ###)
#'
#'
#'             # check basic data formats
#'             stopifnot(all(names(rprec) %in% colnames(data)))
#'             if (is.null(mask)) {
#'                 mask <- rep(TRUE, nrow(data))
#'             }
#'
#'
#'
#'             ## store initialization values
#'             # store immutable objects
#'             private$.y <- y
#'
#'             private$.update_data(data, id_col = id_col)
#'
#'
#'             private$.rprec <- rprec
#'
#'             private$.Phi_inv <- Phi_inv
#'
#'
#'             # calculate X, V, W, d, and n
#'
#'
#'             private$.update(x_function, mask)
#'
#'
#'         },
#'
#'
#'         print = function() {
#'             # TODO:
#'             # - model formula?
#'             # - number of datapoints
#'             # - size of mask
#'             # - clusters?
#'
#'         },
#'
#'
#'         #' @description Split the `HatMatrix` object into multiple sub-models based
#'         #'     on blocks of the model design matrix `X`, as described in [[CITE HERE]]
#'         #' @param subsets A character vector of sub-models to create. Entries
#'         #'     be a subset of all possible output from `.fun`.
#'         #' @param .fun A function which takes as input the column names of
#'         #'    `obj$X` and produces as output an equally-sized character vector,
#'         #'    with entries denoting the block each column belongs to.
#'         #' @examples
#'         split_W = function(subsets, .fun) {
#'
#'             # returns a list of HatMatrix objects, named by subsets,
#'             # and itself, named by "full"
#'             newx_ <- function(set, .fun, xfun) {
#'                 function(data, mask) {
#'                     X <- xfun(data, mask)
#'                     col_sets <- .fun(colnames(X))
#'                     X[, col_sets == lazyeval::lazy_eval(set), drop = FALSE]
#'                 }
#'             }
#'
#'
#'
#'             sig_i <- self$Sigma_inv(private$.X, private$.rprec)
#'             sig_subsets <- .fun(colnames(private$.X))
#'             contains_random <- purrr::map_lgl(subsets, function(set) {
#'                 any(sig_i[sig_subsets == set] > 0)
#'             })
#'
#'             if (!all(contains_random)) {
#'                 stop(sprintf("Subsets created by `.fun` must contain random effects. Subsets that do not are: '%s'",
#'                              paste(subsets[!contains_random], collapse = "' '")))
#'             }
#'
#'
#'             c(list("full" = self),
#'               lapply(subsets, function(set) {
#'                   # adjust x_function to save only a subset of columns
#'
#'                   # create a new HatMatrix object with same properties otherwise
#'                   HatMatrix$new(
#'                       pryr::unenclose(newx_(lazyeval::lazy(set), .fun, self$x_function)),
#'                       data = self$data, rprec = self$rprec,
#'                       y = self$y, Phi_inv = self$Phi_inv,
#'                       glm = self$glm, mask = self$mask
#'                   )
#'               }) %>% purrr::set_names(subsets)
#'             )
#'
#'
#'
#'
#'
#'
#'         },
#'
#'
#'         #' @description
#'         #'
#'         Sigma_inv = function(X = self$x, rprec = self$rprec) {
#'             if (checkmate::check_matrix(rprec, "numeric", any.missing = F)) {
#'                 self$rprec
#'             } else {
#'                 sigi_diag <- rep(0, ncol(X))
#'                 for (partition in names(rprec)) {
#'                     sigi_diag[stringr::str_detect(colnames(X), partition)] <- 1/rprec[[partition]]^2
#'
#'                 }
#'
#'                 diag(sigi_diag)
#'             }
#'         }
#'
#'
#'
#'     ),
#'
#'     active = list(
#'         #' @field y The supplied (or derived if using `to_slice.x`)
#'         #'     vector of responses
#'         y = pryr::unenclose(read_only(.y)),
#'
#'         #' @field yhat The supplied (or derived if using `to_slice.x`)
#'         #'     vector of fitted values
#'         yhat = pryr::unenclose(read_only(.yhat)),
#'
#'         #' @field data The supplied (or derived if using `to_slice.x`)
#'         #'     data
#'         data = pryr::unenclose(read_only(.data)),
#'
#'         #' @field W_rowid The supplied (or derived if using `to_slice.x`)
#'         #'     data
#'         W_rowid = pryr::unenclose(read_only(.W_rowid)),
#'
#'         #' @field W_rowid The supplied (or derived if using `to_slice.x`)
#'         #'     data
#'         W_colid = pryr::unenclose(read_only(.W_colid)),
#'
#'         #' @field d The total number of coefficient parameters
#'         d = pryr::unenclose(read_only(.d)),
#'
#'         #' @field n The total number of rows in dataset
#'         n = pryr::unenclose(read_only(.n)),
#'
#'         #' #' @field clusters The supplied vector of cluster labels.
#'         #' clusters = pryr::unenclose(read_only(.clusters)),
#'
#'         #' @field W The \eqn{N \times N} (or \eqn{N \times C} if `clusters` is
#'         #'     suppplied) matrix of factors
#'         W = function(val) {
#'             if (missing(val)) {
#'                 W <- private$.get_W()
#'
#'                 if (length(self$clusters) > 0) {
#'                     W <- cluster_W(W, self$clusters)
#'                 }
#'                 dimnames(W) <- NULL
#'                 W
#'
#'             } else {
#'                 message("Value is read-only")
#'             }
#'         },
#'
#'         #' @field
#'         #' Return the supplied (or derived if using `to_HatMatrix.x`) posterior means
#'         #' of coefficient standard deviations
#'         rprec = pryr::unenclose(read_only(.rprec)),
#'
#'         #' @field
#'         #' Return the supplied (or derived if using `to_s lice.x`) vector of
#'         #' data noise variances
#'         Phi_inv = pryr::unenclose(read_only(.Phi_inv)),
#'
#'           #' @field
#'         #' Return the derived design matrix as a result of calling
#'         #' \code{xfunction(data, mask)}
#'         X =  pryr::unenclose(read_only(.X)),
#'
#'         #' @field
#'         #' Return the derived conditional variance
#'         V =  pryr::unenclose(read_only(.V)),
#'
#'         #' @field
#'         #' Return the supplied mask
#'         mask = pryr::unenclose(read_only(.mask)),
#'
#'
#'         x_function = pryr::unenclose(read_only(.x_function)),
#'
#'         allow_cv = pryr::unenclose(read_only(.allow_cv))
#'     )
#'
#' )
#'
#'
#'
#'
