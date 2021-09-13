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
