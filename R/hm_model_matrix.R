#' #' R6 class for XfrMatrix object
#' #'
#' #' @description
#' #' XfrMatrix [R6] object created by [xfr_matrix].
#' #'
#' #' @seealso xfr_matrix
#' #'
#' XfrMatrix <- R6::R6Class(
#'   "XfrMatrix",
#'   private = list(),
#'   public = list(
#'     Xf = NULL,
#'     Xr = NULL,
#'     initialize = function(Xf, Xr = NULL) {
#'       checkmate::check_matrix(Xf,
#'                               mode = "numeric",
#'                               any.missing = F, min.rows = 1,
#'                               min.cols = 1
#'       )
#'       if (!is.null(Xr)) {
#'         checkmate::check_matrix(Xr,
#'                                 mode = "numeric",
#'                                 any.missing = F, min.rows = 1,
#'                                 min.cols = 1
#'         )
#'       }
#'
#'       self$Xf <- Xf
#'       self$Xr <- Xr
#'     }
#'   )
#' )
#'
#' #' Create a XfrMatrix object
#' #'
#' #' Create a [XfrMatrix] object.
#' #'
#' #' @template xf-arg
#' #' @template xr-arg
#' #'
#' #' @return
#' #'
#' #' @examples
#' xfr_matrix <- function(Xf, Xr = NULL) {
#'   XfrMatrix$new(Xf, Xr)
#' }
#'
#'
#'
#'
#'
#'
#'
#' hmm_from <- function(obj) {
#'   UseMethod("hmm_from")
#' }
#'
#'
#'
#' #' Title
#' #'
#' #' Internal function
#' #'
#' #' @param formula
#' #' @param ftype
#' #' @param ...
#' #' @param family
#' #' @param env
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' @importFrom stats gaussian
#' hmm_from.formula <- function(formula,
#'                              family = gaussian(),
#'                              ftype = c(
#'                                "lm", "glm", "lmer", "glmer",
#'                                "gamm4"
#'                              ),
#'                              env = parent.frame(),
#'                              ...) {
#'   mf <- match.call()
#'   mf$env <- mf$ftype <- NULL
#'
#'   # special use cases for formula changes here
#'   if (ftype == "gamm4" & "random" %in% names(mf)) {
#'     fake.formula <- as.character(mgcv::interpret.gam(formula)$fake.formula)
#'     mf$formula <- paste(fake.formula[2], fake.formula[1], fake.formula[3],
#'       "+", mf$random[2],
#'       collapse = " "
#'     )
#'   }
#'
#'   has_random <- !(ftype %in% c("lm", "glm"))
#'
#'   # call appropriate model frame-building function
#'   mf[[1L]] <- quote(
#'     ifelse(has_random, lme4::glFormula, model.matrix)
#'   )
#'
#'   mm <- eval(mf, env)
#'
#'   Xf <- ifelse(has_random, mm$X, mm)
#'   Xr <- ifelse(has_random, t(as.matrix(mm$reTrms)), mm)
#'
#'   hm_model_matrix(Xf, Xr)
#' }
#'
#' hmm_from.lm <- function(lmobj) {
#'   # validate is lmobj
#'   # lmobj <- lm(cyl ~ mpg, data = mtcars, subset = 1:10)
#'
#'   if ("model" %in% names(lmobj)) {
#'     return(lmobj$model)
#'   }
#'
#'   # build call to
#'   mf <- lmobj$call
#'   m <- match(c("formula", "data", "subset", "na.action", "drop.unused.levels", "xlev", "contrasts"),
#'     table = names(mf), nomatch = 0L
#'   )
#'   mf <- mf[c(1L, m)]
#'   mf[[1L]] <- quote(hmm_from.formula)
#'   mf$ftype <- "lm"
#'
#'   env <- attr(lmobj$terms, ".Environment")
#'   xcall <- match.call(hmm_from.formula, mf, envir = env)
#'   eval(xcall)
#' }
#'
#'
#' hmm_from.merMod <- function(fit) {
#'   call <- fit@call
#'   data <- fit@frame
#' }
#'
#'
#' hmm_from.stanreg <- function(fit) {
#'   # formula <- fit$formula
#'   # call <- fit$call
#'   #
#'   # m <- match(formalArgs(lme4::glFormula),
#'   #            table = names(call), nomatch = 0L)
#'   # mc <- call[c(1L,build_x)]
#'   # mc$data <- fit$data
#'   # eval(match.call(build_x, mc))
#' }
