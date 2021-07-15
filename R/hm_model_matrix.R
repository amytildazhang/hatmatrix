
new_hm_model_matrix <- function(Xf = matrix(), Xr = NULL) {
    checkmate::check_matrix(Xf, mode = "numeric",
                            any.missing = F, min.rows = 1,
                            min.cols = 1)
    if (!is.null(Xr)) {
        checkmate::check_matrix(Xr, mode = "numeric",
                                any.missing = F, min.rows = 1,
                                min.cols = 1)

    }
    structure(
        cbind(Xf, Xr),
        class = "hm_model_matrix",
        fixef_idx = 1:ncol(Xf)
    )
}

hmm_to_xf <- function(hmm) {
    fixef_idx <- attr(hmm, "fixef_idx")
    unclass(hmm)[, fixef_idx]
}

hmm_to_xr <- function(hmm) {
    fixef_idx <- attr(hmm, "fixef_idx")
    unclass(hmm)[, -fixef_idx]
}

validate_hm_model_matrix <- function(hm_mm) {
    #TODO check that class is hm_model_matrix
    nfixeff <- max(attr(hm_mm, "fixef_idx"))
    checkmate::assert_matrix(hm_mm, mode = "numeric",
                             any.missing = F,
                             min.rows = T,
                             min.cols = nfixeff)
    Xf <- unclass(hm_mm)
    Xf <- Xf[, 1:nfixeff]
    stopifnot(Matrix::rankMatrix(Xf) == nfixeff)


}


error_hm_model_matrix <- function(error) {
    #TODO
    stop()
}



hm_model_matrix <- function(Xf = NULL, Xr = NULL) {
    hm_mm <- new_hm_model_matrix(Xf, Xr)

    tryCatch(validate_hm_model_matrix(hm_mm),
             error = error_hm_model_matrix)

}





hmm_from <- function(obj) {
    UseMethod("hmm_from")
}



#' Title
#'
#' Internal function
#'
#' @param formula
#' @param ftype
#' @param ...
#' @param family
#' @param env
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom stats gaussian
hmm_from.formula <- function(formula,
                           family = gaussian(),
                           ftype = c("lm", "glm", "lmer", "glmer",
                                     "gamm4"),
                           env = parent.frame(),
                           ...) {
    mf <- match.call()
    mf$env <- mf$ftype <- NULL

    # special use cases for formula changes here
    if (ftype == "gamm4" & "random" %in% names(mf)) {
        fake.formula <- as.character(mgcv::interpret.gam(formula)$fake.formula)
        mf$formula <- paste(fake.formula[2], fake.formula[1], fake.formula[3],
                            "+", mf$random[2], collapse = " ")
    }

    has_random <- !(ftype %in% c("lm", "glm"))

    # call appropriate model frame-building function
    mf[[1L]] <- quote(
        ifelse(has_random, lme4::glFormula, model.matrix)
    )

    mm <- eval(mf, env)

    Xf <- ifelse(has_random, mm$X, mm)
    Xr <-  ifelse(has_random, t(as.matrix(mm$reTrms)), mm)

    hm_model_matrix(Xf, Xr)
}

hmm_from.lm <- function(lmobj) {
    # validate is lmobj
    # lmobj <- lm(cyl ~ mpg, data = mtcars, subset = 1:10)

    if ("model" %in% names(lmobj)) {
        return(lmobj$model)
    }

    # build call to
    mf <- lmobj$call
    m <- match(c("formula", "data", "subset", "na.action", "drop.unused.levels", "xlev", "contrasts"),
               table = names(mf), nomatch = 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(hmm_from.formula)
    mf$ftype <- "lm"

    env <- attr(lmobj$terms, ".Environment")
    xcall <- match.call(hmm_from.formula, mf, envir = env)
    eval(xcall)
}


hmm_from.merMod <- function(fit) {
    call <- fit@call
    data <- fit@frame

}


hmm_from.stanreg <- function(fit) {
    # formula <- fit$formula
    # call <- fit$call
    #
    # m <- match(formalArgs(lme4::glFormula),
    #            table = names(call), nomatch = 0L)
    # mc <- call[c(1L,build_x)]
    # mc$data <- fit$data
    # eval(match.call(build_x, mc))
}
