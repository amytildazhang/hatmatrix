#' Function to calculate noise variance
#'
#' @param Xbeta
#' @param family
#' @param link
#'
#' @return Vector 
#' @export
#'
#' @examples
phi <- function(Xbeta, family, inverse = TRUE) {
    # validate family object has class family or extended.family
    # validate Xbeta

    mu <- family$linkinv(Xbeta)
    var <- family$variance(mu)

    linkFun <- .get_link_fun(family$link)
    deltaMultiplier <- linkFun(mu, deriv = 1) # scaling due to Delta method

    if (inverse) {
        exp(log(deltaMultiplier - var))
    } else {
        exp(log(var - deltaMultiplier))
    }

}



.get_link_fun <- function(link) {
    linkFun <- switch(link,
           "log" = VGAM::loglink,
           "cloglog" = VGAM::clogloglink,
           "logit" = VGAM::logitlink,
           "probit" = VGAM::probitlink,
           "cauchit" = VGAM::cauchitlink,
           "identity" = VGAM::identitylink,
           "inverse" = VGAM::reciprocallink,
           "sqrt" = pryr::unenclose(.set_powerlink(0.5))
    )
    if (is.null(linkFun)) {
        if (stringr::str_detect(link, "1/mu")) {
            pow <- stringr::str_replace(stringr::str_extract(link, "\\^d+"), "\\^", "")
            linkFun <- pryr::unenclose(.set_powerlink(as.numeric(pow)))
        }
    }

    linkFun
}

.set_powerlink <- function(power) {
    function(theta, ....) VGAM::powerlink(theta, power, ...)
}




phi_from <- function(x) {
    UseMethod("phi_from")
}


phi_from.stanreg <- function(x) {
    #TODO
}

