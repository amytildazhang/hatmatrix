
# test output of link functions, e.g.
# linkFun <- .get_link_fun("cloglog"); linkFun(mu, deriv = 1); linkFun(mu)

test_that("Output of link functions work", {
  expect_equal(
      phi(qnorm(0.3), family = stats::binomial(link = "probit"), inverse = FALSE),
      0.3*(1 - 0.3)/dnorm(qnorm(0.3))^2
  )

    bifam <- stats::binomial()
    expect_equal(
        phi(bifam$linkfun(0.3), family = stats::binomial(link = "logit"), inverse = FALSE),
        1/(0.3*(1 - 0.3))
    )

    # test scaling
    expect_equal(
        phi(qnorm(0.3), family = stats::binomial(link = "probit"), inverse = FALSE, scale = 10),
        0.3*(1 - 0.3)/(10*dnorm(qnorm(0.3))^2)
    )

})
