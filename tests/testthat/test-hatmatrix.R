
Xf <- model.matrix(mpg ~ disp, data = mtcars)
Xr <- model.matrix(mpg~-1 + factor(cyl), data = mtcars)
Sigma_inv <- diag(rep(3, ncol(Xr)))
Phi_inv <- rep(2, nrow(Xf))

P <- ncol(Xf) + ncol(Xr)
check_sig <- matrix(0, P, P)
check_sig[(1 + ncol(Xf)):P, (1 + ncol(Xf)):P] <- Sigma_inv


thm <- HatMatrixCalculator$new(Xf, Phi_inv, Xr, Sigma_inv)
X <- cbind(Xf, Xr)
W <- X %*% solve(t(X) %*% diag(Phi_inv) %*% X + check_sig) %*% t(X) %*% diag(Phi_inv)


test_that("Internal consistency of HatMatrix", {
    expect_snapshot(thm)

    expect_equal(thm$Sigma_inv, Sigma_inv)
    expect_equal(thm$expanded_sigma(),check_sig)
    expect_equal(thm$Xf, Xf)
    expect_equal(thm$Xr, Xr)

    expect_null(thm$i)
    expect_null(thm$j)
    expect_equal(W, thm$calc())
    expect_equal(as.integer(1:nrow(Xf)), thm$i)
    expect_equal(thm$i, thm$j)

})


subset_mask <- c(rep(T, 20), rep(F, nrow(X) - 20))
subset <- which(subset_mask)
X_t <- X[subset, ]
W_t <- X %*% solve(t(X_t) %*% diag(Phi_inv[subset]) %*% X_t + check_sig) %*% t(X_t) %*% diag(Phi_inv[subset])





test_that("thm$calc()", {
    # scalar inputs for i or j
    expect_equal(W[3, 4, drop = F], (thm$calc(3, 4)))
    expect_equal(W[3, , drop = F], (thm$calc(3)))
    expect_equal(W[, 3, drop = F], (thm$calc(j = 3)))

    # different arguments for subset
    expect_equal(W[3, 4, drop = F], (thm$calc(3, 4, subset = rep(TRUE, nrow(Xf)))))
    expect_equal(W[3, 4, drop = F], (thm$calc(3, 4, subset = 1:nrow(Xf))))

    # vector input for i, j
    i <- 1:5; j <- 6:10;
    check_w <- W[i, j]
    expect_equal(check_w, thm$calc(i, j))
    expect_equal(check_w,  hatmatrix(Xf, Phi_inv, Xr, Sigma_inv, i = i, j = j))
    expect_equal(W_t[, j, drop = F], thm$calc(subset = subset, j = j))
    expect_equal(W_t[, j, drop = F], thm$calc(subset = subset_mask, j = j))

})

test_that("wrapper functions around thm$calc() are consistent", {
    expect_equal(hatmatrix(Xf, Phi_inv, Xr, Sigma_inv), thm$calc())

})
