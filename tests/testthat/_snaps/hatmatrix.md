# Internal consistency of HatMatrix

    Code
      thm
    Output
      <HatMatrix>
        Public:
          Phi_inv: active binding
          Sigma_inv: active binding
          V: active binding
          W: active binding
          Xf: active binding
          Xr: active binding
          calc: function (i = 1:self$n, j = 1:self$n, return_V = FALSE, subset = NULL) 
          clone: function (deep = FALSE) 
          diag_phi: function (idx) 
          expanded_sigma: function () 
          i: active binding
          initialize: function (Xf, Phi_inv, Xr = NULL, Sigma_inv = NULL) 
          j: active binding
          n: active binding
        Private:
          .Phi_inv: 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2  ...
          .Sigma_inv: 3 0 0 0 3 0 0 0 3
          .V: NA
          .W: NA
          .Xf: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  ...
          .Xr: 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1  ...
          .i: NULL
          .j: NULL
          .n: 32
          .save_values: function (Xf, Xr, Sigma_inv, Phi_inv) 
          .subset: NULL
          .subset_idx: function (subset) 

