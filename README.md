# Poisson Fixed Effects Robust
Evan Wright

The package contains a function which computes standard errors
following Wooldridge (1999) for Poisson regression with
fixed effects, and a hypothesis test of the conditional mean
assumption (3.1). The standard errors are robust to within-group
autocorrelation of errors.

Example usage:
```
data("ex.dt.good")
pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), group.name = "id",
               data = ex.dt.good)
```

The standard errors should match those of the Stata command
`xtpoisson y x1 x2, fe vce(robust)`

Wooldridge, Jeffrey M. (1999): "Distribution-free estimation of some nonlinear
    panel data models," _Journal of Econometrics_, 90, 77-97.