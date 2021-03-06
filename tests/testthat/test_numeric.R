context("Numeric outputs")
library(poisFErobust)

test_that("Numeric output is accurate for good data", {
  out.list <- pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), 
                             group.name = "id", allow.set.key = TRUE, 
                             data = poisFErobust::ex.dt.good)
  expect_equivalent(out.list[[1]], c(0.9899730, 0.9917526), 
                    tolerance = 0.000001)
  expect_equivalent(out.list[[2]], c(0.03112512, 0.02481941), 
                    tolerance = 0.000001)
  expect_equivalent(out.list[[3]], 0.6996001, 
                    tolerance = 0.000001)
})

test_that("Numeric output is accurate for bad data", {
  out.list <- pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), 
                             group.name = "id", allow.set.key = TRUE,
                             data = poisFErobust::ex.dt.bad)
  expect_equivalent(out.list[[1]], c(0.4800735, 2.9866911), 
                    tolerance = 0.000001)
  expect_equivalent(out.list[[2]], c(0.2864666, 1.2743953), 
                    tolerance = 0.000001)
  expect_equivalent(out.list[[3]], 0.02213269, 
                    tolerance = 0.000001)
})

test_that("NA observations or non-factor group.name should raise error", {
  datawithNA <- copy(poisFErobust::ex.dt.good)
  datawithNA[1,1] = NA
  expect_error(pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), 
                              group.name = "id", allow.set.key = TRUE,
                              data = datawithNA),
               "NA detected.*")
  datawithnonfactor <- copy(poisFErobust::ex.dt.good)
  datawithnonfactor$id <- as.integer(datawithnonfactor$id)
  expect_error(pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), 
                              group.name = "id", allow.set.key = TRUE,
                              data = datawithnonfactor),
               "group.*")
})