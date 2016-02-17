library(arjtools)
context("Correct numerical results for WAAI and CAAI")

seg <- data.frame(
  start = c(526235, 3615271, 3616449, 25721244, 25727848, 33119458, 37434182, 37490345, 54896943, 56048173),
  end = c(3614856, 3616415, 25721029, 25721247, 33117320, 37432452, 37488368, 54895597, 56042460, 79916906),
  values = c(0.4816, -0.9807, 0.4455, -1.2194, 0.2684, -0.2726, 0.286, -0.3156, 0.2485, -0.3172),
  nprobes = c(872, 2, 13162, 2, 4418, 2486, 46, 9279, 971, 14513)
)

test_that("WAAI computes correct result", {
  expect_equal(waai(seg$values, seg$nprobes), -0.01071436471)
})

test_that("CAAI computes correct result", {
  expect_equal(caai(seg$values, seg$start, seg$end), 0.4568845879)
})
